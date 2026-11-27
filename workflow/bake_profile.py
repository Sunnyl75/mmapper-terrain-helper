"""Build matching MMapper background, elevation, and hillshade for a specific map snapshot."""
import argparse
from collections import Counter, defaultdict, deque
import csv
from hashlib import file_digest
import json
from pathlib import Path
import tempfile

import numpy as np
from osgeo import gdal
from scipy.ndimage import distance_transform_edt, map_coordinates

from bake_hillshade import bake
from register import Fit, blocked_validation, choose, match, robust_controls, sampled_field, summary
from room_io import CARDINAL, OUTDOOR, STEPS, TERRAIN, read_map

WET={TERRAIN.index(v) for v in ("WATER","RAPIDS","SHALLOW")}
ROAD=TERRAIN.index("ROAD")
RATE=2


def digest(path):
    with path.open("rb") as stream:
        return file_digest(stream,"sha256").hexdigest()


class Pull:
    def __init__(self,field,bounds,center):
        self.field=field
        self.bounds=bounds
        self.center=center

    def delta(self,xy):
        if not self.field.any():
            return np.broadcast_to(self.center,xy.shape).copy()
        west,south,_,_=self.bounds
        indices=[(xy[:,1]-south)*RATE,(xy[:,0]-west)*RATE]
        return self.center+np.column_stack([
            map_coordinates(self.field[...,i],indices,order=1,mode="nearest") for i in range(2)])

    def forward(self,original):
        # Damp inverse iteration for crossings of regional seams.
        guess=np.asarray(original,dtype=float)-self.center
        for _ in range(20):
            residual=original-guess-self.delta(guess)
            guess+=.65*residual
            if np.max(abs(residual),initial=0)<.025:
                break
        return guess


def working_dem(path):
    src=gdal.Open(str(path))
    west,dx,rx,north,ry,dy=src.GetGeoTransform()
    if dx<=0 or dy!=-dx or rx or ry:
        raise ValueError("Canonical raster must have square, north-up pixels")
    east=west+src.RasterXSize*dx
    south=north+src.RasterYSize*dy
    bounds=np.rint([west,south,east,north]).astype(int)
    if not np.allclose(bounds,[west,south,east,north]) or src.GetRasterBand(1).DataType!=gdal.GDT_UInt16:
        raise ValueError("Canonical raster must be flat-grid UInt16 with whole-room bounds")
    return src,bounds,dx


def resample_dem(src,path,dx,pull):
    west,_,_,north=pull.bounds
    out=gdal.GetDriverByName("GTiff").Create(str(path),src.RasterXSize,src.RasterYSize,1,gdal.GDT_UInt16,
        options=["TILED=YES","COMPRESS=DEFLATE","PREDICTOR=2","NUM_THREADS=2"])
    out.SetProjection(src.GetProjection())
    out.SetGeoTransform((int(west),dx,0,int(north),0,-dx))
    out.SetMetadata(src.GetMetadata())
    out.SetMetadataItem("ADAPTATION","River-refined Arda elevation reprojected to this MMapper map")
    target=out.GetRasterBand(1)
    raw=src.ReadAsArray()
    if not pull.field.any():
        target.WriteArray(raw)
        out=None
        return
    height,width=raw.shape
    source_w,_,_,source_n,_,_=src.GetGeoTransform()
    x=west+(np.arange(width)+.5)*dx
    cols=np.arange(width)
    for top in range(0,height,64):
        rows=np.arange(top,min(top+64,height))
        y=north-(rows+.5)*dx
        xx,yy=np.meshgrid(x,y)
        points=np.column_stack((xx.ravel(),yy.ravel()))
        delta=pull.delta(points)-pull.center
        touched=np.any(delta!=0,axis=1).reshape(xx.shape)
        changed=raw[rows].copy()
        iy,ix=np.nonzero(touched)
        d=delta.reshape((*xx.shape,2))[iy,ix]/dx
        if len(iy):
            cy=rows[iy]-d[:,1]
            cx=cols[ix]+d[:,0]
            if cy.min()<-.5 or cy.max()>height-.5 or cx.min()<-.5 or cx.max()>width-.5:
                raise ValueError("Regional correction exceeds source elevation coverage")
            changed[iy,ix]=map_coordinates(raw,[cy,cx],order=1,mode="nearest")
        target.WriteArray(changed,0,int(top))
    out=None


def baseline_grid(rooms,bounds):
    west,south,east,north=bounds
    shape=(north-south,east-west)
    base=np.zeros(shape,np.uint8)
    road=np.zeros(shape,np.uint8)
    for r in rooms:
        x,y,z=r.xyz
        if z or r.terrain not in OUTDOOR or not (west<=x<east and south<=y<north):
            raise ValueError("Reference background must contain only its aligned outdoor room grid")
        if base[y-south,x-west]:
            raise ValueError("Duplicate cell in reference background")
        base[y-south,x-west]=TERRAIN.index(r.terrain)
        road[y-south,x-west]=r.flags if r.terrain=="ROAD" else 0
    if not base.all():
        raise ValueError("Reference background does not fill the elevation footprint")
    return base,road


def neighbors(cell,shape):
    x,y=cell
    return [(x+dx,y+dy) for dx,dy in STEPS if 0<=x+dx<shape[1] and 0<=y+dy<shape[0]]


def short_route(image,blocked,start,end,allowed):
    # Search only near the span of this inherited edge. Existing network
    # cells may be visited; missing cells may be filled only when inferred.
    bound=min(6,abs(start[0]-end[0])+abs(start[1]-end[1])+2)
    states=deque([(start,())])
    seen={start}
    while states:
        p,route=states.popleft()
        if p==end:
            return route
        if len(route)>=bound:
            continue
        choices=neighbors(p,image.shape)
        choices.sort(key=lambda q:image[q[1],q[0]] not in allowed)
        for q in choices:
            if q in seen:
                continue
            tile=image[q[1],q[0]]
            if tile not in allowed and (blocked[q[1],q[0]] or tile in WET or tile==ROAD):
                continue
            seen.add(q)
            states.append((q,route+(q,)))
    return None


def repair_network(source,image,blocked,pull,source_bounds):
    west,south,_,_=source_bounds
    fx,fy,_,_=pull.bounds
    counters={}
    for kind,codes,fill in (("water",WET,TERRAIN.index("WATER")),("road",{ROAD},ROAD)):
        mask=np.isin(source,list(codes))
        yy,xx=np.nonzero(mask)
        projected=np.floor(pull.forward(np.column_stack((xx+west+.5,yy+south+.5)))).astype(int)
        mapped={tuple(a): tuple(b) for a,b in zip(np.column_stack((xx,yy)),projected-[fx,fy])
                if 0<=b[0]-fx<image.shape[1] and 0<=b[1]-fy<image.shape[0]}
        repaired=unresolved=skipped=0
        for old,start in mapped.items():
            for dx,dy in ((1,0),(0,1)):
                end=mapped.get((old[0]+dx,old[1]+dy))
                if end is None or start==end:
                    continue
                if image[start[1],start[0]] not in codes or image[end[1],end[0]] not in codes:
                    skipped+=1
                    continue
                length=abs(start[0]-end[0])+abs(start[1]-end[1])
                if length<=1:
                    continue
                route=short_route(image,blocked,start,end,codes) if length<=4 else None
                if route is None:
                    unresolved+=1
                    continue
                for x,y in route:
                    if image[y,x] not in codes:
                        image[y,x]=fill
                        repaired+=1
        counters[kind]={"inferred_cells_added_for_short_gaps":repaired,
                        "separated_inherited_edges_not_bridged":unresolved,
                        "inherited_edges_with_reclassified_endpoints":skipped}
    return counters


def background(src_grid,src_flags,source_bounds,target,pull):
    west,south,east,north=pull.bounds
    sy,sx=np.indices((north-south,east-west))
    world=np.column_stack((sx.ravel()+west+.5,sy.ravel()+south+.5))
    xy=np.floor(world+pull.delta(world)).astype(int)-source_bounds[:2]
    if np.any(xy.min(axis=0)<0) or np.any(xy.max(axis=0)>=np.array(src_grid.shape[::-1])):
        raise ValueError("A background cell lies outside canonical coverage")
    image=src_grid[xy[:,1],xy[:,0]].reshape(sy.shape).copy()
    flags=src_flags[xy[:,1],xy[:,0]].reshape(sy.shape).copy()
    blocked=np.zeros(sy.shape,bool)
    known=np.zeros(sy.shape,np.uint8)
    known_flags=np.zeros(sy.shape,np.uint8)
    excluded=np.zeros(sy.shape,bool)
    outside=[]
    seen=set()
    for room in target:
        x,y,z=room.xyz
        if z:
            continue
        if not (west<=x<east and south<=y<north):
            if room.outdoor:
                outside.append({"xy":[x,y],"room":room.name,"server_id":room.server_id})
            continue
        # Mappers can contain overlaps. Require review rather than choosing a room by file order.
        if (x,y) in seen:
            raise ValueError(f"Duplicate target surface room at {(x,y)}")
        seen.add((x,y))
        q=(y-south,x-west)
        blocked[q]=True
        if room.outdoor:
            image[q]=known[q]=TERRAIN.index(room.terrain)
            known_flags[q]=room.flags if room.terrain=="ROAD" else 0
        else:
            excluded[q]=True
    to_fill=excluded & np.isin(image,list(WET|{ROAD}))
    if to_fill.any():
        _,nearest=distance_transform_edt(excluded | np.isin(image,list(WET|{ROAD})),return_indices=True)
        image[to_fill]=image[tuple(nearest[:,to_fill])]

    repairs=repair_network(src_grid,image,blocked,pull,source_bounds) if pull.field.any() else {}
    is_road=image==ROAD
    for k,(dx,dy) in enumerate(STEPS):
        neighbor=np.roll(is_road,(-dy,-dx),axis=(0,1))
        if dx==1: neighbor[:,-1]=False
        elif dx==-1: neighbor[:,0]=False
        if dy==1: neighbor[-1,:]=False
        elif dy==-1: neighbor[0,:]=False
        flags[is_road & ~blocked & neighbor]|=1<<k
    flags[blocked]=known_flags[blocked]
    flags[~is_road]=0
    return image,flags,known,known_flags,outside,repairs


def write_background(path,grid,flags,bounds):
    width=grid.shape[1]
    west,south,east,north=bounds
    with path.open("w",encoding="utf-8") as f:
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n<map type="mmapper2xml" version="1.0.0">\n')
        for y in range(north-1,south-1,-1):
            for x in range(west,east):
                ix,iy=x-west,y-south
                token=int(grid[iy,ix])
                key=(north-1-y)*width+ix
                fields=[f'<room id="{key}"><terrain>{TERRAIN[token]}</terrain><coord x="{x}" y="{y}" z="0"/>']
                for i,((dx,dy),direction) in enumerate(zip(STEPS,CARDINAL)):
                    fields.append(f'<exit dir="{direction}">')
                    if west<=x+dx<east and south<=y+dy<north:
                        fields.append(f'<to>{key+dx-dy*width}</to>')
                    if flags[iy,ix]&(1<<i):
                        fields.append("<exitflag>ROAD</exitflag>")
                    fields.append("</exit>")
                f.write("".join(fields)+"</room>\n")
        f.write("</map>\n")


def build(a, paths):
    reference,target=read_map(a.reference),read_map(a.target)
    if a.bounds is None:
        src,source_bounds,dx=working_dem(a.dem)
    else:
        source_bounds=np.array(a.bounds,dtype=int)
    pairs,topology=match(reference,target,a.controls)
    if len(pairs)<30:
        raise ValueError("Need more matched outdoor controls in the target geometry")
    left,right,exclude=robust_controls(pairs,topology)
    # Strong, compact detached map areas (for example a standalone mini-map)
    # provide game cells but aren't useful controls for the regional surface.
    median=Fit(right[~exclude],left[~exclude],"translation").center
    detached=(np.linalg.norm(left-right-median,axis=1)>a.max_regional_shift)&np.array([p.evidence!="manual" for p in pairs])
    eligible=~(exclude|detached)
    scores,errors,fold=blocked_validation(right[eligible],left[eligible])
    selection=choose(scores)
    fit=Fit(right[eligible],left[eligible],selection)
    center=fit.center.astype(int)
    bounds=source_bounds-np.tile(center,2)
    field,checks=sampled_field(fit,bounds,RATE) if selection!="translation" else (
        np.zeros(((bounds[3]-bounds[1])*RATE+1,(bounds[2]-bounds[0])*RATE+1,2),dtype=np.float32),
        {"minimum_inverse_jacobian":1.,"minimum_local_regional_scale":1.,"field_grid_samples_damped":0})
    pull=Pull(field,bounds,center)
    low=summary(np.linalg.norm(right[eligible]+pull.delta(right[eligible])-left[eligible],axis=1))
    report={"profile":a.profile,"reference":a.reference.name,"target":a.target.name,
            "canonical_bounds":source_bounds.tolist(),"target_bounds":bounds.tolist(),
            "profile_review":"Regions outside shared evidence or with high residuals require visual review.",
            "source_files":{str(p.name):digest(p) for p in (a.reference,a.target,a.dem,a.background) if p},
            "matched_control_count":len(pairs),"eligible_geographic_controls":int(eligible.sum()),
            "match_sources":dict(Counter(p.evidence for p,e in zip(pairs,eligible) if e)),
            "rejected_isolated_pairs":int(exclude.sum()),"detached_geometry_pairs":int((detached & ~exclude).sum()),
            "regional_shift_limit_rooms":a.max_regional_shift,"selected_model":selection,
            "target_to_reference_global_pull_xy":center.tolist(),"spatially_held_out_scores":scores,
            "final_damped_field_fitted_residual":low,"deformation":checks}
    if a.controls:
        report["manual_control_sha256"]=digest(a.controls)
    regions=defaultdict(list)
    values=np.full(len(pairs),np.nan)
    values[eligible]=errors[selection]
    move_residual=np.linalg.norm(right+pull.delta(right)-left,axis=1)
    for xy,v in zip(right[eligible],values[eligible]):
        regions[tuple(np.floor(xy/32).astype(int)*32)].append(v)
    report["regions_32_rooms"]=[{"xy":list(map(int,k)),**summary(np.array(v))}
                                  for k,v in sorted(regions.items())]
    with paths["controls"].open("w",newline="",encoding="utf-8") as f:
        writer=csv.writer(f)
        writer.writerow(("reference_x","reference_y","target_x","target_y","evidence","shared_cardinal_neighbors",
                         "room_name","status","held_out_error_rooms","fitted_field_error_rooms"))
        for i,p in enumerate(pairs):
            status="used" if eligible[i] else "isolated" if exclude[i] else "detached"
            writer.writerow((*p.reference,*p.target,p.evidence,int(topology[i]),p.name,status,
                             round(values[i],4) if eligible[i] else "",round(move_residual[i],4)))
    np.savez_compressed(paths["field"],offset=field,bounds=bounds,center=center,rate=RATE)
    if not a.audit_only:
        grid,flags=baseline_grid(read_map(a.background),source_bounds)
        image,roads,known,known_flags,outside,repairs=background(grid,flags,source_bounds,target,pull)
        if not (np.array_equal(image[known>0],known[known>0]) and
                np.array_equal(roads[known==ROAD],known_flags[known==ROAD])):
            raise ValueError("Generated grid changed a mapped outdoor sector or road exit")
        report["outside_canonical_coverage"]=outside
        report["adjustments_to_inferred_connectivity"]=repairs
        report["outdoor_target_sectors_preserved"]=int((known>0).sum())
        report["target_road_cells_preserved"]=int((known==ROAD).sum())
        write_background(paths["background"],image,roads,bounds)
        resample_dem(src,paths["dem"],dx,pull)
        bake(paths["dem"],paths["shade"],a.pixels_per_room)
        assets=(paths["background"],paths["dem"],paths["shade"],paths["shade"].with_suffix(".pgw"))
        report["output_sha256"]={asset.name:digest(asset) for asset in assets}
    paths["report"].write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--profile",required=True,choices=("bundled","faine","nazgum"))
    for key in ("reference","target","output"):
        p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--dem",type=Path)
    p.add_argument("--background",type=Path)
    p.add_argument("--bounds",nargs=4,type=int,metavar=("WEST","SOUTH","EAST","NORTH"))
    p.add_argument("--controls",type=Path)
    p.add_argument("--pixels-per-room",type=int,default=11)
    p.add_argument("--max-regional-shift",type=float,default=10)
    p.add_argument("--audit-only",action="store_true")
    a=p.parse_args()
    if a.bounds and (not a.audit_only or a.bounds[0]>=a.bounds[2] or a.bounds[1]>=a.bounds[3]):
        p.error("--bounds requires --audit-only and ordered whole-room edges")
    if (a.dem is None and a.bounds is None) or (not a.audit_only and a.background is None):
        p.error("Supply --dem and --background or use --audit-only with --bounds")
    if a.output.exists():
        p.error("Choose a fresh output folder for this map snapshot")
    if a.pixels_per_room<=0 or a.max_regional_shift<=0:
        p.error("Pixel density and maximum regional shift must be positive")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    gdal.UseExceptions()
    with tempfile.TemporaryDirectory(prefix=".profile-build-",dir=a.output.parent) as tmp:
        out=Path(tmp)/"profile"
        out.mkdir()
        prefix=a.profile
        paths={name:out/filename for name,filename in {
            "background":prefix+"-background.mm2xml","dem":prefix+"-elevation.tif",
            "shade":prefix+".hillshade.png","controls":"matched-controls.csv",
            "field":"regional-displacement.npz","report":"profile-report.json"}.items()}
        data=build(a,paths)
        out.rename(a.output)
    print(json.dumps({k:data[k] for k in
          ("profile","eligible_geographic_controls","selected_model","target_bounds","deformation")},indent=2))


if __name__=="__main__":
    main()
