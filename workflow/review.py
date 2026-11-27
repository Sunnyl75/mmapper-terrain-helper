"""Render map and regression views and check paired GIS trial exports."""
import argparse
import html
import json
from pathlib import Path
from collections import defaultdict

import numpy as np
from osgeo import gdal,ogr
from PIL import Image, ImageDraw
from scipy.ndimage import label
from scipy.spatial import cKDTree

from bake_profile import Pull, baseline_grid
from rebuild import CODE, COLORS, CROSS, WET, Fit, layer_raster, target_mask
from room_io import read_map
from mm2_background import verify_mm2
from river_reaches import Elevation, layer_paths, side_intersection

VIEWS = {
    "Old Ford and the Anduin": (468,-83,506,-50),
    "Ford of Bruinen": (365,-91,387,-70),
    "Lake Evendim": (143,-39,184,-6),
    "Long Lake": (621,-43,645,-7),
    "Dunland fords": (294,-222,327,-197),
    "Mitheithel and lower Bruinen": (245,-173,350,-101),
}
CHECKS = {
    "Dunland fords": [(311,-203),(311,-213)],
    "Ford of Bruinen to Tharbad": [(376,-80),(249,-171)],
}
NAMED = {"Bruinen": ("bruinen","loudwater"), "Mitheithel": ("hoarwell","mitheithel"),
         "Baranduin": ("baranduin","brandywine"), "Glanduin": ("glanduin",),
         "Isen": ("isen",), "Anduin": ("anduin",), "Celebrant": ("celebrant","silverlode")}


def reach_review(root, identity, survey, before_dem, before_gis):
    selected = json.loads((survey/"reach-fit.json").read_text())
    before=Elevation(before_dem)
    after=Elevation(root/"elevation.tif")
    old=layer_paths(before_gis,bounds=identity.bounds)
    new=layer_paths(root/"aligned-gis.gpkg",bounds=identity.bounds)
    paths=defaultdict(lambda:[[],[]])
    for group,features in enumerate((old,new)):
        for p in features:
            paths[p.fid][group].append(p)
    measures=defaultdict(lambda:defaultdict(list))
    missing=defaultdict(int)
    for row in selected["selected_stations"]:
        q=np.array(row["xy"])
        n=np.array(row["normal"])
        x,side=identity.profile.forward(np.stack([q,q+n*.4]))
        n=(side-x)/np.linalg.norm(side-x)
        group=row["source_id"]
        if row["evidence"]=="Old Ford survey":
            group="Old Ford"
        picks=[]
        for ds,traces,hint in zip((before,after),paths[row["source_id"]],(row["trial_trough"],0.)):
            trough=ds.measure(x,n,hint,2.8)
            nearest=min(((*p.tree.query(x),i) for i,p in enumerate(traces)),default=None)
            gis=side_intersection(traces[nearest[2]],x,n,nearest[1]) if nearest and nearest[0]<6 else None
            picks.append((trough,gis))
        for i,kind in enumerate(("dem_trough","gis_centreline")):
            a,b=picks[0][i],picks[1][i]
            if a is None or b is None:
                missing[str(group)+":"+kind]+=1
                continue
            measures[str(group)][kind].append((abs(a),abs(b)))
    return {key:{kind:{"paired_stations":len(vals),
                      "mean_trial_rooms":float(np.mean(vals,axis=0)[0]),
                      "mean_reach_fit_rooms":float(np.mean(vals,axis=0)[1]),
                      "p90_trial_rooms":float(np.quantile(vals,.9,axis=0)[0]),
                      "p90_reach_fit_rooms":float(np.quantile(vals,.9,axis=0)[1]),
                      "stations_without_a_paired_measurement":missing[key+":"+kind]}
                 for kind,vals in values.items()} for key,values in measures.items()}


def map_image(grid, bounds, shade, size):
    rgb=Image.fromarray(COLORS[grid[::-1]],"RGB").resize(shade.size,Image.Resampling.NEAREST)
    return Image.alpha_composite(rgb.convert("RGBA"),shade).convert("RGB").resize(size,Image.Resampling.LANCZOS)


def coordinate_crop(grid, bounds, shade, box, size):
    west,south,east,north=bounds
    w,s,e,n=box
    density=shade.width/(east-west)
    return map_image(grid[s-south:n-south,w-west:e-west],box,
                     shade.crop(tuple(round(i*density) for i in (w-west,north-n,e-west,north-s))),size)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("target","profile-fit","output"):
        p.add_argument("--"+arg,type=Path,required=True)
    p.add_argument("--previous-background",type=Path)
    p.add_argument("--previous-hillshade",type=Path)
    p.add_argument("--previous-dem",type=Path)
    p.add_argument("--previous-gis",type=Path)
    p.add_argument("--reach-fit",type=Path)
    p.add_argument("--exclude-outdoor",type=Path)
    a=p.parse_args()
    import csv
    hidden=set()
    if a.exclude_outdoor:
        with a.exclude_outdoor.open(newline="") as f:
            hidden={(int(row["x"]),int(row["y"])) for row in csv.DictReader(f)}
    out=a.output
    report=json.loads((out/"build-report.json").read_text())
    sector_file=np.load(out/"sectors-review.npz")
    profile=np.load(a.profile_fit)
    f=Fit(np.load(out/"reference-pull.npz")["offset"],
          Pull(profile["offset"],profile["bounds"],profile["center"]))
    bounds=tuple(sector_file["bounds"])
    west,south,east,north=bounds
    sector=sector_file["sector"]
    dem=gdal.Open(str(out/"elevation.tif"))
    vectors=ogr.Open(str(out/"aligned-gis.gpkg"))
    width,height=east-west,north-south
    name=report["profile"]
    shade=Image.open(out/f"{name}.hillshade.png").convert("RGBA")
    step,_,_,dy,lx,ly=map(float,(out/f"{name}.hillshade.pgw").read_text().split())
    dem_w,dd,_,dem_n,_,ddd=dem.GetGeoTransform()
    if not np.allclose([step*shade.width,-dy*shade.height,lx-step/2,ly-dy/2,
                        dem_w,dem_n,dem.RasterXSize*dd,-dem.RasterYSize*ddd],
                       [width,height,west,north,west,north,width,height]):
        raise ValueError("XML, hillshade and DEM must share the same room bounds")
    rooms=read_map(a.target)
    binary=out/f"{name}-background-river-reaches.mm2"
    roundtrip, road_flags=baseline_grid(read_map(binary),bounds)
    if not np.array_equal(roundtrip,sector):
        raise ValueError("MM2 sector grid differs from the rendered terrain")
    known,search,known_flags=target_mask(rooms,bounds,hidden)
    if not np.array_equal(sector[search],known[search]) or not np.array_equal(road_flags[search],known_flags[search]):
        raise ValueError("An eligible MUME room or main road was overwritten")
    report["mm2_verification"]=verify_mm2(binary,sector,road_flags,bounds)
    if a.previous_dem and a.previous_gis and a.reach_fit:
        report["river_reach_station_audit"]=reach_review(out,f,a.reach_fit,a.previous_dem,a.previous_gis)
    initial=np.array([(x,y) for y in range(south+1,north,6) for x in range(west+1,east,6)])
    err=np.linalg.norm(initial-f.target(f.native(initial)),axis=1)
    report["gis_dem_round_trip_error_rooms"]={"median":float(np.median(err)),"p99":float(np.quantile(err,.99)),"max":float(err.max())}
    if np.quantile(err,.99)>.12:
        raise ValueError("The GIS forward transform disagrees with the elevation inverse transform")
    labels,_=label(np.isin(sector,WET),CROSS)
    for title,coords in CHECKS.items():
        axy=f.profile.forward(np.array(coords)+.5)
        target_coords=np.floor(axy).astype(int)
        marks=[int(labels[y-south,x-west]) for x,y in target_coords]
        report.setdefault("regression_locations",{})[title]={
            "target_room_checkpoints":target_coords.tolist(),"wet_labels":marks,
            "same_water_component":len(set(marks))==1 and 0 not in marks}
    rivers=vectors.GetLayerByName("Rivers")
    long_lake_points=[]
    for fid in (1957,907,2287):
        rivers.SetAttributeFilter(f"source_id = {fid}")
        feature=rivers.GetNextFeature()
        points=np.asarray(feature.GetGeometryRef().GetGeometryRef(0).GetPoints())[:,:2]
        if fid==2287:
            y=f.profile.forward(np.array([[637.5,-41.]]))[0,1]
            point=points[np.argmin(abs(points[:,1]-y))]
        else:
            point=points[-1]
        long_lake_points.append(np.floor(point).astype(int).tolist())
    marks=[int(labels[y-south,x-west]) for x,y in long_lake_points]
    report["regression_locations"]["Long Lake inlet to outlet"]={
        "target_room_checkpoints":long_lake_points,"wet_labels":marks,
        "same_water_component":len(set(marks))==1 and 0 not in marks}
    checks={}
    for feature,names in NAMED.items():
        rivers.SetAttributeFilter(f"name = '{feature}'")
        paths=[np.asarray(g.GetPoints())[:,:2] for item in rivers for g in item.GetGeometryRef()]
        if not paths:
            continue
        xyz=np.vstack(paths)
        considered=[r.xyz[:2] for r in rooms if r.outdoor and CODE[r.terrain] in WET and
                     any(a in r.name.casefold() for a in names)]
        if not considered:
            continue
        distance=cKDTree(xyz).query(np.array(considered)+.5)[0]
        checks[feature]={"named_foreground_water_rooms":len(distance),
                         "median_rooms_to_warped_gis":float(np.median(distance)),
                         "p90_rooms_to_warped_gis":float(np.quantile(distance,.9))}
    report["named_surface_water_alignment"]=checks
    out.joinpath("build-report.json").write_text(json.dumps(report,indent=2)+"\n")
    preview=map_image(sector,bounds,shade,(1454,584))
    preview.save(out/"overview.png",optimize=True)
    before=old_shade=None
    if a.previous_background and a.previous_hillshade:
        before,_=baseline_grid(read_map(a.previous_background),bounds)
        old_shade=Image.open(a.previous_hillshade).convert("RGBA")
    lake_results={}
    for lake_name in ("Lake Evendim","Long Lake"):
        layer=vectors.GetLayerByName("lakes")
        layer.SetAttributeFilter(f"name = '{lake_name}'")
        feature=layer.GetNextFeature()
        xmin,xmax,ymin,ymax=feature.GetGeometryRef().GetEnvelope()
        w,s=max(west,int(np.floor(xmin))-2),max(south,int(np.floor(ymin))-2)
        e,n=min(east,int(np.ceil(xmax))+2),min(north,int(np.ceil(ymax))+2)
        # Count subpixels within the GIS shoreline that the room export
        # displays as wet, and separately check the flat DEM basin core.
        rate=round(1/dd)
        footprint=layer_raster(layer,(w,s,e,n),rate).astype(bool)
        z=dem.ReadAsArray((w-west)*rate,(north-n)*rate,(e-w)*rate,(n-s)*rate)[::-1]
        heights,counts=np.unique(z[footprint],return_counts=True)
        basin_labels,_=label(abs(z.astype(np.int32)-int(heights[counts.argmax()]))<=1)
        overlap,counts=np.unique(basin_labels[footprint&(basin_labels>0)],return_counts=True)
        core=basin_labels==overlap[counts.argmax()]
        lake_results[lake_name]={}
        for image,key in ((before,"previous"),(sector,"gis_first")):
            if image is None:
                continue
            water=np.isin(image[s-south:n-south,w-west:e-west],WET)
            mask=np.repeat(np.repeat(water,rate,axis=0),rate,axis=1)
            lake_results[lake_name][key]={
                "percent_of_GIS_lake_in_wet_cells":round(float(mask[footprint].mean())*100,2),
                "percent_of_DEM_basin_core_in_wet_cells":round(float(mask[core].mean())*100,2)}
    report["shoreline_review"]=lake_results
    out.joinpath("build-report.json").write_text(json.dumps(report,indent=2)+"\n")
    sections=[]
    for i,(title,box) in enumerate(VIEWS.items(),1):
        w,s,e,n=box
        p1,p2=f.profile.forward(np.array([[w,s],[e,n]]))
        box=tuple(map(int,np.rint([*p1,*p2])))
        w,s,e,n=box
        size=(520,max(190,min(540,int(520*(n-s)/(e-w)))))
        now=coordinate_crop(sector,bounds,shade,box,size)
        if before is not None:
            old=coordinate_crop(before,bounds,old_shade,box,size)
            img=Image.new("RGB",(2*size[0]+16,size[1]+28),"#faf8f3")
            img.paste(old,(0,28))
            img.paste(now,(size[0]+16,28))
            draw=ImageDraw.Draw(img)
            draw.text((10,7),"GIS-FIRST TRIAL",fill="#252723")
            draw.text((size[0]+27,7),"RIVER-REACH FIT",fill="#252723")
        else:
            img=now
        filename=f"detail-{i}.png"
        img.save(out/filename,optimize=True)
        sections.append((title,filename))
    anchors=" · ".join(f'<a href="#detail-{i}">{html.escape(t)}</a>' for i,(t,_) in enumerate(sections,1))
    content="".join(f'<section id="detail-{i}"><h2>{html.escape(t)}</h2><img src="{p}" alt="{html.escape(t)}"></section>'
                    for i,(t,p) in enumerate(sections,1))
    rows="".join(f'<li>{html.escape(k)} · {v["gis_first"]["percent_of_DEM_basin_core_in_wet_cells"]:.1f}% of the DEM basin core in wet rooms</li>'
                 for k,v in lake_results.items())
    out.joinpath("review.html").write_text(f'''<!doctype html><html lang="en"><meta charset="utf-8">
<title>{name.title()} terrain alignment</title>
<style>body{{max-width:1100px;margin:2rem auto;padding:0 1rem;font:16px system-ui;color:#292d25;background:#f6f4ed}}img{{max-width:100%}}section{{margin:2.5rem 0}}a{{color:#165242}}</style>
<h1>{name.title()} terrain alignment</h1><p>{anchors} · <a href="build-report.json">Build report</a></p>
<ul>{rows}</ul>
<img src="overview.png" alt="Full map with GIS sectors and hillshade">{content}</html>''',encoding="utf-8")
    print(name,report["regression_locations"],report["gis_dem_round_trip_error_rooms"])


if __name__=="__main__":
    main()
