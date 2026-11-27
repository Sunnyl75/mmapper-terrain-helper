"""Build a paired MMapper XML, DEM and hillshade from native Arda GIS and elevation."""
import argparse
import csv
from collections import Counter
import heapq
import json
from pathlib import Path

import numpy as np
from osgeo import gdal, ogr, osr
from scipy.ndimage import (distance_transform_edt, gaussian_filter, label, map_coordinates,
                           maximum_filter, minimum_filter)
from scipy.spatial import cKDTree

from bake_hillshade import bake
from bake_profile import Pull, digest
from fit_reference import (BOUNDS, DESTINATION, ORIGIN, TO_ROOM, TO_NATIVE, geom_grid)
from mm2_background import verify_mm2, write_mm2
from room_io import STEPS, TERRAIN, read_map

RATE = 2
SUPERSAMPLE = 6
CROSS = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])
CODE = {name: i for i,name in enumerate(TERRAIN)}
WET = [CODE[name] for name in ("WATER", "SHALLOW", "RAPIDS")]
LAYERS = ("Coastline", "lakes", "Wetlands", "Rivers", "Roads", "forests", "Hills", "Mountains")
PINS = {
    "Long Lake north": (635.035, -9.4736),
    "Long Lake south": (637.4386, -32.5561),
    "Forest River mouth": (634.7682, -26.2981),
    "Evendim west shore": (148.1, -21.5),
    "Evendim east shore": (170.9, -16.9),
    "Evendim south shore": (160., -32.2),
    "Evendim north shore": (160.5, -10.),
    "Mitheithel–Bruinen GIS junction": (299.125,-126.304),
}
DUNLAND = [
    ("Dunland upper ford", (311.46, -206.12), (311.5, -202.5), 9.0),
    ("Dunland lower ford", (311.51, -213.31), (311.5, -212.5), 8.0),
]
COLORS = np.array([
    (35,35,35), (36,36,36), (36,36,36), (130,143,83), (70,112,68),
    (151,141,89), (121,112,105), (73,143,146), (39,100,145), (52,126,162),
    (37,95,139), (191,159,109), (110,135,77), (33,33,33), (36,36,36)
], dtype=np.uint8)


def sample_field(field, xy, bounds):
    west,south,_,_=bounds
    inds = [(xy[:,1]-south)*RATE, (xy[:,0]-west)*RATE]
    return np.column_stack([map_coordinates(field[...,k],inds,order=1,mode="constant")
                            for k in range(2)])


def gauss_normals(a, b, radius):
    d2 = np.sum((a[:,None,:]-b[None,:,:])**2,axis=2)
    # Compact at twice the radius, with zero slope at the cutoff.
    q = np.clip(d2/(4*radius**2),0,1)
    return np.exp(-d2/(2*radius**2))*(1-q)**2


def make_reference_fit(control_path, joins_path, ref):
    west,south,east,north = BOUNDS
    w,h = (east-west)*RATE+1,(north-south)*RATE+1
    diff = np.zeros((h,w,2), np.float32)
    weight = np.zeros((h,w),np.float32)
    data = np.atleast_1d(np.genfromtxt(control_path,delimiter=",",names=True))
    target = np.column_stack([data["target_x"],data["target_y"]])
    source = np.column_stack([data["source_x"],data["source_y"]])
    displacement = source-target
    excluded = np.linalg.norm(displacement,axis=1)>7.5
    # Isolate bad correlations without forcing genuine regional offsets to zero.
    near = cKDTree(target).query(target,k=7)[1][:,1:]
    local = np.median(displacement[near],axis=1)
    excluded |= np.linalg.norm(displacement-local,axis=1)>2.3
    quality = np.clip((data["correlation"]-.89)/.10,0,1)
    quality *= np.clip(data["distinctiveness"]/.07,.3,1.)
    quality *= np.clip((18.-data["distance_to_mapped_surface"])/14,.25,1)
    xx, yy = np.meshgrid(west+np.arange(w)/RATE,south+np.arange(h)/RATE)
    fit_xy=np.column_stack([xx.ravel(),yy.ravel()])
    ix,iy = np.rint((target-[west,south])*RATE).astype(int).T
    holdout = ((target[:,0].astype(int)//16+3*(target[:,1].astype(int)//16))%7==0)&~excluded

    def interpolate(available):
        weights=weight.copy()
        delta=diff.copy()
        weights[iy[available],ix[available]]=quality[available]
        delta[iy[available],ix[available]]=displacement[available]*quality[available,None]
        sigma=6.*RATE
        count=gaussian_filter(weights,sigma,mode="constant")
        numerator=gaussian_filter(delta,(sigma,sigma,0),mode="constant")
        distance=distance_transform_edt(weights==0)/RATE
        t=np.clip((30-distance)/21,0,1)
        weighted=numerator/np.maximum(count[...,None],.00012)*(t*t*(3-2*t))[...,None]
        return gaussian_filter(weighted,(2.0*RATE,2.0*RATE,0),mode="constant")

    core=interpolate(~excluded)
    validation=interpolate(~excluded&~holdout)
    validation_before=np.linalg.norm(displacement[holdout],axis=1)
    validation_after=np.linalg.norm(displacement[holdout]-sample_field(validation,target[holdout],BOUNDS),axis=1)
    anchors=[]
    outdoor={r.xyz[:2] for r in ref if r.outdoor}
    with joins_path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if not row["bundled_room_x"]:
                continue
            room = int(row["bundled_room_x"]),int(row["bundled_room_y"])
            if room not in outdoor:
                raise ValueError(f"Reference map is missing control {row['join']}")
            dst=(float(row["matched_x"]),float(row["matched_y"]))
            src=(float(row["source_x"]),float(row["source_y"]))
            anchors.append((row["join"],src,dst,15.))
    anchors.extend(DUNLAND)
    anchors.extend((name,point,point,9.) for name,point in PINS.items())
    names,source_join,dest_join,radii=zip(*anchors)
    dest_join=np.array(dest_join)
    source_join=np.array(source_join)
    radii=np.array(radii)
    anchor_delta=source_join-dest_join
    correction=anchor_delta-sample_field(core,dest_join,BOUNDS)
    basis=gauss_normals(dest_join,dest_join,radii)
    coefficients=np.linalg.solve(basis+np.eye(len(anchors))*.025,correction)
    shift=np.zeros_like(core)
    for start in range(0,len(fit_xy),24000):
        stop=start+24000
        shift.reshape(-1,2)[start:stop]=gauss_normals(fit_xy[start:stop],dest_join,radii)@coefficients
    distance_edge=np.minimum.reduce([xx-west,east-xx,yy-south,north-yy])
    fade=np.clip(distance_edge/13,0,1)
    field=(core+shift)*(fade*fade*(3-2*fade))[...,None]
    dy,dx=np.gradient(field,1/RATE,1/RATE,axis=(0,1))
    jac=(1+dx[...,0])*(1+dy[...,1])-dx[...,1]*dy[...,0]
    if jac.min()<.18:
        raise ValueError(f"Reference fit has excessive deformation, determinant {jac.min():.3f} at {(fit_xy[jac.argmin()]).tolist()}")
    errors=np.linalg.norm(sample_field(field,dest_join,BOUNDS)-anchor_delta,axis=1)
    info={
        "archive_terrain_controls_retained":int((~excluded).sum()),
        "holdout_count":int(holdout.sum()),
        "holdout_median_rooms_similarity_fit":float(np.median(validation_before)),
        "holdout_median_rooms_with_local_fit":float(np.median(validation_after)),
        "holdout_p90_rooms_similarity_fit":float(np.quantile(validation_before,.9)),
        "holdout_p90_rooms_with_local_fit":float(np.quantile(validation_after,.9)),
        "inverse_jacobian_minimum":float(jac.min()),
        "inverse_jacobian_maximum":float(jac.max()),
        "maximum_local_shift_rooms":float(np.linalg.norm(field,axis=2).max()),
        "reviewed_controls":{name:{"source_base_room":src.tolist(),"reference_room":dst.tolist(),
                                 "residual_rooms":float(error)}
                             for name,src,dst,error in zip(names,source_join,dest_join,errors)}
    }
    return field.astype(np.float32),info


class Fit:
    def __init__(self, field, profile, base=None, surveyed=None):
        self.reference = Pull(field, BOUNDS,np.zeros(2))
        self.profile = profile
        self.bounds = tuple(map(int,profile.bounds))
        self.source_registration={}
        self.endpoints=None
        if base is not None:
            self.base=Pull(base,BOUNDS,np.zeros(2))
        if surveyed:
            grouped={}
            for row in surveyed:
                point=np.array(row["xy"])
                normal=np.array(row["normal"])
                centre=point+row["trial_gis"]*normal
                offset=np.clip(row["trial_trough"]-row["trial_gis"],-2.8,2.8)*normal
                grouped.setdefault(row["source_id"],[]).append((centre,offset))
            for fid,pairs in grouped.items():
                centres,offsets=map(np.array,zip(*pairs))
                self.source_registration[fid]=(cKDTree(centres),offsets)

    def native(self, xy):
        ref=xy+self.profile.delta(xy)
        src=ref+self.reference.delta(ref)
        return (src-DESTINATION) @ TO_NATIVE+ORIGIN

    def target(self, xy):
        source=(xy-ORIGIN) @ TO_ROOM+DESTINATION
        return self.profile.forward(self.reference.forward(source))

    def river_target(self, xy, fid):
        if fid not in self.source_registration:
            return self.target(xy)
        world=(xy-ORIGIN) @ TO_ROOM+DESTINATION
        trial=self.base.forward(world)
        stations,offsets=self.source_registration[fid]
        d,i=stations.query(trial,k=min(len(offsets),5))
        if d.ndim==1:
            d,i=d[:,None],i[:,None]
        weights=np.exp(-.5*(d/1.9)**2)
        shift=np.sum(offsets[i]*weights[...,None],axis=1)/np.maximum(weights.sum(axis=1),1.e-6)[:,None]
        t=np.clip((4.-d[:,0])/2.5,0.,1.)
        if self.endpoints is not None:
            p,_=self.endpoints.query(trial)
            t*=np.clip((p-1.)/2.1,0.,1.)
        trial=trial+shift*(t*t*(3-2*t))[:,None]
        source=trial+self.base.delta(trial)
        return self.profile.forward(self.reference.forward(source))


def warp_dem(src, fit, path, rate):
    west,south,east,north=fit.bounds
    d=1/rate
    w,h=(east-west)*rate,(north-south)*rate
    x=west+(np.arange(w)+.5)*d
    sw,sdx,_,sn,_,sdy=src.GetGeoTransform()
    dst=gdal.GetDriverByName("GTiff").Create(str(path),w,h,1,gdal.GDT_UInt16,
        options=["TILED=YES","COMPRESS=DEFLATE","PREDICTOR=2","NUM_THREADS=2"])
    dst.SetGeoTransform((west,d,0,north,0,-d))
    dst.SetProjection(local_srs().ExportToWkt())
    dst.SetMetadata({"SOURCE_URL":"https://github.com/bburns/Arda/releases/tag/dem-32k-v1",
                     "ATTRIBUTION":"Pete Lewis; ME-DEM by monks and Redrobes",
                     "ADAPTATION":"Native Arda 32k jointly registered with the GIS vectors",
                     "COORDINATE_UNITS":"one MMapper room per unit"})
    gdal.SetCacheMax(220_000_000)
    for start in range(0,h,64):
        rows=np.arange(start,min(start+64,h))
        xx,yy=np.meshgrid(x,north-(rows+.5)*d)
        native=fit.native(np.column_stack((xx.ravel(),yy.ravel())))
        col=(native[:,0]-sw)/sdx-.5
        row=(native[:,1]-sn)/sdy-.5
        l,t=int(np.floor(col.min())),int(np.floor(row.min()))
        r,b=int(np.floor(col.max()))+2,int(np.floor(row.max()))+2
        if min(l,t)<0 or r>src.RasterXSize or b>src.RasterYSize:
            raise ValueError("The fitted MMapper extent exceeds the native DEM")
        sample=src.ReadAsArray(l,t,r-l,b-t)
        values=map_coordinates(sample,[row-t,col-l],order=1,mode="nearest",output=np.uint16)
        dst.WriteArray(values.reshape(xx.shape),0,start)
    dst.FlushCache()


def local_srs():
    srs=osr.SpatialReference()
    srs.SetLocalCS("MUME room grid")
    srs.SetLinearUnits("room",1.0)
    return srs


def densify(points, step=315.):
    dist=np.linalg.norm(np.diff(points,axis=0),axis=1)
    parts=np.maximum(1,np.ceil(dist/step)).astype(int)
    starts=np.repeat(points[:-1],parts,axis=0)
    deltas=np.repeat(np.diff(points,axis=0)/parts[:,None],parts,axis=0)
    offsets=np.arange(parts.sum())-np.repeat(np.r_[0,parts.cumsum()[:-1]],parts)
    return np.vstack([starts+deltas*offsets[:,None],points[-1]])


def deform_geometry(geom, transform):
    kind=ogr.GT_Flatten(geom.GetGeometryType())
    result=ogr.Geometry(kind)
    if geom.GetGeometryName()=="LINEARRING":
        result=ogr.Geometry(ogr.wkbLinearRing)
    if geom.GetPointCount():
        points=np.asarray(geom.GetPoints())[:,:2]
        if len(points)>1:
            points=densify(points)
        for x,y in transform(points):
            result.AddPoint_2D(x,y)
    else:
        for part in geom:
            result.AddGeometry(deform_geometry(part,transform))
    return result


def warp_vectors(src_path, fit, output):
    source=ogr.Open(str(src_path))
    driver=ogr.GetDriverByName("GPKG")
    out=driver.CreateDataSource(str(output))
    west,south,east,north=fit.bounds
    coords=np.array([[west,south],[west,north],[east,south],[east,north]])
    base=(coords-DESTINATION+fit.profile.center)@TO_NATIVE+ORIGIN
    lo,hi=base.min(axis=0)-38_000,base.max(axis=0)+38_000
    if fit.source_registration:
        rivers=source.GetLayerByName("Rivers")
        rivers.SetSpatialFilterRect(*lo,*hi)
        endpoints=[]
        for feature in rivers:
            for line in feature.GetGeometryRef():
                points=np.asarray(line.GetPoints())[:,:2]
                if len(points)>1:
                    endpoints.extend(points[[0,-1]])
        if endpoints:
            original=(np.array(endpoints)-ORIGIN) @ TO_ROOM+DESTINATION
            fit.endpoints=cKDTree(fit.base.forward(original))
    for name in LAYERS:
        src=source.GetLayerByName(name)
        src.SetSpatialFilterRect(*lo,*hi)
        layer=out.CreateLayer(name,local_srs(),geom_type=ogr.wkbUnknown,
                              options=["SPATIAL_INDEX=NO"])
        for k,pattern in (("source_id",ogr.OFTInteger),("name",ogr.OFTString),("source_type",ogr.OFTString)):
            layer.CreateField(ogr.FieldDefn(k,pattern))
        layer.StartTransaction()
        definition=layer.GetLayerDefn()
        for f in src:
            geom=f.GetGeometryRef()
            if geom is None or geom.IsEmpty():
                continue
            fid=f.GetFID()
            transform=(lambda points:fit.river_target(points,fid)) if name=="Rivers" else fit.target
            target=deform_geometry(geom,transform)
            feature=ogr.Feature(definition)
            feature.SetField("source_id",f.GetFID())
            if f.GetField("name"):
                feature.SetField("name",f.GetField("name"))
            if f.GetFieldIndex("type")>=0 and f.GetField("type"):
                feature.SetField("source_type",f.GetField("type"))
            feature.SetGeometry(target)
            if layer.CreateFeature(feature):
                raise ValueError(f"Could not write feature {name} {f.GetFID()}")
        layer.CommitTransaction()
    out.FlushCache()


def layer_raster(layer,bounds,rate=SUPERSAMPLE):
    west,south,east,north=bounds
    ds=gdal.GetDriverByName("MEM").Create("",(east-west)*rate,(north-south)*rate,1,gdal.GDT_Byte)
    ds.SetGeoTransform((west,1/rate,0,north,0,-1/rate))
    ds.SetProjection(local_srs().ExportToWkt())
    layer.ResetReading()
    gdal.RasterizeLayer(ds,[1],layer,burn_values=[1])
    return ds.ReadAsArray()[::-1].copy()


def room_coverage(data,rate=SUPERSAMPLE):
    h,w=data.shape
    return data.reshape(h//rate,rate,w//rate,rate).mean(axis=(1,3))


def road_water_features(layer,bounds,locked,sectors,road=False):
    west,south,east,north=bounds
    shape=(north-south,east-west)
    cells=np.zeros(shape,bool)
    link_flags=np.zeros(shape,np.uint8)
    nodes={}
    directions={(0,1):0,(0,-1):1,(1,0):2,(-1,0):3}
    layer.ResetReading()
    # A line passing briefly through a corner still enters that room. A
    # diagonal corner tie receives one orthogonal connector.
    for f in layer:
        if road and f.GetField("source_type")=="TERTIARY":
            continue
        geom=f.GetGeometryRef()
        fid=f.GetField("source_id")
        samples=[]
        for line in geom:
            p=np.asarray(line.GetPoints())[:,:2]
            if len(p)<2:
                continue
            if not f.GetField("name") and np.linalg.norm(np.diff(p,axis=0),axis=1).sum() < (3 if road else 2.2):
                continue
            dense=densify(p,.13)
            line_cells=np.floor(dense).astype(int)-[west,south]
            keep=np.r_[True,np.any(np.diff(line_cells,axis=0),axis=1)]
            points=[tuple(p) for p in line_cells[keep]]
            previous=None
            for x,y in points:
                if not 0<=x<shape[1] or not 0<=y<shape[0]:
                    previous=None
                    continue
                cells[y,x]=True
                here=(x,y)
                samples.append(here)
                if previous is not None:
                    px,py=previous
                    chain=[previous,here]
                    if abs(px-x)+abs(py-y)==2:
                        alternatives=[(px,y),(x,py)]
                        allowed={CODE["ROAD"]} if road else set(WET)
                        middle=min(alternatives,key=lambda p:(
                            locked[p[1],p[0]] and sectors[p[1],p[0]] not in allowed,
                            np.min(np.sum((np.asarray(p)+.5+np.array([west,south])-dense)**2,axis=1)),
                        ))
                        cells[middle[1],middle[0]]=True
                        chain.insert(1,middle)
                    for (ax,ay),(bx,by) in zip(chain,chain[1:]):
                        key=(bx-ax,by-ay)
                        if key in directions:
                            k=directions[key]
                            link_flags[ay,ax]|=1<<k
                            link_flags[by,bx]|=1<<(k^1)
                previous=here
        if samples:
            nodes[fid]=samples
    return cells,link_flags,nodes


def connect_seams(sector, locked, river, seed, max_radius=5):
    """Connect mapped water at open room edges to nearby GIS river cells."""
    water=np.isin(sector,WET)
    near=distance_transform_edt(~river)
    additions=0
    unresolved=[]
    h,w=sector.shape
    frontier=seed&water
    for sy,sx in np.argwhere(frontier):
        if near[sy,sx] > max_radius:
            continue
        allowed=[]
        for dx,dy in STEPS:
            x,y=sx+dx,sy+dy
            if 0<=x<w and 0<=y<h and not locked[y,x] and sector[y,x]!=CODE["ROAD"]:
                allowed.append((x,y))
        if not allowed or any(river[y,x] and water[y,x] for x,y in allowed):
            continue
        visited={(sx,sy)}
        todo=[]
        for x,y in allowed:
            heapq.heappush(todo,(1+.13*near[y,x],1,x,y,((x,y),)))
        path=None
        while todo:
            cost,steps,x,y,route=heapq.heappop(todo)
            if (x,y) in visited or steps>max_radius:
                continue
            visited.add((x,y))
            if river[y,x] and water[y,x]:
                path=route
                break
            for dx,dy in STEPS:
                a,b=x+dx,y+dy
                if (a,b) not in visited and 0<=a<w and 0<=b<h and not locked[b,a] and sector[b,a]!=CODE["ROAD"]:
                    heapq.heappush(todo,(cost+1+.13*near[b,a],steps+1,a,b,route+((a,b),)))
        if path is not None:
            for x,y in path:
                if not water[y,x]:
                    sector[y,x]=CODE["WATER"]
                    water[y,x]=True
                    additions+=1
        else:
            unresolved.append((int(sx),int(sy)))
    return additions,unresolved


def target_mask(rooms,bounds,hidden):
    west,south,east,north=bounds
    shape=(north-south,east-west)
    known=np.zeros(shape,np.uint8)
    locked=np.zeros(shape,bool)
    flags=np.zeros(shape,np.uint8)
    for r in rooms:
        x,y,z=r.xyz
        if z or not (west<=x<east and south<=y<north) or (x,y) in hidden:
            continue
        q=y-south,x-west
        if r.outdoor:
            locked[q]=True
            known[q]=CODE[r.terrain]
            if r.terrain=="ROAD":
                flags[q]=r.flags
    return known,locked,flags


def build_sectors(dem_path, vectors_path, fit, rooms, hidden):
    bounds=fit.bounds
    west,south,east,north=bounds
    w,h=east-west,north-south
    src=ogr.Open(str(vectors_path))
    dem=gdal.Open(str(dem_path))
    known,locked,known_flags=target_mask(rooms,bounds,hidden)
    z=dem.ReadAsArray(buf_xsize=w*SUPERSAMPLE,buf_ysize=h*SUPERSAMPLE,
                      resample_alg=gdal.GRIORA_Average)[::-1]
    elev=room_coverage(z)
    relief=maximum_filter(elev,size=5)-minimum_filter(elev,size=5)
    woods=src.GetLayerByName("forests")
    woods.SetAttributeFilter("source_type IS NULL OR source_type <> 'Forest Clearing'")
    forest_samples=layer_raster(woods,bounds)
    woods.SetAttributeFilter("source_type = 'Forest Clearing'")
    forest_samples[layer_raster(woods,bounds).astype(bool)]=0
    forest=room_coverage(forest_samples)
    brush=distance_transform_edt(~(forest>.45))<=1.2
    hills=room_coverage(layer_raster(src.GetLayerByName("Hills"),bounds))
    mountain_axis=distance_transform_edt(~layer_raster(src.GetLayerByName("Mountains"),bounds).astype(bool))/SUPERSAMPLE
    mountains=(elev>19800)&((relief>5600)|(elev>28700)) | ((room_coverage((mountain_axis<2.2))>.35)&(relief>4200)&(elev>15750))
    sector=np.full((h,w),CODE["FIELD"],np.uint8)
    sector[brush&(relief<7000)&(forest>.06)]=CODE["BRUSH"]
    sector[(hills>.28)|((elev>14800)&(relief>4100))]=CODE["HILLS"]
    sector[forest>=.40]=CODE["FOREST"]
    sector[mountains]=CODE["MOUNTAINS"]
    wetland=room_coverage(layer_raster(src.GetLayerByName("Wetlands"),bounds))
    sector[wetland>=.33]=CODE["SHALLOW"]
    shore=layer_raster(src.GetLayerByName("lakes"),bounds)
    lake=room_coverage(shore)
    coast=layer_raster(src.GetLayerByName("Coastline"),bounds).astype(bool)
    # This GIS source contains coast lines. The offshore side is identified
    # from connected low native elevation at the map boundary.
    low=z<=5620
    barriers=maximum_filter(coast,size=3)
    basins,_=label(low & ~barriers,CROSS)
    edge=np.r_[basins[0],basins[-1],basins[:,0],basins[:,-1]]
    edges_z=np.r_[z[0],z[-1],z[:,0],z[:,-1]]
    ids=np.unique(edge[(edge>0)&(edges_z<5000)])
    ocean=room_coverage(np.isin(basins,ids))>=.36
    sector[(lake>=.36)|ocean]=CODE["WATER"]
    river,river_flags,river_nodes=road_water_features(src.GetLayerByName("Rivers"),bounds,locked,known)
    road,flags,_=road_water_features(src.GetLayerByName("Roads"),bounds,locked,known,True)
    sector[river]=CODE["WATER"]
    sector[road & ~np.isin(sector,WET)]=CODE["ROAD"]
    sector[locked]=known[locked]
    flags[sector!=CODE["ROAD"]]=0
    flags[locked]=known_flags[locked]
    neighbors_of_inferred=~locked
    for i,(dx,dy) in enumerate(STEPS):
        rotated=np.roll(neighbors_of_inferred,(-dy,-dx),axis=(0,1))
        if i==0:
            open_edge=np.zeros_like(locked)
        if dx==1: rotated[:,-1]=False
        elif dx==-1: rotated[:,0]=False
        elif dy==1: rotated[-1,:]=False
        else: rotated[0,:]=False
        open_edge |= rotated
    seeds=locked&np.isin(known,WET)&open_edge
    n,unresolved=connect_seams(sector,locked,river,seeds)
    flags[sector!=CODE["ROAD"]]=0
    bgwet=np.isin(sector,WET)
    labels,nwater=label(bgwet,CROSS)
    disconnect=[]
    for fid, ps in river_nodes.items():
        labels_id={int(labels[y,x]) for x,y in ps if labels[y,x]}
        if len(labels_id)>1:
            disconnect.append(fid)
    report={
        "known_outdoor_rooms":int(locked.sum()),
        "mapped_main_roads":int((locked&(known==CODE["ROAD"])).sum()),
        "outside_outdoor_rooms":sum(r.outdoor and not(west<=r.xyz[0]<east and south<=r.xyz[1]<north) for r in rooms),
        "explicitly_hidden_surface_coordinates":len(hidden),
        "terrain_counts":dict(Counter(TERRAIN[i] for i in sector.ravel())),
        "river_features_sampled":len(river_nodes),
        "candidate_river_fragments_for_review":disconnect,
        "water_components":nwater,
        "mapped_seam_repair_cells":n,
        "unresolved_near_seam_candidates":[[x+west,y+south] for x,y in unresolved],
    }
    return sector,flags,report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ["dem","gis","reference","target","profile-fit","output"]:
        p.add_argument("--"+arg,required=True,type=Path)
    here=Path(__file__).parent
    p.add_argument("--terrain-controls",type=Path,default=here/"reference-terrain-controls.csv")
    p.add_argument("--joins",type=Path,default=here/"reference-joins.csv")
    p.add_argument("--exclude-outdoor",type=Path,help="CSV with x,y columns in target rooms")
    p.add_argument("--reach-fit",type=Path,help="Folder containing the approved reference reach fit")
    p.add_argument("--dem-pixels",type=int,default=16)
    p.add_argument("--hillshade-pixels",type=int,default=11)
    p.add_argument("--resume",action="store_true",help="Reuse the DEM and vectors in this output directory")
    args=p.parse_args()
    gdal.UseExceptions()
    ogr.UseExceptions()
    out=args.output
    out.mkdir(exist_ok=True,parents=True)
    target=read_map(args.target)
    reference=read_map(args.reference)
    profile=np.load(args.profile_fit)
    if int(profile["rate"])!=RATE:
        raise ValueError("Unsupported profile displacement grid")
    profile_report=json.loads((args.profile_fit.parent/"profile-report.json").read_text())
    if digest(args.target)!=profile_report["source_files"][profile_report["target"]]:
        raise ValueError("The chosen coordinate profile belongs to a different playable map")
    base,check=make_reference_fit(args.terrain_controls,args.joins,reference)
    field=base
    selected=[]
    if args.reach_fit:
        field=np.load(args.reach_fit/"combined-reference-pull.npz")["offset"]
        report=json.loads((args.reach_fit/"reach-fit.json").read_text())
        selected=report["selected_stations"]
        check["river_reach_fit"]=report["deformation"]
    fit=Fit(field,Pull(profile["offset"],profile["bounds"],profile["center"]),base,selected)
    np.savez_compressed(out/"reference-pull.npz",offset=field,bounds=BOUNDS,rate=RATE)
    dem=out/"elevation.tif"
    gis=out/"aligned-gis.gpkg"
    if not args.resume or not dem.exists():
        print("Warping original 32k elevation",flush=True)
        warp_dem(gdal.Open(str(args.dem)),fit,dem,args.dem_pixels)
    if not args.resume or not gis.exists():
        if gis.exists():
            ogr.GetDriverByName("GPKG").DeleteDataSource(str(gis))
        print("Warping original GIS layers",flush=True)
        warp_vectors(args.gis,fit,gis)
    hidden=set()
    if args.exclude_outdoor:
        with args.exclude_outdoor.open(newline="") as f:
            hidden={(int(row["x"]),int(row["y"])) for row in csv.DictReader(f)}
    print("Rasterizing GIS and applying foreground sectors",flush=True)
    sector,flags,stats=build_sectors(dem,gis,fit,target,hidden)
    name=profile_report["profile"]
    binary=out/f"{name}-background-river-reaches.mm2"
    write_mm2(binary,sector,flags,fit.bounds)
    binary_check=verify_mm2(binary,sector,flags,fit.bounds)
    np.savez_compressed(out/"sectors-review.npz",sector=sector,bounds=fit.bounds)
    shade=out/f"{name}.hillshade.png"
    if not args.resume or not shade.exists():
        bake(dem,shade,args.hillshade_pixels)
    result={
        "profile":name,"bounds":fit.bounds,"reference_sha256":digest(args.reference),
        "playable_sha256":digest(args.target),"profile_fit_sha256":digest(args.profile_fit),
        "geo_package_sha256":digest(args.gis),"native_vrt_sha256":digest(args.dem),
        "terrain_controls_sha256":digest(args.terrain_controls),"joins_sha256":digest(args.joins),
        "fit":check,"background":stats,
        "mm2_verification":binary_check,
        "outputs_sha256":{a.name:digest(a) for a in [dem,gis,binary,shade,shade.with_suffix(".pgw")]},
    }
    (out/"build-report.json").write_text(json.dumps(result,indent=2)+"\n")
    print(name,json.dumps({k:v for k,v in stats.items() if not isinstance(v,list)}),flush=True)


if __name__=="__main__":
    main()
