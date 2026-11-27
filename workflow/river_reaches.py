"""Fit bounded, cross-channel displacement along MUME-supported river reaches."""
import argparse
from collections import defaultdict
from dataclasses import dataclass
import csv
import json
from pathlib import Path
import re

import numpy as np
from osgeo import gdal, ogr
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import distance_transform_edt, gaussian_filter, gaussian_filter1d, map_coordinates
from scipy.signal import find_peaks
from scipy.spatial import cKDTree

from room_io import read_map
from fit_reference import BOUNDS

RATE = 2
CROSS = np.arange(-6., 6.001, .08)
ALIASES = {
    "Anduin": r"anduin|great river", "Bruinen": r"bruinen|loudwater",
    "Mitheithel": r"mitheithel|hoarwell", "Baranduin": r"baranduin|brandywine",
    "Celebrant": r"celebrant|silverlode", "Gwathló": r"gwathl[oó]|greyflood",
    "Lune": r"lune|lh[uû]n", "Gladden River": r"gladden",
    "River Withywindle": r"withywindle", "River Shirebourn": r"shirebourn",
}
EXCLUDE = re.compile(r"lake|pool|pond|swamp|marsh|bog|bay|harbour|dock|side channel|branch|tributary|shore|bank|mouth|delta|weir|confluence|joins|meeting|crossing|ford|bridge",re.I)
LEGACY = Path(__file__).parent


@dataclass
class PathGeometry:
    fid: int
    name: str
    xy: np.ndarray
    s: np.ndarray
    normal: np.ndarray
    tree: cKDTree


@dataclass
class Station:
    s: float
    xy: np.ndarray
    normal: np.ndarray
    trough: float
    gis: float
    source: str


def gradual(x):
    x=np.clip(x,0,1)
    return x*x*(3-2*x)


def layer_paths(gis,bounds=BOUNDS):
    ds=ogr.Open(str(gis))
    out=[]
    for feat in ds.GetLayerByName("Rivers"):
        name=feat.GetField("name") or ""
        for g in feat.GetGeometryRef():
            line=np.array(g.GetPoints())[:,:2]
            if len(line)<2:
                continue
            dist=np.r_[0,np.cumsum(np.linalg.norm(np.diff(line,axis=0),axis=1))]
            if dist[-1]<5:
                continue
            s=np.linspace(0,dist[-1],int(np.ceil(dist[-1]/.23))+1)
            xy=np.column_stack([np.interp(s,dist,p) for p in line.T])
            w,so,e,n=bounds
            if not np.any((w<xy[:,0])&(xy[:,0]<e)&(so<xy[:,1])&(xy[:,1]<n)):
                continue
            tangent=np.gradient(gaussian_filter1d(xy,2.8,axis=0,mode="nearest"),axis=0)
            normal=np.column_stack([-tangent[:,1],tangent[:,0]])
            normal/=np.maximum(1.e-8,np.linalg.norm(normal,axis=1))[:,None]
            out.append(PathGeometry(feat.GetField("source_id"),name,xy,s,normal,cKDTree(xy)))
    return out


def set_up_associations(paths,rooms):
    points=np.vstack([p.xy for p in paths])
    owners=np.repeat(np.arange(len(paths)),[len(p.xy) for p in paths])
    tree=cKDTree(points)
    names={p.name for p in paths if p.name}
    patterns={n:re.compile(r"\b(?:"+ALIASES.get(n,re.escape(n))+r")\b",re.I) for n in names}
    by_name=defaultdict(list)
    for i,p in enumerate(paths):
        by_name[p.name].append(i)
    groups=defaultdict(list)
    status=defaultdict(int)
    eligible=[r for r in rooms if r.outdoor and r.terrain in ("WATER","RAPIDS")
              and not EXCLUDE.search(r.name)]
    candidates={(r.xyz[0],r.xyz[1]):r for r in eligible}
    neighbour_seen={}
    for room in eligible:
        x,y,z=room.xyz
        xy=np.array([x+.5,y+.5])
        matched=[n for n,p in patterns.items() if p.search(room.name)]
        if len(matched)>1:
            status["ambiguous_name"]+=1
            continue
        if matched:
            opts=[]
            for idx in by_name[matched[0]]:
                dist, ix=paths[idx].tree.query(xy)
                opts.append((dist,idx,ix))
            distance,i,ix=min(opts)
            if distance>5.4:
                status["named_but_far_from_matching_path"]+=1
                continue
            evidence="room name"
        else:
            dist,idx=tree.query(xy,k=48,distance_upper_bound=4.6)
            valid=np.isfinite(dist)
            if not valid.any():
                continue
            seen={}
            for dd,v in zip(dist[valid],idx[valid]):
                seen.setdefault(owners[v],dd)
            ordered=sorted(seen,key=seen.get)
            if seen[ordered[0]]>2.7 or (len(ordered)>1 and seen[ordered[1]]-seen[ordered[0]]<.85):
                status["ambiguous_nearby_channels"]+=1
                continue
            i=ordered[0]
            _,ix=paths[i].tree.query(xy)
            evidence="nearby connected water rooms"
        neighbour_seen[(x,y)]=(i,ix,xy,evidence)
    for (x,y),(i,ix,xy,evidence) in neighbour_seen.items():
        if evidence!="room name":
            # Require an actual reciprocal room link inside the candidate reach.
            r=candidates[(x,y)]
            any_neighbor=False
            for dx,dy,key,back in ((0,1,"north","south"),(0,-1,"south","north"),
                                   (1,0,"east","west"),(-1,0,"west","east")):
                q=(x+dx,y+dy)
                other=candidates.get(q)
                if other and q in neighbour_seen and neighbour_seen[q][0]==i and (
                        other.id in r.exits.get(key,()) and r.id in other.exits.get(back,())):
                    any_neighbor=True
            if not any_neighbor:
                status["isolated_geometry_only"]+=1
                continue
        groups[i].append((ix,xy,evidence))
        status[evidence]+=1
    status["water_rooms_considered"]=len(eligible)
    return groups,status


class Elevation:
    def __init__(self,path):
        src=gdal.Open(str(path))
        self.transform=src.GetGeoTransform()
        self.z=src.ReadAsArray(buf_type=gdal.GDT_Float32)
        gaussian_filter(self.z,.18/self.transform[1],output=self.z)

    def measure(self,xy,n,hint,limit):
        west,step,_,north,_,dy=self.transform
        line=xy+CROSS[:,None]*n
        if any((line[:,0]<west)|(line[:,0]>=west+self.z.shape[1]*step)
               |(line[:,1]>north)|(line[:,1]<=north+self.z.shape[0]*dy)):
            return None
        profile=map_coordinates(self.z,[(line[:,1]-north)/dy-.5,(line[:,0]-west)/step-.5],
                                order=1,mode="nearest")
        candidates=[]
        for j in find_peaks(-profile)[0]:
            if not 9<=j<len(profile)-9:
                continue
            relief=min(profile[max(j-15,0):j].max(),profile[j+1:j+16].max())-profile[j]
            delta=abs(CROSS[j]-hint)
            if relief>=120 and delta<=limit:
                candidates.append((delta- min(relief,1500)/4200,float(CROSS[j])))
        return min(candidates)[1] if candidates else None


def side_intersection(path,xy,n,ix):
    # Intersect a straight local tangent with the same across-channel line
    # used for the DEM sample. Reject a nearly parallel coincidence.
    t=np.array([n[1],-n[0]])
    part=path.xy[max(0,ix-8):ix+9]-xy
    across=part@n
    along=part@t
    sign_changes=np.flatnonzero(np.signbit(along[1:])!=np.signbit(along[:-1]))
    if len(sign_changes):
        j=min(sign_changes,key=lambda k:abs(across[k]))
        a=abs(along[j])/(abs(along[j])+abs(along[j+1]))
        return float(across[j]*(1-a)+across[j+1]*a)
    return float((path.xy[ix]-xy)@n)


def surveyed_stations(paths,elevation,survey,old_ford):
    by_name=defaultdict(list)
    for i,p in enumerate(paths):
        by_name[p.name].append(i)
    result=defaultdict(list)
    for path,old in ((survey,False),(old_ford,True)):
        with path.open(newline="") as f:
            for row in csv.DictReader(f):
                if old:
                    if not row["dem_trough_baseline_x"]:
                        continue
                    name="Anduin"
                    xy=np.array([float(row["official_channel_x"]),float(row["cross_section_y"])])
                    normal=np.array([1.,0.])
                    hint=float(row["dem_trough_baseline_x"])-xy[0]
                else:
                    name=row["river"]
                    xy=np.array([float(row["target_x"]),float(row["target_y"])])
                    normal=np.array([float(row["nx"]),float(row["ny"])])
                    hint=float(row["pull_offset"])
                possible=[(*paths[i].tree.query(xy),i) for i in by_name[name]]
                if not possible:
                    continue
                distance,idx,i=min(possible)
                if distance>5.9:
                    continue
                trough=elevation.measure(xy,normal,hint,3.1)
                if trough is None:
                    continue
                p=paths[i]
                if np.dot(normal,p.normal[idx])<0:
                    normal=-normal
                    trough=-trough
                gis=side_intersection(p,xy,normal,idx)
                result[i].append(Station(float(p.s[idx]),xy,normal,trough,gis,
                                         "Old Ford survey" if old else "previous cross-section survey"))
    return result


def automatic_stations(path, hits, elevation, fixed):
    by_bin=defaultdict(list)
    for ix,xy,evidence in hits:
        by_bin[int(path.s[ix]/1.5)].append((ix,xy,evidence))
    result=[]
    for k,group in sorted(by_bin.items()):
        i=int(np.median([h[0] for h in group]))
        s=float(path.s[i])
        if any(abs(s-z.s)<1.3 for z in fixed):
            continue
        normal=path.normal[i]
        point=path.xy[i]
        lateral=float((np.median([a[1] for a in group],axis=0)-point)@normal)
        xy=point+lateral*normal
        measure=elevation.measure(xy,normal,-lateral,2.45)
        if measure is None or abs(measure)>5.7:
            continue
        evidence="room name" if any(h[2]=="room name" for h in group) else "nearby connected water rooms"
        result.append(Station(s,xy,normal,measure,-lateral,evidence))
    return result


def estimate(s,values,q,method):
    if method=="cubic":
        return PchipInterpolator(s,values)(np.clip(q,s[0],s[-1]))
    steps=np.linspace(s[0],s[-1],max(2,int((s[-1]-s[0])*10)+1))
    return np.interp(q,steps,gaussian_filter1d(np.interp(steps,s,values),1.25/.10,mode="nearest"))


def approve(stations):
    s=np.array([r.s for r in stations])
    elevation=np.array([r.trough for r in stations])
    gis=np.array([r.gis for r in stations])
    surveyed=any("survey" in a.source for a in stations)
    named=sum(r.source=="room name" for r in stations)
    required=6 if surveyed else 7 if named>=len(stations)/2 else 10
    if len(s)<required or s[-1]-s[0]<7.5:
        return None,"insufficient contiguous support"
    # Test contiguous blocks to avoid treating adjacent samples as independent.
    blocks=((s-s[0])//3.3).astype(int)
    if len(set(blocks))<3:
        return None,"fewer than three separate station blocks"
    outcomes={}
    for method in ("cubic","smooth"):
        errors=[]
        for k in range(3):
            train=blocks%3!=k
            test=(blocks%3==k)&(s>min(s[train]))&(s<max(s[train]))
            if not test.any() or train.sum()<4:
                continue
            predicted=estimate(s[train],elevation[train],s[test],method)
            errors.extend(zip(abs(elevation[test]),abs(elevation[test]-predicted)))
        if len(errors)>=4:
            before,after=np.array(errors).T
            outcomes[method]=(float(before.mean()),float(after.mean()),len(before),
                              float(np.quantile(before,.9)),float(np.quantile(after,.9)))
    if not outcomes:
        return None,"not enough held-out stations"
    method=min(outcomes,key=lambda m:outcomes[m][1])
    before,after,count,before90,after90=outcomes[method]
    if before-after<.075 or after>before*.88 or after90>before90+.12:
        return None,"held-out fit did not improve reliably"
    return {"method":method,"block_count":len(set(blocks)),"heldout_stations":count,
            "heldout_mean_before":before,"heldout_mean_after":after,
            "heldout_p90_before":before90,"heldout_p90_after":after90},None


def deformations(selected,base):
    west,south,east,north=BOUNDS
    shape=base.shape[:2]
    delta=np.zeros_like(base)
    weights=np.zeros(shape,np.float32)
    for path,stations,choice in selected:
        s=np.array([r.s for r in stations])
        xy=np.array([r.xy for r in stations])
        normals=np.array([r.normal for r in stations])
        mode=choice["method"]
        offsets=np.array([r.trough for r in stations])
        ref=np.linspace(s[0],s[-1],int((s[-1]-s[0])*3)+2)
        pts=np.column_stack([np.interp(ref,s,p) for p in xy.T])
        direction=np.column_stack([np.interp(ref,s,q) for q in normals.T])
        corrections=estimate(s,offsets,ref,mode)[:,None]*direction
        pad=9
        lo=np.maximum([0,0],np.floor((pts.min(0)-pad-[west,south])*RATE).astype(int))
        hi=np.minimum(shape[::-1],np.ceil((pts.max(0)+pad-[west,south])*RATE).astype(int)+1)
        xx,yy=np.meshgrid(west+np.arange(lo[0],hi[0])/RATE,south+np.arange(lo[1],hi[1])/RATE)
        q=np.c_[xx.ravel(),yy.ravel()]
        dist,k=cKDTree(pts).query(q)
        # Near the channel use the complete correction. Round off the end caps
        # and preserve adjoining terrain outside the corridor.
        w=gradual((7.5-dist)/6.4)
        cap=np.minimum(ref[k]-ref[0],ref[-1]-ref[k])
        end=.25+.75*gradual(cap/3.5)
        weight=w*end
        patch=(slice(lo[1],hi[1]),slice(lo[0],hi[0]))
        weights[patch]+=weight.reshape(xx.shape)
        delta[patch]+=(corrections[k]*weight[:,None]).reshape((*xx.shape,2))
    delta/=np.maximum(1,weights)[...,None]
    delta=gaussian_filter(delta,(.60*RATE,.60*RATE,0),mode="constant")
    yy,xx=np.indices(shape)
    edge=np.minimum.reduce([yy,xx,shape[0]-1-yy,shape[1]-1-xx])/RATE
    delta*=gradual(edge/4)[...,None]
    return delta


def minimum_det(field):
    dy,dx=np.gradient(field,.5,.5,axis=(0,1))
    return (1+dx[...,0])*(1+dy[...,1])-dx[...,1]*dy[...,0]


def compose(base, shift):
    yy,xx=np.indices(base.shape[:2],dtype=np.float32)
    return shift+np.stack([
        map_coordinates(base[...,k],[yy+RATE*shift[...,1],xx+RATE*shift[...,0]],
                        order=1,mode="nearest") for k in range(2)],axis=-1)


def bound_deformation(base, shift):
    live=shift.copy()
    for n in range(10):
        combined=compose(base,live)
        determinant=minimum_det(combined)
        danger=determinant<.22
        if not danger.any():
            return live,combined,{"minimum_combined_inverse_determinant":float(determinant.min()),
                "minimum_river_inverse_determinant":float(minimum_det(live).min()),
                "local_damping_iterations":n}
        distances=distance_transform_edt(~danger)/RATE
        live*=(.76+.24*gradual(distances/5.0))[...,None]
    raise ValueError("Could not confine the river correction to a nonfolding field")


def run(args):
    room=read_map(args.map)
    paths=layer_paths(args.gis)
    groups,status=set_up_associations(paths,room)
    ds=Elevation(args.dem)
    fixed=surveyed_stations(paths,ds,args.surveyed,args.old_ford)
    selected,untouched=[],[]
    for i,p in enumerate(paths):
        surveyed=fixed.get(i,[])
        added=automatic_stations(p,groups.get(i,[]),ds,surveyed)
        by_station={}
        for station in added+surveyed:
            by_station[round(station.s/.42)]=station
        all_=sorted(by_station.values(),key=lambda r:r.s)
        segments=[]
        for station in all_:
            if not segments or station.s-segments[-1][-1].s>5.5:
                segments.append([])
            segments[-1].append(station)
        if not segments:
            untouched.append({"source_id":p.fid,"name":p.name,"reason":"no reliable surface channel stations"})
        for seg in segments:
            validate,reason=approve(seg)
            if validate:
                selected.append((p,seg,validate))
            else:
                untouched.append({"source_id":p.fid,"name":p.name,"reason":reason,"stations":len(seg)})
    base=np.load(args.reference_fit)["offset"]
    field=deformations(selected,base)
    field,combined,deformation=bound_deformation(base,field)
    vectors=[]; records=[]
    for p,segment,scores in selected:
        for item in segment:
            vectors.append({"source_id":p.fid,"xy":item.xy.tolist(),"normal":item.normal.tolist(),
                            "trial_trough":item.trough,"trial_gis":item.gis,"station":item.s,
                            "evidence":item.source})
        records.append({"source_id":p.fid,"name":p.name or "Unnamed GIS channel",
                        "mapped_stations":len(segment),"start":segment[0].xy.tolist(),
                        "end":segment[-1].xy.tolist(), **scores})
    args.output.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(args.output/"river-pull.npz",offset=field,bounds=BOUNDS,rate=RATE)
    np.savez_compressed(args.output/"combined-reference-pull.npz",offset=combined,bounds=BOUNDS,rate=RATE)
    report={"source_features_scanned":len({p.fid for p in paths}),"associations":dict(status),
            "selected_reaches":records,"unchanged_or_unresolved_reaches":untouched,
            "selected_stations":vectors,"deformation":deformation}
    (args.output/"reach-fit.json").write_text(json.dumps(report,indent=2)+"\n")
    n=defaultdict(int)
    for r in records:
        n[r["name"]]+=1
    print("Scanned",report["source_features_scanned"],"reaches accepted",len(records),
          "distinct source features",len({r["source_id"] for r in records}),"stations",len(vectors))
    print(dict(n))
    return report


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("map","dem","gis","reference-fit","output"):
        p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--surveyed",type=Path,default=LEGACY/"surveyed-river-cross-sections.csv")
    p.add_argument("--old-ford",type=Path,default=LEGACY/"surveyed-old-ford.csv")
    gdal.UseExceptions()
    ogr.UseExceptions()
    run(p.parse_args())
