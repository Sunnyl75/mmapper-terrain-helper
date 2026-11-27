"""Match reference rooms and select an inverse coordinate fit using spatial holdouts."""
from collections import Counter, defaultdict
from dataclasses import dataclass
import csv

import numpy as np
from scipy.spatial import cKDTree
from scipy.ndimage import distance_transform_edt

from room_io import fingerprint

MODELS=("translation","local-4","local-8","local-15")
END=0xffffffff


@dataclass
class Pair:
    reference: tuple
    target: tuple
    evidence: str
    ref_id: int
    target_id: int
    name: str


def unique(rooms,key):
    hits=defaultdict(list)
    for room in rooms:
        value=key(room)
        if value is not None:
            hits[value].append(room)
    return {k:v[0] for k,v in hits.items() if len(v)==1}


def known_id(r):
    return r.server_id if 0<r.server_id<END else None


def match(reference,target,manual):
    rs,ts=unique(reference,known_id),unique(target,known_id)
    ids=[(rs[k],ts[k],"server-id") for k in sorted(rs.keys() & ts.keys())
         if rs[k].outdoor and ts[k].outdoor]
    used_a={a.id for a,_,_ in ids}
    used_b={b.id for _,b,_ in ids}
    rf,tf=unique(reference,fingerprint),unique(target,fingerprint)
    text=[]
    for key in sorted(rf.keys() & tf.keys()):
        a,b=rf[key],tf[key]
        if (a.outdoor and b.outdoor and a.id not in used_a and b.id not in used_b
            and not (known_id(a) and known_id(b) and known_id(a)!=known_id(b))):
            text.append((a,b,"unique-name-description"))

    # A matching exit to another paired room aids the ruling on isolated outliers.
    pairs=ids+text
    mapping={a.id:b.id for a,b,_ in pairs}
    topology=[]
    for a,b,_ in pairs:
        topology.append(sum(bool(set(b.exits.get(direction,())) &
                                 {mapping[i] for i in linked if i in mapping})
                            for direction,linked in a.exits.items()
                            if direction in ("north","south","east","west")))
    out=[Pair(a.xyz[:2],b.xyz[:2],why,a.id,b.id,b.name) for a,b,why in pairs]
    # Require one-to-one coordinates, including for manually reviewed pairs.
    if manual:
        ac={r.xyz[:2]:r for r in reference if r.outdoor}
        bc={r.xyz[:2]:r for r in target if r.outdoor}
        with manual.open(newline="",encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                a=tuple(int(row["reference_"+k]) for k in "xy")
                b=tuple(int(row["target_"+k]) for k in "xy")
                reason=row["reason"].strip()
                if a not in ac or b not in bc or not reason:
                    raise ValueError("Manual controls need existing outdoor rooms and a reason")
                if known_id(ac[a]) and known_id(bc[b]) and known_id(ac[a])!=known_id(bc[b]):
                    raise ValueError("Manual control contradicts known server IDs")
                keep=[i for i,p in enumerate(out) if p.reference!=a and p.target!=b]
                out=[out[i] for i in keep]
                topology=[topology[i] for i in keep]
                out.append(Pair(a,b,"manual",ac[a].id,bc[b].id,reason))
                topology.append(4)
    counts_a=Counter(r.reference for r in out)
    counts_b=Counter(r.target for r in out)
    retained=[i for i,r in enumerate(out) if counts_a[r.reference]==counts_b[r.target]==1]
    return [out[i] for i in retained],np.array([topology[i] for i in retained])


def robust_controls(pairs,topology):
    left=np.array([r.reference for r in pairs],dtype=float)+.5
    right=np.array([r.target for r in pairs],dtype=float)+.5
    offsets=left-right
    distance,index=cKDTree(right).query(right,k=min(13,len(right)))
    closest=offsets[index[:,1:]]
    med=np.median(closest,axis=1)
    deviation=np.linalg.norm(offsets-med,axis=1)
    neighbor_spread=np.median(np.linalg.norm(closest-med[:,None],axis=2),axis=1)
    agree=np.sum(np.linalg.norm(closest-offsets[:,None],axis=2)<=1.0,axis=1)
    limit=np.maximum(3.0,3.5*neighbor_spread+1)
    manual=np.array([p.evidence=="manual" for p in pairs])
    remove=(deviation>limit)&(agree<2)&(topology==0)&(~manual)&(distance[:,-1]<24)
    return left,right,remove


class Fit:
    def __init__(self,target,reference,model,global_offset=None):
        self.samples=target
        self.displacement=reference-target
        if global_offset is None:
            bucket=np.floor(target/8).astype(int)
            groups=defaultdict(list)
            for key,delta in zip(map(tuple,bucket),self.displacement):
                groups[key].append(delta)
            self.center=np.rint(np.median([np.median(v,axis=0) for v in groups.values()],axis=0))
        else:
            self.center=np.asarray(global_offset,dtype=float)
        self.width=None if model=="translation" else float(model.partition("-")[2])
        self.tree=cKDTree(target) if self.width else None

    def shift(self,points):
        if self.tree is None:
            return np.broadcast_to(self.center,points.shape).copy()
        n=min(len(self.samples),36)
        distances,indices=self.tree.query(points,k=n)
        residual=self.displacement[indices]-self.center
        weights=np.exp(-.5*((distances-distances[:,0,None])/self.width)**2)/(distances+1)
        weights/=weights.sum(axis=1)[:,None]
        out=np.einsum("ni,nij->nj",weights,residual)
        # Outside measured regions, return smoothly to the map-wide offset.
        fade=np.clip((50-distances[:,0])/30,0,1)
        fade=fade**2*(3-2*fade)
        return self.center+out*fade[:,None]


def summary(errors):
    return {"count":int(len(errors)),"mean":float(np.mean(errors)),
            "median":float(np.median(errors)),"p90":float(np.quantile(errors,.9)),
            "p95":float(np.quantile(errors,.95)),"over_two_rooms":int(np.count_nonzero(errors>2))}


def blocked_validation(target,reference):
    block=np.floor(target/16).astype(int)
    split=(31*block[:,0]+17*block[:,1])%5
    scores={}
    errors_by_method={}
    for model in MODELS:
        errors=np.full(len(target),np.nan)
        for fold in range(5):
            test=split==fold
            if not test.any() or (~test).sum()<10:
                continue
            model_fit=Fit(target[~test],reference[~test],model)
            p=target[test]
            errors[test]=np.linalg.norm(p+model_fit.shift(p)-reference[test],axis=1)
        ok=np.isfinite(errors)
        if ok.sum()!=len(target):
            raise ValueError("Insufficient spatial coverage for blocked validation")
        v=defaultdict(list)
        for i,key in enumerate(map(tuple,block)):
            v[key].append(errors[i])
        scores[model]=summary(errors)
        scores[model]["mean_of_block_means"]=float(np.mean([np.mean(x) for x in v.values()]))
        errors_by_method[model]=errors
    return scores,errors_by_method,split


def choose(scores):
    best=min(MODELS,key=lambda m:scores[m]["mean_of_block_means"])
    return "translation" if (scores["translation"]["mean_of_block_means"]<=scores[best]["mean_of_block_means"]+.04
                             and scores["translation"]["p95"]<=scores[best]["p95"]+.15) else best


def ease(t):
    t=np.clip(t,0,1)
    return t*t*(3-2*t)


def determinant(vx,vy,interval):
    # Components are bilinearly interpolated at runtime. Checking each
    # corner covers the extrema of the affine determinant in every cell.
    dx=(vx[:,1:]-vx[:,:-1])/interval
    ex=(vy[:,1:]-vy[:,:-1])/interval
    dy=(vx[1:]-vx[:-1])/interval
    ey=(vy[1:]-vy[:-1])/interval
    j=min(float(((1+dx[rows])*(1+ey[:,cols])-dy[:,cols]*ex[rows]).min())
          for rows in (slice(1,None),slice(None,-1)) for cols in (slice(1,None),slice(None,-1)))
    return j


def sampled_field(model,bounds,rate=2):
    west,south,east,north=bounds
    shape=((north-south)*rate+1,(east-west)*rate+1)
    source_center=model.center
    x=np.linspace(west,east,shape[1])
    y=np.linspace(south,north,shape[0])
    field=np.empty((*shape,2),dtype=np.float32)
    for start in range(0,len(y),24):
        xx,yy=np.meshgrid(x,y[start:start+24])
        p=np.column_stack((xx.ravel(),yy.ravel()))
        residual=(model.shift(p)-source_center).reshape((*xx.shape,2))
        margin=np.minimum.reduce((xx-west,east-xx,yy-south,north-yy))
        field[start:start+len(yy)]=residual*ease(margin/12)[:,:,None]
    amplitude=np.ones(shape,dtype=np.float32)
    damped=np.zeros(shape,dtype=bool)
    for _ in range(30):
        vx,vy=field[...,0],field[...,1]
        dx=(vx[:,1:]-vx[:,:-1])*rate
        ex=(vy[:,1:]-vy[:,:-1])*rate
        dy=(vx[1:]-vx[:-1])*rate
        ey=(vy[1:]-vy[:-1])*rate
        cells=np.minimum.reduce([(1+dx[rows])*(1+ey[:,cols])-dy[:,cols]*ex[rows]
             for rows in (slice(1,None),slice(None,-1)) for cols in (slice(1,None),slice(None,-1))])
        j=float(cells.min())
        if j>=.30:
            return field,{"minimum_inverse_jacobian":j,"minimum_local_regional_scale":float(amplitude.min()),
                          "field_grid_samples_damped":int(damped.sum())}
        bad=np.zeros(shape,dtype=bool)
        b=cells<.30
        for y,x in ((slice(1,None),slice(1,None)),(slice(None,-1),slice(1,None)),
                    (slice(1,None),slice(None,-1)),(slice(None,-1),slice(None,-1))):
            bad[y,x]|=b
        close=distance_transform_edt(~bad)/rate
        atten=.78+.22*ease(close/10)
        field*=atten[...,None]
        amplitude*=atten
        damped|=atten<1.
    raise ValueError("Regional control field could not be regularized without folds")
