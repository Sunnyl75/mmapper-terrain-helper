"""Recover reviewed matching terrain controls from the archived bundled elevation."""
import argparse
import csv
from pathlib import Path

import numpy as np
from osgeo import gdal
from scipy.ndimage import distance_transform_edt, map_coordinates
from scipy.signal import correlate

from room_io import read_map

BOUNDS = (-6, -271, 721, 21)
ORIGIN = np.array([1272120.99119484, 1154992.61815864])
DESTINATION = np.array([635.035, -9.4736])
SCALE = complex(.0006377545820894734, -4.456752403802747e-7)
TO_ROOM = np.array([[SCALE.real, SCALE.imag], [-SCALE.imag, SCALE.real]])
TO_NATIVE = np.linalg.inv(TO_ROOM)


def source_sample(dem, points, overview=0):
    west, dx, _, north, _, dy = dem.GetGeoTransform()
    native = (points - DESTINATION) @ TO_NATIVE + ORIGIN
    px = (native[:, 0] - west) / dx
    py = (native[:, 1] - north) / dy
    band = dem.GetRasterBand(1)
    if overview:
        band = band.GetOverview(overview - 1)
        px *= band.XSize / dem.RasterXSize
        py *= band.YSize / dem.RasterYSize
    px -= .5
    py -= .5
    left = max(0, int(np.floor(px.min())))
    top = max(0, int(np.floor(py.min())))
    right = min(band.XSize, int(np.floor(px.max())) + 2)
    bottom = min(band.YSize, int(np.floor(py.max())) + 2)
    if min(px.min(), py.min()) < -.5 or px.max() > band.XSize-.5 or py.max() > band.YSize-.5:
        raise ValueError("The fitted grid exceeds the source DEM")
    pixels = band.ReadAsArray(left, top, right-left, bottom-top)
    return map_coordinates(pixels.astype(np.float32), [py-top, px-left], order=1, mode="nearest")


def geom_grid(bounds, rate):
    west, south, east, north = bounds
    x = west + (np.arange(round((east-west)*rate))+.5)/rate
    y = south + (np.arange(round((north-south)*rate))+.5)/rate
    return x, y


def derive(native_path, archived_path, map_path, output):
    native = gdal.Open(str(native_path))
    archived = gdal.Open(str(archived_path))
    west, south, east, north = BOUNDS
    rate = 4
    x, y = geom_grid(BOUNDS, rate)
    xx, yy = np.meshgrid(x, y)
    raw = source_sample(native, np.column_stack((xx.ravel(), yy.ravel())), overview=2).reshape(xx.shape)
    aligned = archived.ReadAsArray(buf_xsize=len(x), buf_ysize=len(y),
                                   resample_alg=gdal.GRIORA_Average)[::-1].astype(np.float32)
    active = np.zeros((north-south, east-west), dtype=bool)
    for room in read_map(map_path):
        x, y, _ = room.xyz
        if room.outdoor and west<=x<east and south<=y<north:
            active[y-south,x-west] = True
    proximity = distance_transform_edt(~active)
    half, shift = 16, 38
    rows = []
    for y in range(south+15, north-14, 8):
        for x in range(west+15, east-14, 8):
            dist = proximity[y-south,x-west]
            if dist > 13:
                continue
            iy, ix = (y-south)*rate, (x-west)*rate
            t = aligned[iy-half:iy+half,ix-half:ix+half]
            p = raw[iy-half-shift:iy+half+shift,ix-half-shift:ix+half+shift]
            if t.std()<130 or np.ptp(t)<480:
                continue
            kernel = np.ones_like(t)
            z = correlate(p, kernel, mode="valid", method="fft")
            sq = correlate(p*p, kernel, mode="valid", method="fft")
            cov = correlate(p, t-t.mean(), mode="valid", method="fft")
            scores = cov / np.sqrt(np.maximum(sq-z*z/t.size,1)*np.sum((t-t.mean())**2))
            dr,dc = np.unravel_index(scores.argmax(),scores.shape)
            score = float(scores[dr,dc])
            masked = scores.copy()
            masked[max(0,dr-5):dr+6,max(0,dc-5):dc+6] = -1
            gap = score-masked.max()
            if score < .92 or gap < .018 or min(dr,dc)==0 or max(dr,dc)==2*shift:
                continue
            sx,sy = x+(dc-shift)/rate, y+(dr-shift)/rate
            rows.append((x,y,sx,sy,score,gap,dist))
    with output.open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["target_x","target_y","source_x","source_y","correlation","distinctiveness","distance_to_mapped_surface"])
        w.writerows(rows)
    print(len(rows), "recovered terrain correspondences")
    print("median and p95 shifts:",np.quantile([np.hypot(a-c,b-d) for a,b,c,d,*_ in rows],[.5,.95]))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("dem","archived","reference","output"):
        p.add_argument(arg,type=Path)
    a=p.parse_args()
    derive(a.dem,a.archived,a.reference,a.output)
