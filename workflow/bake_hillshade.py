#!/usr/bin/env python3
"""Bake the plain MMapper hillshade PNG and PGW from a fitted UInt16 GeoTIFF."""

import argparse
import math
from pathlib import Path

import numpy as np
from osgeo import gdal
from PIL import Image, PngImagePlugin
from scipy.ndimage import gaussian_filter

MAX_SIDE = 8192
MAX_PIXELS = 32 * 1024 * 1024
SOURCE_URL = "https://github.com/bburns/Arda/releases/tag/dem-32k-v1"


def grid(source):
    transform = source.GetGeoTransform(can_return_null=True)
    if transform is None:
        raise ValueError("Expected a GeoTIFF fitted to the MMapper room grid")
    west, pixel_x, row_rotation, north, column_rotation, pixel_y = transform
    if pixel_x <= 0 or not math.isclose(pixel_y, -pixel_x, abs_tol=1e-12):
        raise ValueError("Expected square north-up pixels in MMapper room units")
    if abs(row_rotation) > 1e-12 or abs(column_rotation) > 1e-12:
        raise ValueError("Align the raster to the MMapper axes before baking")
    east = west + pixel_x * source.RasterXSize
    south = north + pixel_y * source.RasterYSize
    if any(abs(edge - round(edge)) > 1e-6 for edge in (west, south, east, north)):
        raise ValueError("GeoTIFF edges must fall on whole room coordinates")
    return round(west), round(north), round(east - west), round(north - south), pixel_x


def make_signed_alpha(band, pixel):
    elevations = band.ReadAsArray(buf_type=gdal.GDT_Float32)
    elevations /= 65535
    blurred = gaussian_filter(elevations, 0.55 / (8 * pixel))
    del elevations
    dy, dx = np.gradient(blurred, pixel, pixel)
    del blurred
    dx *= -8
    dy *= 8
    light = ((-dx + dy + 1.7) / np.sqrt(dx * dx + dy * dy + 1) - 1.7) * (0.85 / np.sqrt(4.89))
    del dx, dy
    shadows = light < 0
    light[shadows] = 0.55 * np.expm1(light[shadows] / 0.55)
    np.clip(light * 255, -127, 127, out=light)
    return np.rint(light, out=light)


def average(alpha, width, height):
    if alpha.shape == (height, width):
        return alpha
    rows, cols = alpha.shape
    source = gdal.GetDriverByName("MEM").Create("", cols, rows, 1, gdal.GDT_Float32)
    source.GetRasterBand(1).WriteArray(alpha)
    output = gdal.Translate("", source, format="MEM", width=width, height=height,
                            resampleAlg="average", outputType=gdal.GDT_Float32)
    return np.rint(output.ReadAsArray())


def bake(dem, output, density):
    if output.suffix.lower() != ".png":
        raise ValueError("The output name must end in .png")
    source = gdal.Open(str(dem))
    west, north, rooms_x, rooms_y, pixel = grid(source)
    band = source.GetRasterBand(1)
    if band.DataType != gdal.GDT_UInt16 or band.GetNoDataValue() is not None:
        raise ValueError("Expected a complete UInt16 Arda elevation band with no declared NoData")
    width, height = rooms_x * density, rooms_y * density
    if min(width, height) < 1 or max(width, height) > MAX_SIDE or width * height > MAX_PIXELS:
        raise ValueError("Choose fewer pixels per room to fit MMapper's PNG limits")
    if min(source.RasterXSize, source.RasterYSize) < 2:
        raise ValueError("At least two source pixels are required on each axis")

    alpha = average(make_signed_alpha(band, pixel), width, height)
    packed = np.empty((height, width, 2), dtype=np.uint8)
    packed[..., 0] = np.where(alpha >= 0, 255, 0)
    packed[..., 1] = np.abs(alpha).astype(np.uint8)
    step = 1 / density
    transform = (west, step, 0, north, 0, -step)
    metadata = PngImagePlugin.PngInfo()
    for key, value in (
        ("Source", source.GetMetadataItem("SOURCE_URL") or SOURCE_URL),
        ("Attribution", "32k heightmap by Pete Lewis (therealJiko), based on the ME-DEM by monks and "
                        "Redrobes. Adapted by aligning the elevation raster to MUME and baking "
                        "transparent hillshade."),
        ("Source license", "https://creativecommons.org/licenses/by/4.0/"),
        ("Title", f"Arda 32k hillshade aligned to MUME, {density} pixels per room"),
        ("Description", "Transparent black shadows and white highlights. North-up. "
                        f"Resampled at {density} pixels per room along each axis."),
        ("MUME geotransform", ",".join(format(value, ".15g") for value in transform)),
    ):
        metadata.add_text(key, value)
    Image.fromarray(packed).save(output, pnginfo=metadata)
    world = (step, 0, 0, -step, west + step / 2, north - step / 2)
    output.with_suffix(".pgw").write_text("".join(f"{value:.15g}\n" for value in world), encoding="ascii")
    print(f"{output}: {width} × {height}, {density} pixels per room")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("aligned_arda_32k_tif", type=Path)
    parser.add_argument("output_png", type=Path)
    parser.add_argument("--pixels-per-room", type=int, default=11)
    args = parser.parse_args()
    if args.pixels_per_room <= 0:
        parser.error("pixels per room must be positive")
    gdal.UseExceptions()
    try:
        bake(args.aligned_arda_32k_tif, args.output_png, args.pixels_per_room)
    except (OSError, ValueError, RuntimeError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
