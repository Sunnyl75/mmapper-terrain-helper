# Build a matching background

Use Python 3.11 or later with GDAL, NumPy, SciPy and Pillow. A conda environment works on macOS and Linux.

```sh
conda create -n mmapper-gis -c conda-forge python=3.12 gdal numpy scipy pillow
conda activate mmapper-gis
```

Get the pinned Arda checkout and place the four 32k GeoTIFF quadrants next to its `data/rasters/32k/dem.vrt`. The helper uses `data/vectors/vectors.gpkg` from the same checkout. The source inputs stay local. [Pinned Arda checkout](https://github.com/bburns/Arda/tree/295d735320affed232a0eadffeafb0482b660971) · [Download the 32k quadrants](https://github.com/bburns/Arda/releases/tag/dem-32k-v1) · [Source raster checksums](workflow/arda-32k-sha256.txt)

Get the exact playable maps with `python3 tools/fetch_foregrounds.py`. If you already have the offline bundle, use `python3 tools/fetch_foregrounds.py --verify-only`. Once the dependencies and source files are present, the following commands run offline.

## Build one of the included snapshots

From the repository directory, set the location of your Arda checkout and choose the profile.

```sh
arda="$PWD/local-sources/Arda"
profile=faine
reference="maps/foreground/arda-MMapper-26.06.0.xml"

case "$profile" in
  bundled) playable="$reference" ;;
  faine)   playable="maps/foreground/faine-map.mm2" ;;
  nazgum)  playable="maps/foreground/nazgum-latest.mm2" ;;
esac

profile_fit="profiles/$profile/regional-displacement.npz"
output="build/$profile"

python workflow/rebuild.py \
  --dem "$arda/data/rasters/32k/dem.vrt" \
  --gis "$arda/data/vectors/vectors.gpkg" \
  --reference "$reference" --target "$playable" \
  --profile-fit "$profile_fit" --reach-fit fit --output "$output"

python workflow/review.py --target "$playable" \
  --profile-fit "$profile_fit" --output "$output"
```

The output contains the background `.mm2`, transparent PNG and its PGW, elevation GeoTIFF, aligned GeoPackage, and a local HTML review. The export keeps exact public surface sectors and directional main-road flags. Use `--exclude-outdoor path/to/exclusions.csv` with both commands when applying a reviewed outdoor secret list. Copy [the template](workflow/secret-outdoor-template.csv) and record each excluded room by its target x and y coordinates.

## Fit a newer playable snapshot

Choose the profile name and set `playable` to the new file. Use the same Bundled reference for these geographic controls. Generate a new fit and set `profile_fit` to the resulting NPZ before running both commands above.

```sh
python workflow/bake_profile.py --audit-only --bounds -6 -271 721 21 \
  --profile "$profile" --reference "$reference" --target "$playable" \
  --output "build/$profile-fit"

profile_fit="build/$profile-fit/regional-displacement.npz"
```

The helper reads MMapper XML and MM2 schemas 36–42. Export other foreground formats as MMapper XML. Use validated room pairs for a materially shifted area. The documented room-pair format is in [reviewed-pairs-template.csv](workflow/reviewed-pairs-template.csv) and is supplied to the fit with `--controls`. The default footprint is 727 × 292 rooms. Extra playable rooms outside it are reported.

## Refit the reference rivers

The packaged `fit/reach-fit.json` records accepted and unresolved reaches, survey stations and withheld checks. To recompute the river pass, build a Bundled output to `build/bundled-base` with `--reach-fit` omitted. Then run:

```sh
python workflow/river_reaches.py --map "$reference" \
  --dem build/bundled-base/elevation.tif \
  --gis build/bundled-base/aligned-gis.gpkg \
  --reference-fit build/bundled-base/reference-pull.npz \
  --output build/new-reach-fit
```

Use `--reach-fit build/new-reach-fit` for the finished build. The helper uses original elevation and vectors, with local corrections for supported river reaches and exactly overlaid public room sectors. See the report for remaining uncertain waterways.
