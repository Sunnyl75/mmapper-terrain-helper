#!/usr/bin/env bash
set -euo pipefail

elevation_zip="${1:-$HOME/Downloads/MMapper-River-Reach-Elevation-and-Vectors.zip}"
tag="${2:-v0.1.0}"
repo="Sunnyl75/mmapper-terrain-helper"
case "$elevation_zip" in
    /*) ;;
    *) elevation_zip="$PWD/$elevation_zip" ;;
esac
cd "$(dirname "$0")/.."

if [[ ! -r "$elevation_zip" ]]; then
    printf 'Cannot read elevation archive: %s\n' "$elevation_zip" >&2
    exit 1
fi

out="$(mktemp -d "${TMPDIR:-/tmp}/mmapper-terrain.XXXXXXXX")"
printf 'Release assets are staged in %s\n' "$out"
assets=()
for profile in bundled faine nazgum; do
    src="maps/profiles/$profile"
    background="$out/$profile-background.zip"
    elevation="$out/$profile-elevation.tif"

    profile_files=(
        "$src/$profile-background-river-reaches.mm2"
        "$src/$profile.hillshade.png"
        "$src/$profile.hillshade.pgw"
        "$src/build-report.json"
        ACKNOWLEDGEMENTS.md
    )
    for file in "${profile_files[@]}"; do
        if [[ ! -s "$file" ]]; then
            printf 'Missing release input: %s\n' "$file" >&2
            exit 1
        fi
    done
    zip -q -j "$background" "${profile_files[@]}"
    unzip -p "$elevation_zip" "$profile/elevation.tif" > "$elevation"
    [[ -s "$elevation" ]]
    assets+=("$background" "$elevation")
done

gh release create "$tag" "${assets[@]}" \
    --repo "$repo" --target main --title "MMapper terrain profiles $tag" \
    --notes-file - <<'NOTES'
Choose the Bundled, Faine, or Nazgum background ZIP. It contains a terrain-only .mm2, its matching hillshade PNG and PGW, and the build report. Use the separate elevation GeoTIFF in QGIS.

| Foreground | Where to look for updates |
| --- | --- |
| Bundled | [MUME bundled map](https://github.com/MUME/arda) · [MMapper releases](https://github.com/MUME/MMapper/releases) |
| Faine | [Faine's map in Nazgum's map directory](https://github.com/nazgum/mume/tree/master/maps) |
| Nazgum | [Nazgum's map directory](https://github.com/nazgum/mume/tree/master/maps) |

Use the exact [foreground snapshots](https://github.com/Sunnyl75/mmapper-terrain-helper/blob/main/foreground-sources.json) with these downloads. The source pages above may carry newer maps that need a new alignment.

Elevation by Pete Lewis (therealJiko), adapted from the Arda 32k heightmap under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), based on ME-DEM by monks and Redrobes. Backgrounds use ME-GIS and public MUME terrain. [Source credits and data terms](https://github.com/Sunnyl75/mmapper-terrain-helper/blob/main/ACKNOWLEDGEMENTS.md).
NOTES
