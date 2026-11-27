# MMapper Terrain Helper

Build matching outdoor backgrounds and hillshades for MUME’s [MMapper](https://github.com/MUME/MMapper). The helper aligns Arda elevation and GIS layers to a playable map and fits rivers along supported reaches. It exports a compact terrain-only `.mm2` background and a hillshade with a world file.

## Use a ready-made profile

Check [Releases](https://github.com/Sunnyl75/mmapper-terrain-helper/releases) for versioned backgrounds, hillshades, and aligned elevation GeoTIFFs. The hillshade PNG and its PGW go with the matching background. Use the elevation GeoTIFF in GIS software.

| Profile | Foreground used for this build | Page to check for updates |
| --- | --- | --- |
| Bundled | `arda-MMapper-26.06.0.xml` | [MUME's bundled map repository](https://github.com/MUME/arda) · [MMapper releases](https://github.com/MUME/MMapper/releases) |
| Faine | `faine-map.mm2` | [Faine's map in Nazgum's map directory](https://github.com/nazgum/mume/tree/master/maps) |
| Nazgum | `nazgum-latest.mm2` | [Nazgum's map directory](https://github.com/nazgum/mume/tree/master/maps) |

These are upstream pages that can show renamed or updated files. For an existing background, use the exact playable snapshot recorded in `foreground-sources.json`.

1. Open the matching playable map in a build with **Background Map** and **Hillshade image** preferences.
2. Enable **Show background map** and choose the profile’s background `.mm2`.
3. Choose the profile’s `*.hillshade.png`. Keep its `.pgw` beside it and set both offsets to zero. Open `review.html` beside the background to see comparisons.

Fetch the exact playable maps from their original hosts with `python3 tools/fetch_foregrounds.py`. Add `--verify-only` to check local copies offline. The pinned commits and SHA-256 digests are in [foreground-sources.json](foreground-sources.json). Use the corresponding fit for each map snapshot. Faine’s archived copy has no server IDs, so Fog of War room reveal requires a Faine foreground with server IDs.

## Scheduled updates and other playable maps

See [AUTOMATION.md](AUTOMATION.md) for weekly upstream checks, draft releases, and running a supplied map through the closest starting profile.

## Rebuild a profile

Follow [BUILD.md](BUILD.md) to install the Python dependencies, supply the local Arda inputs and generate the background, aligned DEM, GIS and hillshade. Maps remain flat on screen. The background keeps public outdoor terrain and main roads. Supply a reviewed exclusion list if an outdoor area needs to remain concealed for Fog of War.

See [ACKNOWLEDGEMENTS.md](ACKNOWLEDGEMENTS.md) for source credits and data terms. This is a fan project.
