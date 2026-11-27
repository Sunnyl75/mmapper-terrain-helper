# Scheduled terrain updates

Open **Actions → Check and rebuild terrain maps**. The workflow checks the Bundled, Faine and Nazgum sources each Monday and can also run on demand. It compares SHA-256 of the actual foreground with the original supported snapshots and skips update drafts or releases it has already built. It follows an explicit upstream Git rename. A source with an ambiguous move is reported in the run summary.

## Use a different playable map

1. Open **Run workflow** and paste the public **Raw** GitHub URL of the map into **map_url**. Supported inputs are MMapper XML and MM2 schemas 36–42, up to 70 MiB. Use a commit-pinned URL if the source ref contains a slash.
2. Leave **style** on **auto**. The helper compares room matches, their shared cardinal connections and their coordinate offsets across all three presets. An exact snapshot reuses its verified fit. For a changed map it refines against the closest preset and checks the composed transform against a direct fit to the Bundled reference. You can choose a preset yourself if you know the source.
3. Click **Run workflow**. For an upstream check alone, leave **map_url** blank and select **check_only**.

A verified build creates a draft release for that precise foreground and links to the original source. The draft contains the terrain-only background, matching hillshade and PGW, aligned elevation, and provenance. The run summary or artifacts show the selected style and any fit concerns. The job checks the MM2 read-back, foreground terrain and road flags, known river continuity, and fit residuals. A failed check leaves the existing releases available and the report in the workflow run. Preview the result in MMapper before publishing the draft.

The scanner needs shared outdoor rooms to infer a layout. It stops with an explanation when none of the presets has enough reliable correspondences. The existing manual control workflow in [BUILD.md](BUILD.md) remains available for maps with major changes or detached geometry.

## Fit a local file

With the dependencies from [BUILD.md](BUILD.md) installed, run the following from the repository directory. Keep private maps local.

```sh
python3 tools/fetch_foregrounds.py
python3 tools/terrain_updates.py fit \
  --target "/path/to/playable.mm2" --output build/custom-fit --style auto
```

Then use *build/custom-fit/regional-displacement.npz* for the rebuild and verification in BUILD.md. Set a different output directory for a second run. The reports record the style selection and any fallback.

GitHub disables scheduled workflows in inactive public repositories after 60 days. You can re-enable the workflow from its Actions page.
