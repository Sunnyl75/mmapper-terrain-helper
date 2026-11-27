# Publish on GitHub

Extract the starter archive into a fresh directory. From its root, check the staged paths. The default commit includes the helper, exact source links, profile fits and acknowledgements. The `maps/` directory is kept offline by Git.

```sh
git init -b main
git add -- .gitignore .gitattributes README.md ACKNOWLEDGEMENTS.md BUILD.md PUBLISH.md \
  foreground-sources.json tools workflow profiles fit
git status --short
git commit -m "Add MMapper terrain helper and verified map profiles"

gh auth status || gh auth login --hostname github.com --git-protocol https --web
gh auth setup-git
gh repo create Sunnyl75/mmapper-terrain-helper --public --source=. --remote=origin --push
```

The command uses GitHub CLI. On macOS with Homebrew, install it with `brew install gh` when needed.

## Publish the finished backgrounds and DEMs

Use the extracted starter bundle. It has the exact generated backgrounds and hillshades under `maps/profiles/`. Place `MMapper-River-Reach-Elevation-and-Vectors.zip` in Downloads. The helper extracts the aligned elevation GeoTIFF from that archive for each profile.

```sh
cd "$HOME/Downloads/MMapper-Terrain-Helper" &&
git pull --ff-only &&
bash tools/publish_outputs.sh "$HOME/Downloads/MMapper-River-Reach-Elevation-and-Vectors.zip"
```

The script creates a GitHub Release with three background ZIPs and three elevation GeoTIFFs. Each background ZIP contains the terrain-only MM2, matching hillshade and PGW, build report, and acknowledgements. The release notes point to the upstream foreground pages. The [source manifest](foreground-sources.json) records the exact playable snapshots matched to the release. See [acknowledgements and data terms](ACKNOWLEDGEMENTS.md).

The default release tag is `v0.1.0`. For a future release, update the output files and pass the new tag as the second argument.
