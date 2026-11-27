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

## Add a map release if redistribution is authorised

The offline `maps/` directory has the exact foreground files and matching background MM2s, hillshades, world files, reports and previews. The MUME bundled-map repository grants no broader permission over its game data, and the older ME-GIS vectors have uncertain licence terms. See [ACKNOWLEDGEMENTS.md](ACKNOWLEDGEMENTS.md). The source manifest lets users fetch the exact foreground maps directly from their existing upstream hosts.

With any necessary permissions in place, attach the map ZIP to the repository’s first release. Keep the acknowledgements with it.

```sh
python3 tools/fetch_foregrounds.py --verify-only &&
mkdir -p dist &&
zip -qr dist/MMapper-Terrain-Maps-v0.1.0.zip maps README.md ACKNOWLEDGEMENTS.md BUILD.md foreground-sources.json -x '*/.DS_Store' &&
gh release create v0.1.0 dist/MMapper-Terrain-Maps-v0.1.0.zip \
  --title "MMapper Terrain Helper v0.1.0" --generate-notes
```
