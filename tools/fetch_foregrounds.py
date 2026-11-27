"""Fetch exact public upstream foreground snapshots or verify offline copies."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
from tempfile import NamedTemporaryFile
from urllib.error import URLError
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(entry, path):
    temp = None
    try:
        with NamedTemporaryFile(dir=path.parent, suffix=".download", delete=False) as file:
            temp = Path(file.name)
            with urlopen(entry["url"], timeout=60) as source:
                shutil.copyfileobj(source, file)
        if sha256(temp) != entry["sha256"]:
            raise ValueError(f"Upstream download failed SHA-256 verification for {path.name}")
        temp.replace(path)
    finally:
        if temp is not None:
            temp.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true", help="Use local files and make no network requests")
    parser.add_argument("--output", type=Path, default=ROOT / "maps" / "foreground")
    args = parser.parse_args()
    records = json.loads((ROOT / "foreground-sources.json").read_text(encoding="utf-8"))
    if not args.verify_only:
        args.output.mkdir(parents=True, exist_ok=True)
    try:
        for profile, entry in records.items():
            path = args.output / entry["filename"]
            if not path.exists():
                if args.verify_only:
                    raise ValueError(f"Missing {path}")
                fetch(entry, path)
            if sha256(path) != entry["sha256"]:
                raise ValueError(f"SHA-256 mismatch at {path}. Use an empty output directory to fetch a clean copy.")
            print(f"{profile}: {path.name} verified")
    except (OSError, URLError, ValueError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
