#!/usr/bin/env python3
"""Discover foreground updates and fit a supplied map to a reference profile."""

import argparse
from collections import defaultdict
from hashlib import file_digest, sha256
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from urllib.error import HTTPError
from urllib.parse import quote, unquote, urlsplit
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
PROFILES = ("bundled", "faine", "nazgum")
SOURCES = {
    "bundled": ("MUME/arda", "arda.xml"),
    "faine": ("nazgum/mume", "maps/faine-map.mm2"),
    "nazgum": ("nazgum/mume", "maps/nazgum-latest.mm2"),
}
BOUNDS = (-6, -271, 721, 21)
MAP_LIMIT = 70 * 1024 * 1024


def json_file(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def api(path, accept="application/vnd.github+json"):
    headers = {"Accept": accept, "X-GitHub-Api-Version": "2022-11-28"}
    token = os.environ.get("GH_TOKEN")
    if token and path.startswith(f"repos/{os.environ.get('GITHUB_REPOSITORY')}/"):
        headers["Authorization"] = f"Bearer {token}"
    with urlopen(Request("https://api.github.com/" + path, headers=headers), timeout=30) as stream:
        return json.load(stream)


def raw_url(repo, revision, path):
    return f"https://raw.githubusercontent.com/{repo}/{revision}/{quote(path, safe='/')}"


def map_bytes(url):
    if urlsplit(url).scheme != "https" or urlsplit(url).netloc != "raw.githubusercontent.com":
        raise ValueError("Use an HTTPS Raw GitHub URL")
    with urlopen(url, timeout=75) as stream:
        data = stream.read(MAP_LIMIT + 1)
    if not data or len(data) > MAP_LIMIT:
        raise ValueError("Supply a nonempty MMapper map under 70 MiB")
    return data


def tracked_path(repo, old_path, head, start):
    def exists(path):
        try:
            info = api(f"repos/{repo}/contents/{quote(path, safe='/')}?ref={head}",
                       accept="application/vnd.github.object+json")
        except HTTPError as error:
            if error.code == 404:
                return False
            raise
        return info.get("type") == "file"

    if exists(old_path):
        return old_path
    name = old_path
    # Compare against the original snapshot and accept explicit Git rename
    # records only. A removed path without a clear rename needs a new input.
    files = api(f"repos/{repo}/compare/{start}...{head}").get("files", [])
    candidates = [f["filename"] for f in files if f.get("status") == "renamed"
                  and f.get("previous_filename") == name]
    if len(candidates) == 1 and exists(candidates[0]):
        return candidates[0]
    raise ValueError(f"{repo}/{old_path} moved without a unique rename. Supply its new Raw URL manually.")


def parse_manual(url):
    parsed = urlsplit(url)
    if parsed.scheme != "https" or parsed.netloc != "raw.githubusercontent.com" or parsed.query or parsed.fragment:
        raise ValueError("Supply the public Raw GitHub URL for the map")
    match = re.fullmatch(r"/([A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+)/([^/]+)/(.+)", parsed.path)
    if not match:
        raise ValueError("The Raw URL must identify a repository, commit or branch, and file path")
    repo, ref, path = match.groups()
    path = unquote(path)
    if any(p in ("", ".", "..") for p in path.split("/")) or Path(path).suffix.lower() not in (
            ".mm2", ".xml", ".mm2xml"):
        raise ValueError("Choose a Raw .mm2, .xml, or .mm2xml map file")
    return repo, unquote(ref), path


def find_candidates(args):
    own_repo = os.environ.get("GITHUB_REPOSITORY", "Sunnyl75/mmapper-terrain-helper")
    baseline = json_file(ROOT / "foreground-sources.json")
    existing = set()
    for page in range(1, 11):
        found = api(f"repos/{own_repo}/releases?per_page=100&page={page}")
        existing.update(release["tag_name"] for release in found)
        if len(found) < 100:
            break

    heads = {}
    latest = []
    blocked = []
    if args.map_url:
        repo, ref, path = parse_manual(args.map_url)
        pending = [("custom", repo, ref, path)]
    else:
        pending = [(p, repo, "HEAD", path) for p, (repo, path) in SOURCES.items()]

    for label, repo, ref, path in pending:
        try:
            key = (repo, ref)
            if key not in heads:
                branch = api(f"repos/{repo}")["default_branch"] if ref == "HEAD" else ref
                heads[key] = api(f"repos/{repo}/commits/{quote(branch, safe='')}")["sha"]
            head = heads[key]
            if label in baseline:
                start = urlsplit(baseline[label]["url"]).path.split("/")[3]
                path = tracked_path(repo, path, head, start)
            url = raw_url(repo, head, path)
            digest = sha256(map_bytes(url)).hexdigest()
            if label in baseline and digest == baseline[label]["sha256"]:
                continue
            suffix = Path(path).suffix.lower()
            if suffix not in (".mm2", ".xml", ".mm2xml"):
                raise ValueError(f"The new path {path} is not a supported map")
            tag = f"terrain-{label}-{digest[:16]}"
            if tag not in existing:
                latest.append({"label": label, "tag": tag, "url": url, "sha256": digest, "extension": suffix,
                               "source": f"https://github.com/{repo}/blob/{head}/{quote(path, safe='/')}"})
        except ValueError as error:
            if args.map_url:
                raise
            blocked.append({"profile": label, "reason": str(error)})
    plan = {"include": latest, "review": blocked}
    write_json(args.output, plan)
    print(json.dumps(plan, separators=(",", ":")))


def fetch_candidate(args):
    data = map_bytes(args.url)
    if sha256(data).hexdigest() != args.sha256:
        raise ValueError("Candidate bytes differ from the recorded map")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(data)


def run_workflow(module, *args):
    subprocess.run([sys.executable, str(ROOT / "workflow" / module), *map(str, args)], check=True)


def choose_style(target, presets, requested):
    import numpy as np
    from register import match, robust_controls

    comparisons = {}
    outdoors = sum(r.outdoor for r in target)
    for name, preset in presets.items():
        pairs, topo = match(preset, target, None)
        if len(pairs) < 30:
            continue
        src, dst, rejected = robust_controls(pairs, topo)
        use = ~rejected
        if use.sum() < 30:
            continue
        src, dst = src[use], dst[use]
        offsets = src - dst
        offset = np.rint(np.median(offsets, axis=0))
        residual = np.linalg.norm(offsets - offset, axis=1)
        # Match topology and coordinates after a global move. Normalize
        # coverage against both map sizes so partial inputs can still match.
        overlap = float(len(dst) / min(outdoors, sum(r.outdoor for r in preset)))
        agreement = float(np.mean(residual < .4))
        incidence = float(np.mean(topo[use] >= 2))
        score = .55 * agreement + .30 * overlap + .15 * incidence
        comparisons[name] = {"score": score, "matched_controls": int(len(dst)),
                             "coverage": overlap, "same_layout_fraction": agreement,
                             "linked_controls_fraction": incidence,
                             "preset_from_target_offset": offset.astype(int).tolist()}
    if not comparisons:
        raise ValueError("No preset has enough shared outdoor controls for a fit")
    order = sorted(comparisons, key=lambda name: -comparisons[name]["score"])
    chosen = order[0] if requested == "auto" else requested
    if chosen not in comparisons:
        raise ValueError(f"The selected {chosen} preset has too few matches")
    gap = comparisons[order[0]]["score"] - comparisons[order[1]]["score"] if len(order) > 1 else 1.
    return chosen, comparisons, gap


def fit_target(args):
    sys.path.insert(0, str(ROOT / "workflow"))
    import numpy as np
    from bake_profile import Pull, digest
    from register import determinant, match, robust_controls, summary
    from room_io import read_map

    sources = json_file(ROOT / "foreground-sources.json")
    paths = {name: args.presets / entry["filename"] for name, entry in sources.items()}
    presets = {name: read_map(path) for name, path in paths.items()}
    target = read_map(args.target)
    target_sha = digest(args.target)
    exact = [name for name, entry in sources.items() if entry["sha256"] == target_sha]
    if exact and args.style in ("auto", exact[0]):
        chosen, comparisons, gap = exact[0], {}, 1.
        if not args.choose_only:
            shutil.copytree(ROOT / "profiles" / chosen, args.output)
        detail = {"style": chosen, "style_margin": gap, "method": "exact snapshot preset"}
    else:
        chosen, comparisons, gap = choose_style(target, presets, args.style)
        detail = {"style": chosen, "style_margin": gap, "candidates": comparisons}
        if not args.choose_only:
            # Fit target-to-Bundled independently as a check, then try the
            # fitted target-to-style layout followed by the saved style warp.
            direct = args.output
            run_workflow("bake_profile.py", "--audit-only", "--bounds", *BOUNDS,
                         "--profile", chosen, "--reference", paths["bundled"],
                         "--target", args.target, "--output", direct)
            prior_path = ROOT / "profiles" / chosen / "regional-displacement.npz"
            fit_report = json_file(direct / "profile-report.json")
            detail["method"] = "direct Bundled fit"
            if chosen != "bundled":
                with np.load(prior_path) as p:
                    prior = Pull(p["offset"], p["bounds"], p["center"])
                bridge = direct.parent / "target-to-style"
                try:
                    run_workflow("bake_profile.py", "--audit-only", "--bounds", *map(int, prior.bounds),
                                 "--profile", chosen, "--reference", paths[chosen],
                                 "--target", args.target, "--output", bridge)
                    with np.load(bridge / "regional-displacement.npz") as p:
                        stage = Pull(p["offset"], p["bounds"], p["center"])
                    centre = stage.center + prior.center
                    west, south, east, north = map(int, stage.bounds)
                    rows = np.linspace(south, north, 2*(north-south)+1)
                    cols = np.linspace(west, east, 2*(east-west)+1)
                    field = np.empty((len(rows), len(cols), 2), dtype=np.float32)
                    for first in range(0, len(rows), 40):
                        x, y = np.meshgrid(cols, rows[first:first+40])
                        xy = np.column_stack((x.ravel(), y.ravel()))
                        d = stage.delta(xy)
                        field[first:first+len(y)] = (d+prior.delta(xy+d)-centre).reshape((*x.shape, 2))
                    jacobian = determinant(field[..., 0], field[..., 1], .5)
                    xmin = np.min(cols[None, :] + centre[0] + field[...,0])
                    xmax = np.max(cols[None, :] + centre[0] + field[...,0])
                    ymin = np.min(rows[:, None] + centre[1] + field[...,1])
                    ymax = np.max(rows[:, None] + centre[1] + field[...,1])
                    paired, topo = match(presets["bundled"], target, None)
                    left, right, rejected = robust_controls(paired, topo)
                    use = ~rejected & (np.linalg.norm(left-right-centre, axis=1) <= 10)
                    use &= (right[:,0]>=west) & (right[:,0]<east) & (right[:,1]>=south) & (right[:,1]<north)
                    if use.sum() < 30:
                        raise ValueError("Too few shared Bundled controls for the style composition")
                    proposal = Pull(field, stage.bounds, centre)
                    with np.load(direct / "regional-displacement.npz") as p:
                        fallback = Pull(p["offset"], p["bounds"], p["center"])
                    composed_error = summary(np.linalg.norm(right[use]+proposal.delta(right[use])-left[use],axis=1))
                    direct_error = summary(np.linalg.norm(right[use]+fallback.delta(right[use])-left[use],axis=1))
                    detail["comparison_on_shared_bundled_controls"] = {
                        "style_composed": composed_error, "direct": direct_error}
                    covered = xmin >= BOUNDS[0]-.001 and xmax <= BOUNDS[2]+.001 and (
                        ymin >= BOUNDS[1]-.001 and ymax <= BOUNDS[3]+.001)
                    if covered and jacobian >= .22 and (composed_error["p95"] <= direct_error["p95"] + .25
                                                         and composed_error["mean"] <= direct_error["mean"] + .15):
                        np.savez_compressed(direct/"regional-displacement.npz", offset=field,
                                            bounds=stage.bounds, center=centre, rate=2)
                        detail["method"] = "target-to-style fit composed with saved preset"
                        fit_report["target_bounds"] = list(map(int,stage.bounds))
                        fit_report["target_to_reference_global_pull_xy"] = centre.tolist()
                        fit_report["final_damped_field_fitted_residual"] = composed_error
                        detail["candidate_minimum_inverse_jacobian"] = jacobian
                    else:
                        detail["style_fallback"] = "Composed field exceeded deformation, coverage, or error limits"
                except (ValueError, subprocess.CalledProcessError) as error:
                    detail["style_fallback"] = str(error)
            fit_report["style_selection"] = detail
            write_json(direct/"profile-report.json",fit_report)
    print(json.dumps(detail, indent=2))
    if args.selection_out:
        write_json(args.selection_out, detail)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    find = modes.add_parser("check")
    find.add_argument("--map-url")
    find.add_argument("--output", type=Path, required=True)
    find.set_defaults(action=find_candidates)
    get = modes.add_parser("fetch")
    get.add_argument("--url", required=True)
    get.add_argument("--sha256", required=True)
    get.add_argument("--output", type=Path, required=True)
    get.set_defaults(action=fetch_candidate)
    fit = modes.add_parser("fit")
    fit.add_argument("--target", type=Path, required=True)
    fit.add_argument("--presets", type=Path, default=ROOT / "maps/foreground")
    fit.add_argument("--style", default="auto", choices=("auto", *PROFILES))
    fit.add_argument("--output", type=Path, required=True)
    fit.add_argument("--selection-out", type=Path)
    fit.add_argument("--choose-only", action="store_true")
    fit.set_defaults(action=fit_target)
    args = parser.parse_args()
    args.action(args)


if __name__ == "__main__":
    main()
