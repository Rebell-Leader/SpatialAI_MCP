#!/usr/bin/env python3
"""Flag drift between our context/ assumptions and the live OpenProblems repo.

Our `context/ist-preprocessing-pipeline.md` documents the pipeline stages of
`openproblems-bio/task_ist_preprocessing`. That repo moves; stale guidance is
worse than none. This script fetches the upstream `src/` directory listing and
compares it against the stage directories we expect, printing anything added or
removed so we can refresh the docs.

By default it is *report-only* (exit 0) — including when GitHub is unreachable —
so it never blocks CI on network flakiness. Pass --strict to exit non-zero when
drift (or a fetch error) is detected, e.g. for a scheduled alert.

Usage:
    python scripts/check_upstream_sync.py [--strict]
"""

from __future__ import annotations

import argparse
import json
import os
import urllib.error
import urllib.request

UPSTREAM = "openproblems-bio/task_ist_preprocessing"
API_URL = f"https://api.github.com/repos/{UPSTREAM}/contents/src"

# Top-level directories we currently assume exist under the upstream `src/`.
# Keep this in sync with context/ist-preprocessing-pipeline.md. If this script
# reports additions/removals, update both.
EXPECTED_DIRS = {
    "api",
    "base",
    "control_methods",
    "core",
    "data_processors",
    "datasets",
    "methods_calculate_cell_volume",
    "methods_cell_type_annotation",
    "methods_count_aggregation",
    "methods_data_aggregation",
    "methods_expression_correction",
    "methods_normalization",
    "methods_qc_filter",
    "methods_segmentation",
    "methods_transcript_assignment",
    "metrics",
    "workflows",
}


def fetch_upstream_dirs() -> set[str]:
    """Return the set of directory names under the upstream src/."""
    headers = {
        "Accept": "application/vnd.github+json",
        "User-Agent": "spatialai-upstream-sync-check",
    }
    # Use a token when available (CI) to avoid unauthenticated rate limits.
    token = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    if token:
        headers["Authorization"] = f"Bearer {token}"
    req = urllib.request.Request(API_URL, headers=headers)
    with urllib.request.urlopen(req, timeout=20) as resp:
        entries = json.load(resp)
    return {e["name"] for e in entries if e.get("type") == "dir"}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit non-zero on drift or fetch error (default: report-only).",
    )
    args = parser.parse_args(argv)

    print(f"Checking upstream layout: {UPSTREAM}/src\n")

    try:
        actual = fetch_upstream_dirs()
    except (urllib.error.URLError, TimeoutError, ValueError) as exc:
        print(f"⚠️  Could not fetch upstream ({exc}). Skipping drift check.")
        return 1 if args.strict else 0

    added = sorted(actual - EXPECTED_DIRS)
    removed = sorted(EXPECTED_DIRS - actual)

    if not added and not removed:
        print(f"✅ In sync. {len(actual)} stage directories match expectations.")
        return 0

    print("⚠️  Drift detected between upstream src/ and our assumptions:\n")
    if added:
        print("  NEW upstream (document these in context/):")
        for name in added:
            print(f"    + {name}")
    if removed:
        print("\n  GONE from upstream (our context/ may reference stale stages):")
        for name in removed:
            print(f"    - {name}")
    print(
        "\nUpdate context/ist-preprocessing-pipeline.md and EXPECTED_DIRS in "
        "this script, then re-run."
    )
    return 1 if args.strict else 0


if __name__ == "__main__":
    raise SystemExit(main())
