#!/usr/bin/env python3
"""Score a produced Viash component against the case-study rubric.

Grades *static* conformance to OpenProblems / Viash / SpatialData conventions
offline (0-11). Build and test (`viash ns build`, `viash test`) need a real
Viash/Docker runner and are reported as runtime checks, optionally executed with
--run-viash.

Usage:
    python case-study/grade.py path/to/component_dir        # human summary
    python case-study/grade.py path/to/config.vsh.yaml --json
    python case-study/grade.py path/to/component_dir --run-viash <task_repo>
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

import yaml

PASS_THRESHOLD = 9  # static points required (with build+test) for task success


@dataclass
class Check:
    key: str
    points: int
    passed: bool
    detail: str


@dataclass
class Report:
    config_path: str
    checks: list[Check] = field(default_factory=list)
    runtime: dict = field(default_factory=dict)

    @property
    def score(self) -> int:
        return sum(c.points for c in self.checks if c.passed)

    @property
    def max_score(self) -> int:
        return sum(c.points for c in self.checks)

    @property
    def static_pass(self) -> bool:
        return self.score >= PASS_THRESHOLD


def _find_config(target: Path) -> Path:
    if target.is_file():
        return target
    candidates = sorted(target.rglob("config.vsh.yaml"))
    if not candidates:
        raise FileNotFoundError(f"No config.vsh.yaml found under {target}")
    return candidates[0]


def _script_text(config: dict, config_path: Path) -> str:
    """Concatenate the text of all script resources referenced by the config.

    Looks in both the modern top-level ``resources`` and the legacy
    ``functionality.resources`` so script-based checks are evaluated honestly
    even when the config uses the obsolete layout.
    """
    resources = list(config.get("resources") or [])
    legacy = config.get("functionality")
    if isinstance(legacy, dict):
        resources += list(legacy.get("resources") or [])
    texts = []
    for res in resources:
        if (
            isinstance(res, dict)
            and "path" in res
            and "script" in str(res.get("type", ""))
        ):
            p = (config_path.parent / res["path"]).resolve()
            if p.is_file():
                texts.append(p.read_text(encoding="utf-8", errors="ignore"))
    return "\n".join(texts)


def _engines(config: dict) -> list[dict]:
    return [e for e in (config.get("engines") or []) if isinstance(e, dict)]


def _merge_values(config: dict) -> str:
    m = config.get("__merge__")
    if isinstance(m, str):
        return m
    if isinstance(m, list):
        return " ".join(str(x) for x in m)
    return ""


def grade(config_path: Path) -> Report:
    raw = config_path.read_text(encoding="utf-8", errors="ignore")
    config = yaml.safe_load(raw) or {}
    script = _script_text(config, config_path)
    engines = _engines(config)
    images = [str(e.get("image", "")) for e in engines if e.get("type") == "docker"]
    report = Report(config_path=str(config_path))

    def add(key, points, passed, detail):
        report.checks.append(Check(key, points, bool(passed), detail))

    # 1. Modern flat config (no legacy functionality:/platforms:)
    legacy = [k for k in ("functionality", "platforms") if k in config]
    add(
        "modern_flat_config",
        1,
        not legacy,
        "flat config" if not legacy else f"legacy keys present: {legacy}",
    )

    # 2. Merges the stage component API
    merge = _merge_values(config)
    has_api = "/src/api/comp_" in merge or bool(re.search(r"comp_\w*", merge))
    add(
        "merges_comp_api",
        1,
        has_api,
        f"__merge__={merge!r}" if merge else "no __merge__ of comp_ API",
    )

    # 3. Builds on a shared OpenProblems base image (not a raw python/r-base image)
    on_base = any(img.startswith("openproblems/base_") for img in images)
    add(
        "openproblems_base_image",
        1,
        on_base,
        f"images={images}" if images else "no docker engine image",
    )

    # 4. Core metadata block
    core = all(
        str(config.get(k, "")).strip() for k in ("label", "summary", "description")
    )
    add(
        "metadata_core",
        1,
        core,
        (
            "label+summary+description present"
            if core
            else "missing label/summary/description"
        ),
    )

    # 5. links + references
    has_links_refs = bool(config.get("links")) and bool(config.get("references"))
    add(
        "metadata_links_references",
        1,
        has_links_refs,
        "links+references present" if has_links_refs else "missing links/references",
    )

    # 6. Pinned (no :latest, base image carries a tag)
    no_latest = ":latest" not in raw
    tagged = all(":" in img for img in images) if images else False
    add(
        "pinned_versions",
        1,
        no_latest and tagged,
        "pinned" if (no_latest and tagged) else "unpinned image or :latest present",
    )

    # 7. Nextflow runner with resource directives
    runners = [r for r in (config.get("runners") or []) if isinstance(r, dict)]
    nf = next((r for r in runners if r.get("type") == "nextflow"), None)
    has_directives = bool(nf and nf.get("directives", {}).get("label"))
    add(
        "nextflow_directives",
        1,
        has_directives,
        (
            "nextflow directives.label set"
            if has_directives
            else "no nextflow resource label"
        ),
    )

    # 8. VIASH START/END block in the script
    has_viash_block = "## VIASH START" in script and "## VIASH END" in script
    add(
        "viash_param_block",
        1,
        has_viash_block,
        "VIASH START/END present" if has_viash_block else "no VIASH START/END block",
    )

    # 9. SpatialData I/O (not anndata-only)
    uses_sd = bool(
        re.search(
            r"\bimport spatialdata\b|spatialdata as sd|read_zarr|\.tables\b", script
        )
    )
    anndata_only = bool(re.search(r"read_h5ad|import anndata", script)) and not uses_sd
    add(
        "spatialdata_io",
        1,
        uses_sd,
        (
            "uses SpatialData I/O"
            if uses_sd
            else (
                "anndata-only / no zarr I/O"
                if anndata_only
                else "no recognizable spatial I/O"
            )
        ),
    )

    # 10. Fails loudly on error
    fails_loud = bool(re.search(r"sys\.exit\(\s*[1-9]|raise\s+\w+", script))
    add(
        "fails_loudly",
        1,
        fails_loud,
        "non-zero exit / raise present" if fails_loud else "no loud failure path",
    )

    # 11. Ships a unit test
    test_res = config.get("test_resources") or []
    has_test = any(
        isinstance(r, dict) and "script" in str(r.get("type", "")) for r in test_res
    )
    add(
        "has_unit_test",
        1,
        has_test,
        "test_resources script present" if has_test else "no viash test",
    )

    return report


def maybe_run_viash(report: Report, component_dir: Path, task_repo: Path) -> None:
    """Optionally run `viash ns build` / `viash test` if viash is available."""
    if not shutil.which("viash"):
        report.runtime = {"viash_available": False, "note": "viash not on PATH"}
        return
    cfg = _find_config(component_dir)
    results = {"viash_available": True}
    for name, cmd in (
        ("build", ["viash", "ns", "build"]),
        ("test", ["viash", "test", str(cfg)]),
    ):
        try:
            proc = subprocess.run(
                cmd, cwd=task_repo, capture_output=True, text=True, timeout=1800
            )
            results[name] = proc.returncode == 0
        except (subprocess.SubprocessError, OSError) as exc:
            results[name] = False
            results[f"{name}_error"] = str(exc)
    report.runtime = results


def render(report: Report) -> str:
    lines = [f"Component: {report.config_path}", ""]
    for c in report.checks:
        mark = "✅" if c.passed else "❌"
        lines.append(f"  {mark} [{c.points}] {c.key}: {c.detail}")
    lines.append("")
    lines.append(f"Static conformance: {report.score}/{report.max_score}")
    lines.append(
        f"Static gate (≥{PASS_THRESHOLD}): {'PASS' if report.static_pass else 'FAIL'}"
    )
    if report.runtime:
        lines.append(f"Runtime: {report.runtime}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", type=Path, help="Component dir or config.vsh.yaml")
    parser.add_argument("--json", action="store_true", help="Emit JSON")
    parser.add_argument(
        "--run-viash",
        type=Path,
        metavar="TASK_REPO",
        help="Also run viash ns build / viash test from this task repo root.",
    )
    args = parser.parse_args(argv)

    try:
        config_path = _find_config(args.target)
    except FileNotFoundError as exc:
        print(f"❌ {exc}", file=sys.stderr)
        return 2

    report = grade(config_path)
    if args.run_viash:
        maybe_run_viash(report, args.target, args.run_viash)

    if args.json:
        print(
            json.dumps(
                {
                    "config_path": report.config_path,
                    "score": report.score,
                    "max_score": report.max_score,
                    "static_pass": report.static_pass,
                    "checks": [vars(c) for c in report.checks],
                    "runtime": report.runtime,
                },
                indent=2,
            )
        )
    else:
        print(render(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
