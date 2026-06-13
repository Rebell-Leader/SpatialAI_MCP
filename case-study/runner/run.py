#!/usr/bin/env python3
"""Drive the case-study arms through a coding-agent harness and grade the output.

For each arm it: prepares an isolated copy of the task repo, installs or strips
the skill, runs the agent headless on the frozen prompt, finds the produced
component, and grades it with ``grade.py``. Results are written as JSON and a
markdown table.

This measures the **zero-intervention** scenario (one prompt, no operator help) —
the cleanest automated signal. The "human validations" metric comes from
*interactive* runs the operator does by hand (see case-study/README.md).

Skill arms get the project playbooks injected into the prompt (config
``inject_skills``, default true), so the skill *content* reaches the model
regardless of whether the agent auto-loads skills. Use ``repeats`` / ``--repeats``
(>=3 recommended) to report a median and beat single-run noise.

Examples:
    python case-study/runner/run.py --config case-study/runner/arms.json
    python case-study/runner/run.py --config arms.json --only C_antigravity_skill
    python case-study/runner/run.py --config arms.json --repeats 3
    python case-study/runner/run.py --config arms.json --dry-run

Requires: the chosen harness CLI (e.g. `agy` for Antigravity, `opencode`) on PATH
with credentials configured. Use --dry-run to validate the plan without running.
"""

from __future__ import annotations

import argparse
import json
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASE_STUDY = HERE.parent
REPO_ROOT = CASE_STUDY.parent
GRADE_PY = CASE_STUDY / "grade.py"
STRIP_SH = CASE_STUDY / "scripts" / "strip_skill.sh"
COPILOT_ASSETS = ("AGENTS.md", "skills", "context")

# Where the mock harness drops the component it "produces".
MOCK_COMPONENT_REL = Path("src") / "methods_normalization" / "log_normalization"
# Fixtures the mock harness copies for a good vs. bad outcome.
_FIXTURE = {
    "good": CASE_STUDY / "examples" / "skill_guided",
    "bad": CASE_STUDY / "examples" / "plain_agent_typical",
}


def _resolve(base: Path, value: str) -> Path:
    """Resolve ``value`` as an absolute path, anchoring relatives at ``base``."""
    p = Path(value).expanduser()
    return p if p.is_absolute() else base / p


def _skill_block(copilot_repo: Path, names: list[str] | None) -> str:
    """Concatenate skill bodies to inject into a skill-arm prompt.

    Agents differ in whether they auto-load skills (Gemini/Antigravity CLI won't
    open SKILL.md files on their own). Injecting the relevant playbooks into the
    prompt guarantees the skill *content* reaches the model, so the experiment
    measures the value of the content, not each agent's discovery mechanism.
    """
    skills_dir = copilot_repo / "skills"
    blocks = []
    for skill_md in sorted(skills_dir.glob("*/SKILL.md")):
        name = skill_md.parent.name
        if names and name not in names:
            continue
        text = skill_md.read_text(encoding="utf-8", errors="ignore")
        if text.startswith("---"):  # strip YAML frontmatter
            parts = text.split("---", 2)
            text = parts[2] if len(parts) == 3 else text
        blocks.append(f"\n\n## Reference playbook: {name}\n{text.strip()}")
    if not blocks:
        return ""
    return (
        "\n\n---\nThe following project reference playbooks describe the required "
        "conventions. Follow them.\n" + "".join(blocks)
    )


def _load_config(path: Path) -> dict:
    cfg = json.loads(path.read_text(encoding="utf-8"))
    cfg.setdefault("component_name", "log_normalization")
    cfg.setdefault("timeout_seconds", 1800)
    cfg.setdefault("run_viash", False)
    # Resolve paths against the repo root (derived from this script's location),
    # NOT the current working directory, so the runner works from any directory.
    copilot_repo = _resolve(REPO_ROOT, cfg.get("copilot_repo", ".")).resolve()
    cfg["copilot_repo"] = str(copilot_repo)
    cfg["workdir"] = str(_resolve(REPO_ROOT, cfg.get("workdir", "case-study/runs")))
    cfg["prompt_file"] = str(
        _resolve(copilot_repo, cfg.get("prompt_file", "case-study/task/prompt.txt"))
    )
    return cfg


def _ignore_heavy(_dir, names):
    return {n for n in names if n in {".git", "node_modules", "__pycache__", ".venv"}}


def prepare_workdir(cfg: dict, arm: dict, copilot_repo: Path, label: str) -> Path:
    """Copy the task repo and install/strip the skill for this arm."""
    workdir = Path(cfg["workdir"]).resolve() / label
    if workdir.exists():
        shutil.rmtree(workdir)
    task_repo = Path(cfg["task_repo"]).expanduser().resolve()
    shutil.copytree(task_repo, workdir, ignore=_ignore_heavy)

    if arm.get("skill"):
        # Make the canonical assets present (opencode/gemini read AGENTS.md /
        # GEMINI.md and the referenced skills/ + context/).
        for asset in COPILOT_ASSETS:
            srcp = copilot_repo / asset
            if srcp.is_dir():
                shutil.copytree(srcp, workdir / asset, dirs_exist_ok=True)
            elif srcp.is_file():
                shutil.copyfile(srcp, workdir / asset)
        # Project the harness-native files (GEMINI.md/.gemini, CLAUDE.md/.claude,
        # AGENTS.md/.codex, .mcp.json, ...).
        target = arm.get("agent_target", "codex")
        subprocess.run(
            [
                sys.executable,
                "-m",
                "installer",
                "--source",
                str(copilot_repo),
                "install",
                "--target",
                target,
                "--dest",
                str(workdir),
            ],
            check=True,
        )
    else:
        subprocess.run(["bash", str(STRIP_SH), str(workdir)], check=True)
    return workdir


def run_harness(cfg: dict, arm: dict, workdir: Path, prompt: str) -> dict:
    """Invoke the configured harness command in the workdir."""
    harness = cfg["harnesses"][arm["harness"]]
    command = [
        part.replace("{prompt}", prompt).replace("{model}", arm.get("model", ""))
        for part in harness["command"]
    ]
    log_path = workdir / "_agent_run.log"
    start = time.time()
    if not shutil.which(command[0]):
        return {
            "ran": False,
            "error": f"harness '{command[0]}' not on PATH",
            "log": None,
            "seconds": 0,
        }
    with open(log_path, "w", encoding="utf-8") as log:
        proc = subprocess.run(
            command,
            cwd=workdir,
            stdout=log,
            stderr=subprocess.STDOUT,
            timeout=cfg["timeout_seconds"],
        )
    return {
        "ran": True,
        "returncode": proc.returncode,
        "log": str(log_path),
        "seconds": round(time.time() - start, 1),
    }


def mock_harness(arm: dict, workdir: Path, mode: str) -> dict:
    """A fake agent for smoke-testing the pipeline without a real model.

    Writes a known component into the workdir so prepare -> find -> grade can be
    exercised offline. Modes:
      * ``good`` / ``bad``  -> always produce that fixture
      * ``skill-aware``     -> produce the good component when skill assets are
        present in the workdir, else the bad one (so the full pipeline reproduces
        the expected pass/fail pattern without burning tokens)
    """
    if mode == "skill-aware":
        skill_present = (workdir / "AGENTS.md").exists() or (
            workdir / "GEMINI.md"
        ).exists()
        quality = "good" if skill_present else "bad"
    else:
        quality = mode  # "good" or "bad"

    dest = workdir / MOCK_COMPONENT_REL
    dest.mkdir(parents=True, exist_ok=True)
    for f in _FIXTURE[quality].iterdir():
        if f.is_file():
            shutil.copyfile(f, dest / f.name)
    return {"ran": True, "mock": mode, "produced": quality, "returncode": 0}


def find_component(workdir: Path, component_name: str) -> Path | None:
    """Find the produced component dir (a config.vsh.yaml named/located by name)."""
    best = None
    for cfg_path in workdir.rglob("config.vsh.yaml"):
        text = cfg_path.read_text(encoding="utf-8", errors="ignore")
        if component_name in text or component_name in str(cfg_path):
            # Prefer one whose `name:` matches exactly.
            if f"name: {component_name}" in text:
                return cfg_path.parent
            best = cfg_path.parent
    return best


def grade_component(
    component_dir: Path, copilot_repo: Path, run_viash: Path | None
) -> dict:
    cmd = [sys.executable, str(GRADE_PY), str(component_dir), "--json"]
    if run_viash:
        cmd += ["--run-viash", str(run_viash)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError:
        return {"error": "grade.py failed", "stderr": proc.stderr, "score": 0}


def run_arm(
    cfg: dict,
    arm: dict,
    copilot_repo: Path,
    dry_run: bool,
    mock: str | None = None,
    repeat_idx: int | None = None,
) -> dict:
    prompt = Path(cfg["prompt_file"]).read_text(encoding="utf-8").strip()
    # Skill arms get the playbooks injected into the prompt, so the skill content
    # reaches the model regardless of whether the agent auto-loads skills.
    if arm.get("skill") and cfg.get("inject_skills", True):
        prompt += _skill_block(copilot_repo, cfg.get("inject_skill_names"))
    label = arm["name"] if repeat_idx is None else f"{arm['name']}__rep{repeat_idx}"
    result = {
        "arm": arm["name"],
        "harness": arm["harness"],
        "model": arm.get("model"),
        "skill": bool(arm.get("skill")),
    }
    if dry_run:
        action = (
            f"install skill (target={arm.get('agent_target', 'codex')})"
            if arm.get("skill")
            else "strip skill"
        )
        cmd = cfg["harnesses"][arm["harness"]]["command"]
        result["plan"] = {
            "workdir": str(Path(cfg["workdir"]) / label),
            "prep": action,
            "inject_skills": bool(arm.get("skill") and cfg.get("inject_skills", True)),
            "command": (
                f"<mock:{mock}>"
                if mock
                else [
                    p.replace("{prompt}", "<PROMPT>").replace(
                        "{model}", arm.get("model", "")
                    )
                    for p in cmd
                ]
            ),
        }
        return result

    workdir = prepare_workdir(cfg, arm, copilot_repo, label)
    if mock:
        run_info = mock_harness(arm, workdir, mock)
    else:
        run_info = run_harness(cfg, arm, workdir, prompt)
    result["run"] = run_info
    component = find_component(workdir, cfg["component_name"])
    if component is None:
        result["graded"] = {"score": 0, "note": "no component produced"}
    else:
        result["component"] = str(component.relative_to(workdir))
        viash_root = workdir if cfg.get("run_viash") else None
        result["graded"] = grade_component(component, copilot_repo, viash_root)
    return result


def run_arm_repeated(
    cfg: dict, arm: dict, copilot_repo: Path, mock: str | None, repeats: int
) -> dict:
    """Run an arm ``repeats`` times and aggregate (median + range + pass rate)."""
    runs = [
        run_arm(
            cfg,
            arm,
            copilot_repo,
            dry_run=False,
            mock=mock,
            repeat_idx=(i if repeats > 1 else None),
        )
        for i in range(repeats)
    ]
    scores = [r["graded"].get("score", 0) for r in runs]
    passes = [bool(r["graded"].get("static_pass")) for r in runs]
    return {
        "arm": arm["name"],
        "harness": arm["harness"],
        "model": arm.get("model"),
        "skill": bool(arm.get("skill")),
        "repeats": repeats,
        "scores": scores,
        "score_median": statistics.median(scores),
        "score_min": min(scores),
        "score_max": max(scores),
        "pass_rate": round(sum(passes) / repeats, 2),
        "runs": runs,
    }


def render_table(results: list[dict]) -> str:
    rows = [
        "| Arm | Harness | Model | Skill | Median /11 | Range | Pass-rate | N |",
        "| --- | --- | --- | :---: | :---: | :---: | :---: | :---: |",
    ]
    for r in results:
        median = r.get("score_median", "—")
        lo, hi = r.get("score_min"), r.get("score_max")
        rng = (
            f"{lo}–{hi}"
            if lo is not None and lo != hi
            else (str(lo) if lo is not None else "—")
        )
        pass_rate = f"{int(round(r['pass_rate'] * 100))}%" if "pass_rate" in r else "—"
        rows.append(
            f"| {r['arm']} | {r['harness']} | {r.get('model','')} | "
            f"{'✅' if r['skill'] else '❌'} | {median}/11 | {rng} | "
            f"{pass_rate} | {r.get('repeats', 1)} |"
        )
    return "\n".join(rows)


def _synthesize_task_repo(path: Path) -> None:
    """Create a minimal stand-in task repo so mock runs need no external clone."""
    (path / "src" / "api").mkdir(parents=True, exist_ok=True)
    (path / "src" / "methods_normalization").mkdir(parents=True, exist_ok=True)
    (path / "_viash.yaml").write_text(
        "name: task_ist_preprocessing\n", encoding="utf-8"
    )
    (path / "src" / "api" / "comp_normalization.yaml").write_text(
        "# stub component API for the normalization stage\n", encoding="utf-8"
    )
    (path / "README.md").write_text("# (synthetic mock task repo)\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--only", action="append", help="Run only these arm name(s).")
    parser.add_argument("--dry-run", action="store_true", help="Plan without running.")
    parser.add_argument(
        "--mock-harness",
        choices=["skill-aware", "good", "bad"],
        metavar="MODE",
        help="Use a fake agent instead of a real CLI (smoke test). "
        "'skill-aware' produces a good component when skill assets are present, "
        "else a bad one.",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=None,
        help="Runs per arm (overrides config 'repeats'; default 1). Use >=3 to "
        "report a median and beat single-run noise.",
    )
    args = parser.parse_args(argv)

    cfg = _load_config(args.config)
    copilot_repo = Path(cfg["copilot_repo"])

    task_repo_val = str(cfg.get("task_repo", "")).strip()
    task_repo = Path(task_repo_val).expanduser() if task_repo_val else None

    # In mock mode, synthesize a stand-in task repo if the configured one is absent
    # so the full prepare->run->grade pipeline runs with no external dependencies.
    if args.mock_harness and not args.dry_run:
        if task_repo is None or not task_repo.exists():
            synth = Path(cfg["workdir"]) / "_mock_task_repo"
            if synth.exists():
                shutil.rmtree(synth)
            _synthesize_task_repo(synth)
            cfg["task_repo"] = str(synth)
            print(f"(mock) synthesized stand-in task repo at {synth}")
    # Real run: fail early and clearly if task_repo is the unedited placeholder.
    elif not args.dry_run:
        if (
            not task_repo_val
            or task_repo_val.startswith("/path/to/")
            or task_repo is None
            or not task_repo.exists()
        ):
            print(
                "❌ 'task_repo' is not set to a real checkout "
                f"(current value: {task_repo_val or '(missing)'}).\n"
                "   Edit your arms config and point it at a local clone of\n"
                "   openproblems-bio/task_ist_preprocessing, for example:\n"
                '       "task_repo": "/home/you/task_ist_preprocessing"\n'
                "   Or smoke-test the pipeline offline first with no clone:\n"
                "       --mock-harness skill-aware",
                file=sys.stderr,
            )
            return 2

    arms = cfg["arms"]
    if args.only:
        arms = [a for a in arms if a["name"] in set(args.only)]
        if not arms:
            print(f"❌ no arms match {args.only}", file=sys.stderr)
            return 2

    if args.dry_run:
        results = [
            run_arm(cfg, arm, copilot_repo, True, args.mock_harness) for arm in arms
        ]
        print(json.dumps(results, indent=2))
        return 0

    repeats = max(1, args.repeats or int(cfg.get("repeats", 1)))
    results = [
        run_arm_repeated(cfg, arm, copilot_repo, args.mock_harness, repeats)
        for arm in arms
    ]

    out_dir = Path(cfg["workdir"]).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "results.json").write_text(
        json.dumps(results, indent=2), encoding="utf-8"
    )
    table = render_table(results)
    (out_dir / "results.md").write_text(table + "\n", encoding="utf-8")
    print(table)
    print(f"\nWrote {out_dir / 'results.json'} and results.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
