"""Command-line interface for the provider-agnostic skill installer.

Usage:
    spatialai-install install --target claude        # one provider
    spatialai-install install --target all           # every provider
    spatialai-install install --target cursor --dest /path/to/project
    spatialai-install list                           # show targets & skills
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import List

from .targets import Skill, Sources, TARGETS


def _find_repo_root(start: Path) -> Path:
    """Walk up from ``start`` to the directory containing AGENTS.md."""
    for candidate in [start, *start.parents]:
        if (candidate / "AGENTS.md").is_file():
            return candidate
    raise FileNotFoundError(
        "Could not locate AGENTS.md. Run inside the SpatialAI_MCP repo, or pass "
        "--source pointing at the repo root."
    )


_FRONTMATTER = re.compile(r"^---\s*\n(.*?)\n---\s*\n(.*)$", re.DOTALL)


def _parse_skill(skill_md: Path) -> Skill:
    """Parse a SKILL.md's YAML frontmatter without a yaml dependency."""
    text = skill_md.read_text(encoding="utf-8")
    match = _FRONTMATTER.match(text)
    name = skill_md.parent.name
    description = ""
    body = text
    if match:
        front, body = match.group(1), match.group(2)
        for line in front.splitlines():
            if ":" in line:
                key, _, value = line.partition(":")
                key, value = key.strip(), value.strip().strip("'\"")
                if key == "name" and value:
                    name = value
                elif key == "description":
                    description = value
    return Skill(name=name, description=description, path=skill_md, body=body)


def _load_sources(source_root: Path) -> Sources:
    skills_dir = source_root / "skills"
    skills: List[Skill] = []
    if skills_dir.is_dir():
        for skill_md in sorted(skills_dir.glob("*/SKILL.md")):
            skills.append(_parse_skill(skill_md))
    return Sources(
        repo_root=source_root,
        agents_md=(source_root / "AGENTS.md").read_text(encoding="utf-8"),
        skills=skills,
    )


def _cmd_list(source_root: Path) -> int:
    src = _load_sources(source_root)
    print("Supported targets:")
    for t in TARGETS.values():
        print(f"  {t.key:<8} {t.label}")
    print(f"\nSkills ({len(src.skills)}) discovered in {source_root / 'skills'}:")
    for s in src.skills:
        print(f"  {s.name:<28} {s.description[:70]}")
    return 0


def _cmd_install(source_root: Path, dest: Path, target_keys: List[str]) -> int:
    src = _load_sources(source_root)
    all_written: List[Path] = []
    for key in target_keys:
        target = TARGETS[key]
        written: List[Path] = []
        target.install(src, dest, written)
        all_written.extend(written)
        rels = [str(p.relative_to(dest)) for p in written]
        print(f"✅ {target.label}: wrote {len(written)} file(s)")
        for r in rels:
            print(f"     {r}")
    print(
        f"\nInstalled {len(target_keys)} target(s), "
        f"{len(all_written)} file(s) into {dest}"
    )
    return 0


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="spatialai-install",
        description="Install OpenProblems spatial co-pilot rules/skills into any "
        "AI coding agent (Claude Code, Codex, Cursor, Copilot).",
    )
    parser.add_argument(
        "--source",
        type=Path,
        default=None,
        help="Repo root holding AGENTS.md + skills/ (default: auto-detect).",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("list", help="List supported targets and discovered skills.")

    p_install = sub.add_parser(
        "install", help="Project rules/skills into one or more agents."
    )
    p_install.add_argument(
        "--target",
        "-t",
        action="append",
        choices=[*TARGETS.keys(), "all"],
        help="Target agent (repeatable). Use 'all' for every provider.",
    )
    p_install.add_argument(
        "--dest",
        "-d",
        type=Path,
        default=Path.cwd(),
        help="Destination project root (default: current directory).",
    )

    args = parser.parse_args(argv)

    try:
        source_root = (
            args.source.resolve()
            if args.source
            else _find_repo_root(Path(__file__).resolve().parent)
        )
    except FileNotFoundError as exc:
        print(f"❌ {exc}", file=sys.stderr)
        return 1

    if args.command == "list":
        return _cmd_list(source_root)

    targets = args.target or ["all"]
    if "all" in targets:
        targets = list(TARGETS.keys())
    # De-duplicate while preserving order.
    seen: set[str] = set()
    targets = [t for t in targets if not (t in seen or seen.add(t))]

    return _cmd_install(source_root, args.dest.resolve(), targets)


if __name__ == "__main__":
    raise SystemExit(main())
