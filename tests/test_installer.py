"""Tests for the provider-agnostic skill installer."""

import json
from pathlib import Path

import pytest

from installer.cli import _load_sources, _parse_skill, main
from installer.targets import TARGETS

REPO_ROOT = Path(__file__).resolve().parent.parent


def test_repo_has_canonical_assets():
    """The repo ships AGENTS.md and at least the four core skills."""
    src = _load_sources(REPO_ROOT)
    assert src.agents_md.strip()
    names = {s.name for s in src.skills}
    assert {
        "spatial-data-validation",
        "viash-component-authoring",
        "nextflow-debugging",
        "openproblems-submission",
    } <= names


def test_every_skill_has_a_description():
    src = _load_sources(REPO_ROOT)
    for skill in src.skills:
        assert skill.description, f"{skill.name} is missing a description"


def test_parse_skill_reads_frontmatter(tmp_path):
    skill_dir = tmp_path / "demo"
    skill_dir.mkdir()
    md = skill_dir / "SKILL.md"
    md.write_text(
        "---\nname: demo-skill\ndescription: A demo.\n---\n\n# Body\ntext\n",
        encoding="utf-8",
    )
    skill = _parse_skill(md)
    assert skill.name == "demo-skill"
    assert skill.description == "A demo."
    assert "# Body" in skill.body


@pytest.mark.parametrize("target_key", list(TARGETS.keys()))
def test_each_target_installs(tmp_path, target_key):
    src = _load_sources(REPO_ROOT)
    written: list[Path] = []
    TARGETS[target_key].install(src, tmp_path, written)
    assert written, f"{target_key} wrote no files"
    for path in written:
        assert path.exists()


def test_claude_outputs(tmp_path):
    src = _load_sources(REPO_ROOT)
    TARGETS["claude"].install(src, tmp_path, [])
    assert (tmp_path / "CLAUDE.md").exists()
    assert "@AGENTS.md" in (tmp_path / "CLAUDE.md").read_text()
    mcp = json.loads((tmp_path / ".mcp.json").read_text())
    assert "openproblems-spatial" in mcp["mcpServers"]
    # Each skill is copied into .claude/skills/<name>/SKILL.md
    for skill in src.skills:
        assert (tmp_path / ".claude" / "skills" / skill.name / "SKILL.md").exists()


def test_cursor_rules_have_frontmatter(tmp_path):
    src = _load_sources(REPO_ROOT)
    TARGETS["cursor"].install(src, tmp_path, [])
    main_rule = (tmp_path / ".cursor" / "rules" / "openproblems.mdc").read_text()
    assert main_rule.startswith("---")
    assert "alwaysApply: true" in main_rule
    for skill in src.skills:
        rule = (tmp_path / ".cursor" / "rules" / f"{skill.name}.mdc").read_text()
        assert "alwaysApply: false" in rule


def test_copilot_is_self_contained(tmp_path):
    src = _load_sources(REPO_ROOT)
    TARGETS["copilot"].install(src, tmp_path, [])
    content = (tmp_path / ".github" / "copilot-instructions.md").read_text()
    # Copilot has no import mechanism, so AGENTS.md content must be embedded.
    assert "OpenProblems Spatial Transcriptomics Co-pilot" in content


def test_cli_install_all(tmp_path):
    rc = main(["install", "--target", "all", "--dest", str(tmp_path)])
    assert rc == 0
    assert (tmp_path / "CLAUDE.md").exists()
    assert (tmp_path / ".cursor" / "mcp.json").exists()
    assert (tmp_path / ".github" / "copilot-instructions.md").exists()


def test_cli_list_runs():
    assert main(["list"]) == 0
