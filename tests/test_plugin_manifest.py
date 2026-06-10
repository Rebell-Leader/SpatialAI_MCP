"""Validate the Claude Code plugin marketplace manifests.

These guard the publishing path: a malformed marketplace.json/plugin.json means
the repo can't be installed via `/plugin marketplace add`.
"""

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PLUGIN_DIR = REPO_ROOT / ".claude-plugin"


def _load(name: str) -> dict:
    return json.loads((PLUGIN_DIR / name).read_text(encoding="utf-8"))


def test_manifests_are_valid_json():
    assert _load("marketplace.json")
    assert _load("plugin.json")


def test_marketplace_required_fields():
    mp = _load("marketplace.json")
    # Required by the marketplace schema: name, owner, plugins.
    assert isinstance(mp["name"], str) and mp["name"]
    assert isinstance(mp["owner"], dict) and mp["owner"].get("name")
    assert isinstance(mp["plugins"], list) and mp["plugins"]
    for entry in mp["plugins"]:
        assert entry["name"], "each plugin entry needs a name"
        assert entry["source"], "each plugin entry needs a source"


def test_plugin_manifest_required_fields():
    plugin = _load("plugin.json")
    # `name` is the only strictly required field.
    assert plugin["name"]
    # keywords, if present, must be a list (a string is a load error).
    assert isinstance(plugin.get("keywords", []), list)


def test_marketplace_entry_matches_plugin():
    mp = _load("marketplace.json")
    plugin = _load("plugin.json")
    names = {e["name"] for e in mp["plugins"]}
    assert plugin["name"] in names, (
        "plugin.json name must match a marketplace entry name"
    )


def test_root_source_plugin_has_discoverable_components():
    """The './' plugin relies on auto-discovery: skills/ and .mcp.json must exist."""
    mp = _load("marketplace.json")
    root_entries = [e for e in mp["plugins"] if e["source"] in ("./", ".")]
    if not root_entries:
        return
    skills = list((REPO_ROOT / "skills").glob("*/SKILL.md"))
    assert skills, "root-source plugin declares no skills via skills/"
    assert (REPO_ROOT / ".mcp.json").is_file(), "root-source plugin needs .mcp.json"
