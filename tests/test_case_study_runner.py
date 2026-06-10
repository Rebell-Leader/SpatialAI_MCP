"""Guard the case-study runner's wiring (no live agent needed)."""

import importlib.util
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
RUN_PY = REPO_ROOT / "case-study" / "runner" / "run.py"
EXAMPLES = REPO_ROOT / "case-study" / "examples"
EXAMPLE_CONFIG = REPO_ROOT / "case-study" / "runner" / "arms.example.json"

_spec = importlib.util.spec_from_file_location("case_study_run", RUN_PY)
run_mod = importlib.util.module_from_spec(_spec)
sys.modules["case_study_run"] = run_mod
_spec.loader.exec_module(run_mod)


def test_dry_run_returns_zero(capsys):
    rc = run_mod.main(["--config", str(EXAMPLE_CONFIG), "--dry-run"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "antigravity" in out
    assert "install skill" in out and "strip skill" in out


def test_find_component(tmp_path):
    comp = tmp_path / "src" / "methods_normalization" / "log_normalization"
    comp.mkdir(parents=True)
    (comp / "config.vsh.yaml").write_text("name: log_normalization\n", encoding="utf-8")
    found = run_mod.find_component(tmp_path, "log_normalization")
    assert found == comp


def test_find_component_missing(tmp_path):
    (tmp_path / "src").mkdir()
    assert run_mod.find_component(tmp_path, "log_normalization") is None


def test_grade_component_on_fixture():
    graded = run_mod.grade_component(
        EXAMPLES / "skill_guided", REPO_ROOT, run_viash=None
    )
    assert graded["score"] == 11
    assert graded["static_pass"] is True


def test_render_table_has_rows():
    results = [
        {
            "arm": "X",
            "harness": "gemini",
            "model": "m",
            "skill": True,
            "repeats": 3,
            "scores": [10, 11, 11],
            "score_median": 11,
            "score_min": 10,
            "score_max": 11,
            "pass_rate": 0.67,
        },
    ]
    table = run_mod.render_table(results)
    assert "| Arm |" in table
    assert "X" in table and "11/11" in table and "10–11" in table


def test_mock_harness_skill_aware(tmp_path):
    # Skill present -> good component.
    (tmp_path / "AGENTS.md").write_text("x", encoding="utf-8")
    info = run_mod.mock_harness({"name": "c"}, tmp_path, "skill-aware")
    assert info["produced"] == "good"
    comp = run_mod.find_component(tmp_path, "log_normalization")
    assert comp is not None
    assert run_mod.grade_component(comp, REPO_ROOT, None)["score"] == 11


def test_mock_harness_no_skill(tmp_path):
    # No skill assets -> bad component.
    info = run_mod.mock_harness({"name": "a"}, tmp_path, "skill-aware")
    assert info["produced"] == "bad"
    comp = run_mod.find_component(tmp_path, "log_normalization")
    assert run_mod.grade_component(comp, REPO_ROOT, None)["score"] <= 3


def test_mock_end_to_end(tmp_path):
    """Full prepare->mock->grade pipeline reproduces the expected pass/fail split."""
    config = {
        "copilot_repo": str(REPO_ROOT),
        "task_repo": str(tmp_path / "does_not_exist"),
        "workdir": str(tmp_path / "runs"),
        "arms": [
            {
                "name": "skill",
                "harness": "gemini",
                "model": "m",
                "skill": True,
                "agent_target": "gemini",
            },
            {"name": "noskill", "harness": "gemini", "model": "m", "skill": False},
        ],
        "harnesses": {"gemini": {"command": ["gemini", "-p", "{prompt}"]}},
    }
    cfg_path = tmp_path / "arms.json"
    cfg_path.write_text(json.dumps(config), encoding="utf-8")

    rc = run_mod.main(["--config", str(cfg_path), "--mock-harness", "skill-aware"])
    assert rc == 0

    results = json.loads((tmp_path / "runs" / "results.json").read_text())
    by_arm = {r["arm"]: r["score_median"] for r in results}
    assert by_arm["skill"] == 11
    assert by_arm["noskill"] == 0


def test_repeats_aggregate_with_median(tmp_path):
    """--repeats runs each arm N times and reports median/range/pass-rate."""
    config = {
        "copilot_repo": str(REPO_ROOT),
        "task_repo": str(tmp_path / "missing"),
        "workdir": str(tmp_path / "runs"),
        "arms": [{"name": "skill", "harness": "g", "skill": True,
                  "agent_target": "gemini"}],
        "harnesses": {"g": {"command": ["g", "-p", "{prompt}"]}},
    }
    cfg_path = tmp_path / "arms.json"
    cfg_path.write_text(json.dumps(config), encoding="utf-8")
    rc = run_mod.main(
        ["--config", str(cfg_path), "--mock-harness", "good", "--repeats", "3"]
    )
    assert rc == 0
    results = json.loads((tmp_path / "runs" / "results.json").read_text())
    r = results[0]
    assert r["repeats"] == 3
    assert r["scores"] == [11, 11, 11]
    assert r["score_median"] == 11
    assert r["pass_rate"] == 1.0


def test_skill_block_injects_playbooks():
    block = run_mod._skill_block(REPO_ROOT, None)
    assert "Reference playbook: viash-component-authoring" in block
    assert "openproblems/base_python" in block  # body content, not just a name


def test_config_paths_are_cwd_independent(tmp_path, monkeypatch):
    """Paths resolve against the repo root, not the current working directory."""
    monkeypatch.chdir(tmp_path)  # run from an unrelated directory
    cfg = run_mod._load_config(EXAMPLE_CONFIG)
    prompt = Path(cfg["prompt_file"])
    assert prompt.is_absolute() and prompt.is_file()
    assert prompt == REPO_ROOT / "case-study" / "task" / "prompt.txt"
    assert Path(cfg["copilot_repo"]) == REPO_ROOT
    assert Path(cfg["workdir"]).is_absolute()


def test_real_run_rejects_placeholder_task_repo(tmp_path, capsys):
    """A real run with the unedited placeholder fails clearly (exit 2)."""
    config = {
        "copilot_repo": str(REPO_ROOT),
        "task_repo": "/path/to/task_ist_preprocessing",
        "workdir": str(tmp_path / "runs"),
        "arms": [{"name": "x", "harness": "gemini", "model": "m", "skill": False}],
        "harnesses": {"gemini": {"command": ["gemini", "-p", "{prompt}"]}},
    }
    cfg_path = tmp_path / "arms.json"
    cfg_path.write_text(json.dumps(config), encoding="utf-8")
    rc = run_mod.main(["--config", str(cfg_path)])
    assert rc == 2
    assert "task_repo" in capsys.readouterr().err
