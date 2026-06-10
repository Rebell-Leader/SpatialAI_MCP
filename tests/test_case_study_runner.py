"""Guard the case-study runner's wiring (no live agent needed)."""

import importlib.util
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
    assert "A_gemini_noskill" in out
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
            "graded": {"score": 11, "max_score": 11, "static_pass": True},
            "component": "src/x",
        },
    ]
    table = run_mod.render_table(results)
    assert "| Arm |" in table
    assert "X" in table and "11/11" in table
