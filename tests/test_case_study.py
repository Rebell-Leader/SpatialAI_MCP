"""Guard the case-study grader: it must discriminate the reference fixtures."""

import importlib.util
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
GRADE_PY = REPO_ROOT / "case-study" / "grade.py"
EXAMPLES = REPO_ROOT / "case-study" / "examples"

# Load grade.py as a module (it lives outside the package, by design). It must be
# registered in sys.modules before exec so its dataclasses can resolve the
# string annotations created by `from __future__ import annotations`.
_spec = importlib.util.spec_from_file_location("case_study_grade", GRADE_PY)
grade_mod = importlib.util.module_from_spec(_spec)
sys.modules["case_study_grade"] = grade_mod
_spec.loader.exec_module(grade_mod)


def _score(example: str) -> int:
    config = grade_mod._find_config(EXAMPLES / example)
    return grade_mod.grade(config).score


def test_plain_agent_fixture_scores_low():
    report = grade_mod.grade(grade_mod._find_config(EXAMPLES / "plain_agent_typical"))
    assert report.score <= 3, f"plain fixture should score low, got {report.score}"
    assert not report.static_pass


def test_skill_guided_fixture_scores_high():
    report = grade_mod.grade(grade_mod._find_config(EXAMPLES / "skill_guided"))
    assert report.score >= 10, f"skill fixture should score high, got {report.score}"
    assert report.static_pass


def test_grader_discriminates():
    assert _score("skill_guided") - _score("plain_agent_typical") >= 8


def test_every_check_runs_on_both_fixtures():
    for example in ("plain_agent_typical", "skill_guided"):
        report = grade_mod.grade(grade_mod._find_config(EXAMPLES / example))
        assert len(report.checks) == 11
        assert report.max_score == 11
