# Results (template — fill after running; do not fabricate)

> Record exact model IDs, the task-repo commit SHA, Viash version, and N (repeats).
> Report median (range) across repeats. Grade blind.

- **Task-repo commit:** `________`
- **Viash version:** `________`
- **Repeats per arm (N):** `____`
- **Max turns budget:** 12

## Headline table

| Arm | Model | Skill | Conformance /11 | Build | Test | Task success | Steps | Human validations |
| --- | --- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| A | _commercial_ | ❌ | | | | | | |
| B | _open-source_ | ❌ | | | | | | |
| C | _same as B_ | ✅ | | | | | | |
| D | _same as A_ | ✅ | | | | | | |

## Hypothesis check

- [ ] **C ≥ A** on conformance and task success
- [ ] **C** uses **fewer steps** than A
- [ ] **C** needs **fewer human validations** than A
- [ ] **C ≫ B** (skill effect on a fixed model is large)

## Per-arm intervention log

For each run, list the human validations and what they corrected, so the
"fewer validations" claim is auditable.

| Arm / seed | Turn | Operator message | Tag |
| --- | --- | --- | --- |
| | | | |

## Notes / caveats

_e.g. flaky build, model refused, ran out of turns, etc._
