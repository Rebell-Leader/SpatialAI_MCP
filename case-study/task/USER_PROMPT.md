# Frozen user prompt (identical for every arm)

> Give this verbatim to the agent as the user's single opening message. Do not
> add hints. Any further operator messages count as "human validations".

---

I'm contributing to the OpenProblems `task_ist_preprocessing` benchmark. Please
add a new method component for the **normalization** stage that performs
**log-normalization**: scale each cell's counts to counts-per-10,000, then apply
`log1p`. Call it `log_normalization`.

It should be a proper component for this task so it builds in the namespace and
passes its test. Put it in the right place under `src/`, and let me know how to
build and test it.

---

## Notes for the operator (not shown to the agent)

- This prompt deliberately does **not** mention Viash syntax, SpatialData, base
  images, the component API, or tests. A domain-aware agent should know these
  from context; a general one must infer or ask.
- Stop the run at the first of: task success (see rubric), the max-turn budget
  (default **12** agent turns), or the operator giving up.
- Record each operator message after this one as a **human validation**, tagged
  with what it corrected (e.g. `legacy-syntax`, `wrong-data-format`, `missing-test`).
