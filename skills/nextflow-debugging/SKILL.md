---
name: nextflow-debugging
description: Diagnose a failing Nextflow pipeline run in an OpenProblems benchmark. Use when a `nextflow run` fails, a process errors out, or the user shares a .nextflow.log / .command.err and asks why it broke.
---

# Nextflow debugging

The server does not yet parse logs for you (`analyze_nextflow_log` is roadmap), so
triage the logs directly. Be systematic.

## Steps

1. **Find the failing task.** From the run output, note the process name and the
   work directory hash (`work/ab/cdef...`). The per-task files there are the
   ground truth:
   - `.command.sh` — the exact command run.
   - `.command.err` / `.command.out` — stderr/stdout.
   - `.command.log`, `.exitcode` — exit status.

2. **Classify the error by exit code / message:**
   - **137 / OOM-killed** → memory limit. Increase the process `memory` directive
     or add a dynamic retry (`memory { 8.GB * task.attempt }`).
   - **127 / "command not found"** → tool missing in the container; fix the
     component's engine/setup or the image.
   - **126 / permission denied** → executable bit / entrypoint issue.
   - **File/path errors** → check channel wiring and input staging; confirm the
     input file actually exists and matches the expected format
     (run `validate_spatial_data`).
   - **Config/DSL errors** → confirm DSL2 syntax; OpenProblems pipelines are DSL2.

3. **Add resilience** where appropriate: `errorStrategy 'retry'`, `maxRetries`,
   and dynamic resources by `task.attempt`. Don't paper over a real bug with
   retries.

4. **Reproduce in isolation.** Re-run a single component with `viash run` on the
   test data before re-running the whole pipeline. Use `-resume` to avoid
   recomputing successful tasks.

5. **Validate inputs** with `validate_spatial_data` / `analyze_spatial_metadata`
   when the failure looks data-shaped (wrong format, missing elements,
   raw-vs-normalized mismatch).

## Output

Report: the failing process, the root-cause class, the evidence (the log line),
and the specific fix — plus whether it's a data problem, an environment problem,
or a code problem.
