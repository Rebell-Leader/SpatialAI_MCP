---
name: viash-component-authoring
description: Author or modify a Viash component (config.vsh.yaml + script) for an OpenProblems method, especially a task_ist_preprocessing stage. Use when the user wants to add/implement a spatial method, scaffold a component, or fix a config.vsh.yaml.
---

# Viash component authoring

A Viash component = a `config.vsh.yaml` (metadata: arguments, resources, engines)
+ a script (Python/R/Bash). Viash compiles it into a CLI, a Docker image, and a
Nextflow module.

## Steps

1. **Identify the stage and its interface.** OpenProblems components conform to a
   `comp_*` API merged in via `__merge__: /src/api/comp_<stage>.yaml`. Find the
   right interface for the pipeline stage you're implementing (see
   `context/ist-preprocessing-pipeline.md`) and match its inputs/outputs exactly.

2. **Write `config.vsh.yaml`.** Use the modern **flat** config (no
   `functionality:` wrapper) with `engines:` / `runners:` (not the legacy
   `platforms:`). Follow OpenProblems conventions: merge the stage API, build on
   the shared base image with a setup partial, include the standard metadata
   block, and label the nextflow runner with resource directives. This mirrors a
   real `task_ist_preprocessing` component:

   ```yaml
   __merge__: /src/api/comp_<stage>.yaml

   name: my_method
   label: "My Method"
   summary: "One-line biological purpose."
   description: "Longer description with method context."
   links:
     documentation: "https://github.com/openproblems-bio/task_ist_preprocessing"
     repository: "https://github.com/openproblems-bio/task_ist_preprocessing"
   references:
     doi: "10.xxxx/xxxxx"        # the method's paper

   arguments:
     - name: "--some_param"
       type: double
       default: 1.0
       description: "Biologically meaningful parameter."
   # Note: --input/--output are usually inherited from the merged comp_<stage> API.

   resources:
     - type: python_script
       path: script.py

   engines:
     - type: docker
       image: openproblems/base_python:1        # shared base, not raw python:*
       __merge__:
         - /src/base/setup_spatialdata_partial.yaml   # extra deps via setup partial
     - type: native

   runners:
     - type: executable
     - type: nextflow
       directives:
         label: [ midtime, lowcpu, lowmem ]     # resource labels
   ```

   Confirm the exact `__merge__` paths (the stage API and the base setup partial)
   and the current Viash schema against the live repo before building — keys
   evolve between major Viash versions.

3. **Write the script** with the dev parameter block:

   ```python
   import spatialdata as sd

   ## VIASH START
   par = {"input": "resources_test/.../dataset.zarr", "output": "output.zarr", "some_param": 1.0}
   ## VIASH END

   sdata = sd.read_zarr(par["input"])
   # ... implement; keep coordinate systems consistent ...
   result.write(par["output"])
   ```

   Fail loudly: log the error and `sys.exit(1)` on failure.

4. **Validate before building.** Run `analyze_workflow_configuration` on the
   `config.vsh.yaml` to catch structural/dependency issues.

5. **Build & test** with the user's local Viash (the server does not run it):
   - `viash ns build`
   - `viash run config.vsh.yaml -- --input <test>.zarr --output tmp/out.zarr`
   - unit test: add a `test_resources` script (e.g. `test.py`) to the config and
     run `viash test config.vsh.yaml` — every OpenProblems component ships a test.
   - integration: `nextflow run main.nf -profile test,docker`

## Don't

- Don't invent `create_spatial_component` / `build_viash_component` tool calls —
  they're roadmap, not implemented. Generate files and run `viash` in the terminal.
