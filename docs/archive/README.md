# Archive — historical material (do not treat as current)

Everything in this directory is **frozen, historical** material from earlier
states of the project — mostly the hackathon era, when the co-pilot was framed
as a **Continue.dev + VSCode** integration with a **Dockerized, Gradio-based,
execution-first** server. The project has since pivoted to a **provider-neutral**
co-pilot (any MCP-capable agent) whose server is **read-only validation/analysis**.

For what the project actually is today, see:
[`README.md`](../../README.md), [`AGENTS.md`](../../AGENTS.md),
[`QUICKSTART.md`](../../QUICKSTART.md), [`context/`](../../context/),
and [`skills/`](../../skills/).

## Contents

| Path | What it was | Why archived |
| --- | --- | --- |
| `kiro-spec-production-mcp-server/` | Original Kiro IDE spec (requirements/design/tasks) | Continue.dev + Docker + execution-engine vision; superseded by AGENTS.md + README. Detailed roadmap requirements kept for reference. |
| `project_artifacts/` | Hackathon strategy doc, pitch script, onboarding notes | Describe the early vision; could mislead. |
| `CONTINUE_DEV_INTEGRATION.md`, `CONTINUE_DEV_SETUP.md` | Continue.dev-specific integration/setup guides | Project is now provider-agnostic; setup is in QUICKSTART.md + the installer. |
| `SETUP.md` | Old setup guide | References a nonexistent `mcp_server.main` module, an `echo_test` tool, and Docker deployment; superseded by QUICKSTART.md. |
| `continue_config_example.json` | Continue.dev `config.json` example | The installer now generates per-agent MCP config. |
| `examples/` | MCP client demos + Continue.dev demo | Reference a nonexistent `mcp_server.main` module and unimplemented tools. |
| `docker/` | Containerized Gradio demo (Dockerfile, compose) | Copies a nonexistent `app.py`/`requirements.txt`, runs a Gradio UI on port 7860; the Docker-first approach was explicitly discarded (see `ARCHITECTURE_PIVOT_SUMMARY.md`). |
| `AGENT_INTEGRATION_GUIDE.md` | Continue.dev agent integration guide | Cited a nonexistent `check_environment` tool. |
| `*_SUMMARY.md` | Development progress narratives | Point-in-time notes, not current docs. |

Nothing here is referenced by the live code, the installer, or the current docs.
