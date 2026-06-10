# Archived: Kiro hackathon spec (production-mcp-server)

These three files are the **original spec from the hackathon** that started this
project, authored in the [Kiro](https://kiro.dev) spec-driven IDE (they lived
under `.kiro/specs/production-mcp-server/`).

They are **historical** and partly **superseded**:

- They describe the original **Continue.dev + Docker + execution-engine** vision
  ("IDE (VSCode) + Continue.dev Extension + Local MCP Server", an Execution
  Engine, a State Manager). The project has since pivoted to a **provider-neutral**
  co-pilot (any MCP-capable agent) whose server is **read-only validation/analysis**.
- The current, authoritative direction lives in [`AGENTS.md`](../../../AGENTS.md)
  and the **Current Status** section of [`README.md`](../../../README.md).

Why they're kept: `requirements.md` and `design.md` still contain useful,
detailed acceptance criteria and architecture for the genuine **roadmap** items
(in-server Nextflow/Viash/Docker execution, workflow state management, and the
OpenProblems build/benchmark/submission tools). When those are eventually
implemented, this is the original requirements record to draw from — with the
Continue.dev-specific framing updated to the provider-agnostic model.

Do not treat these as current. Tool inventories, the architecture diagram, and
the task checkboxes here are frozen at hackathon state.
