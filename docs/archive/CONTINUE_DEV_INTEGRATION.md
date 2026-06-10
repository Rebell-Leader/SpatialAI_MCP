# Continue.dev Integration Guide

This guide covers two approaches for integrating OpenProblems spatial transcriptomics documentation with Continue.dev:

1. **Enhanced MCP Server** (Primary approach - what we've built)
2. **Continue.dev Document Artifacts** (Alternative approach)

## 🎯 Approach 1: Enhanced MCP Server (RECOMMENDED)

Our OpenProblems MCP Server now provides **real, comprehensive documentation** from official sources through the Model Context Protocol.

### Features

✅ **Real-time documentation access** from official sources
✅ **Structured knowledge delivery** via MCP Resources
✅ **File system operations** for local development
✅ **Environment validation** and setup assistance
✅ **Pipeline creation and validation**
✅ **Automated documentation updates**

### Setup

#### 1. Install Dependencies
```bash
pip install -e .
```

#### 2. Download Real Documentation
```bash
openproblems-mcp download-docs
```

This command downloads and caches:
- **Nextflow Documentation** - Complete official docs from nextflow.io
- **Viash Documentation** - Comprehensive guides from viash.io
- **OpenProblems Documentation** - READMEs and guides from GitHub repositories
- **Docker Best Practices** - Bioinformatics-specific containerization patterns
- **Spatial Workflow Templates** - Ready-to-use pipeline templates

#### 3. Configure Continue.dev

Add to your Continue.dev configuration (`~/.continue/config.json`):

```json
{
  "mcpServers": {
    "openproblems": {
      "command": "python",
      "args": ["-m", "mcp_server.main"],
      "cwd": "/path/to/SpatialAI_MCP"
    }
  }
}
```

#### 4. Verify Integration
```bash
openproblems-mcp doctor --check-tools
openproblems-mcp info
```

### Continue.dev Workflow Example

Once configured, Continue.dev agents can:

```typescript
// Agent can access comprehensive documentation
const nextflowDocs = await mcp.readResource("documentation://nextflow");
const spatialTemplates = await mcp.readResource("templates://spatial-workflows");

// Agent can perform file operations
const projectFiles = await mcp.callTool("list_directory", { directory_path: "." });
const pipelineContent = await mcp.callTool("read_file", { file_path: "main.nf" });

// Agent can validate and create pipelines
const validation = await mcp.callTool("validate_nextflow_config", {
  pipeline_path: "main.nf"
});

// Agent can check environment setup
const environment = await mcp.callTool("check_environment", {});
```

### Available MCP Resources

| Resource URI | Content | Size |
|--------------|---------|------|
| `documentation://nextflow` | Complete Nextflow docs | ~50KB+ |
| `documentation://viash` | Complete Viash docs | ~30KB+ |
| `documentation://docker` | Bioinformatics Docker patterns | ~10KB |
| `templates://spatial-workflows` | Spatial pipeline templates | ~15KB |
| `server://status` | Server status and capabilities | ~1KB |

### Available MCP Tools

| Tool | Description | Use Case |
|------|-------------|----------|
| `read_file` | Read file contents | Analyze configs, scripts |
| `write_file` | Create/modify files | Generate pipelines, configs |
| `list_directory` | Navigate project structure | Explore repositories |
| `check_environment` | Validate tool installation | Setup verification |
| `validate_nextflow_config` | Pipeline syntax checking | Quality assurance |
| `run_nextflow_workflow` | Execute pipelines | Testing and deployment |
| `build_docker_image` | Container preparation | Environment setup |
| `analyze_nextflow_log` | Debug pipeline errors | Troubleshooting |

---

## 🔄 Approach 2: Continue.dev Document Artifacts (ALTERNATIVE)

For users who prefer to manage documentation directly in Continue.dev:

### Setup

#### 1. Download Documentation
```bash
openproblems-mcp download-docs
cd data/docs_cache
```

#### 2. Add to Continue.dev Documents

In Continue.dev, add these cached documentation files as document artifacts:

```
data/docs_cache/nextflow_docs.md
data/docs_cache/viash_docs.md
data/docs_cache/openproblems_docs.md
data/docs_cache/docker_docs.md
data/docs_cache/spatial_templates_docs.md
```

#### 3. Configure Continue.dev

Add to `~/.continue/config.json`:

```json
{
  "docs": [
    {
      "title": "Nextflow Documentation",
      "startUrl": "file:///path/to/SpatialAI_MCP/data/docs_cache/nextflow_docs.md"
    },
    {
      "title": "Viash Documentation",
      "startUrl": "file:///path/to/SpatialAI_MCP/data/docs_cache/viash_docs.md"
    },
    {
      "title": "OpenProblems Documentation",
      "startUrl": "file:///path/to/SpatialAI_MCP/data/docs_cache/openproblems_docs.md"
    },
    {
      "title": "Docker Best Practices",
      "startUrl": "file:///path/to/SpatialAI_MCP/data/docs_cache/docker_docs.md"
    },
    {
      "title": "Spatial Pipeline Templates",
      "startUrl": "file:///path/to/SpatialAI_MCP/data/docs_cache/spatial_templates_docs.md"
    }
  ]
}
```

### Pros and Cons

| | MCP Server Approach | Document Artifacts Approach |
|---|---|---|
| **Pros** | • Real-time access<br>• Structured delivery<br>• File operations<br>• Tool execution | • Simple setup<br>• Direct file access<br>• No server dependency |
| **Cons** | • Requires MCP setup<br>• More complex | • Manual updates<br>• No tool execution<br>• Static content |

---

## 🏆 Recommendation: Use Enhanced MCP Server

The **Enhanced MCP Server approach** is recommended because:

1. **Real-time Documentation** - Always up-to-date with official sources
2. **Interactive Capabilities** - Agent can perform actions, not just read docs
3. **Structured Knowledge** - Organized, searchable, contextual information
4. **Complete Workflow** - From documentation to execution
5. **Environment Integration** - Validates setup and provides guidance

### Example Continue.dev Agent Conversation

```
🧬 User: "Help me create a spatial transcriptomics quality control pipeline"

🤖 Agent: Let me help you with that! I'll:
1. Check your environment setup
2. Get the latest Nextflow best practices
3. Use spatial transcriptomics templates
4. Create an optimized pipeline for you

[Agent uses MCP tools to check environment, read documentation, and create pipeline]

✅ Agent: "I've created a spatial QC pipeline following OpenProblems standards.
The pipeline includes:
- Scanpy-based quality control
- Proper Docker containerization
- DSL2 Nextflow syntax
- Resource management
- Output publishing

Would you like me to validate the syntax and explain any part?"
```

---

## 🔧 Maintenance

### Updating Documentation
```bash
# Refresh all documentation
openproblems-mcp download-docs

# Check server status
openproblems-mcp doctor

# Test integration
openproblems-mcp tool check_environment
```

### Monitoring
```bash
# View cached documentation
ls -la data/docs_cache/

# Check server resources
openproblems-mcp info
```

---

## 🚀 Next Steps

1. **Set up the Enhanced MCP Server** using Approach 1
2. **Download real documentation** with `openproblems-mcp download-docs`
3. **Configure Continue.dev** to connect to the MCP server
4. **Test the integration** with spatial transcriptomics workflows
5. **Enjoy AI-assisted bioinformatics development!**

The integration provides computational biologists with **unprecedented AI assistance** for spatial transcriptomics pipeline development, combining the power of Continue.dev with comprehensive, real-time bioinformatics knowledge.
