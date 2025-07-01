#!/usr/bin/env python3
"""
Hugging Face Spaces Demo for OpenProblems Spatial Transcriptomics MCP Server

This is a demo version adapted for HF Spaces deployment that showcases
the MCP server capabilities in a user-friendly Gradio interface.
"""

import gradio as gr
import json
import os
from typing import Dict, Any, List


class MockMCPServer:
    """Mock MCP server for HF Spaces demo (without external tool dependencies)."""

    def __init__(self):
        self.tools_info = {
            "check_environment": "Check if bioinformatics tools are available",
            "validate_nextflow_config": "Validate Nextflow pipeline syntax",
            "run_nextflow_workflow": "Execute Nextflow workflows",
            "run_viash_component": "Run Viash components",
            "build_docker_image": "Build Docker containers",
            "analyze_nextflow_log": "Analyze pipeline execution logs",
            "read_file": "Read file contents",
            "write_file": "Write files",
            "list_directory": "List directory contents",
            "list_available_tools": "List all MCP tools",
            "echo_test": "Test MCP connectivity"
        }

        self.resources_info = {
            "server://status": "MCP server status and capabilities",
            "documentation://nextflow": "Nextflow best practices",
            "documentation://viash": "Viash component guidelines",
            "documentation://docker": "Docker optimization tips",
            "templates://spatial-workflows": "Spatial transcriptomics templates"
        }

    def check_environment(self, tools_to_check: str = "nextflow,viash,docker,java") -> str:
        """Mock environment check for HF Spaces."""
        tools = [tool.strip() for tool in tools_to_check.split(",")]

        # Simulate environment check results
        results = {
            "environment_check": {
                "timestamp": "2024-01-20T10:30:00Z",
                "platform": "Hugging Face Spaces (Ubuntu 20.04)",
                "python_version": "3.10.14"
            },
            "tools_status": {},
            "recommendations": []
        }

        # Mock results for demo
        for tool in tools:
            if tool == "docker":
                results["tools_status"][tool] = {
                    "available": False,
                    "version": None,
                    "status": "Not available in HF Spaces environment",
                    "required_for": "Container-based workflows"
                }
                results["recommendations"].append(f"For production: Install {tool} on your local system")
            else:
                results["tools_status"][tool] = {
                    "available": False,
                    "version": None,
                    "status": "Demo environment - tools not installed",
                    "install_command": f"Install with: curl -s https://get.{tool}.io | bash" if tool in ["nextflow", "viash"] else "sudo apt install openjdk-17-jre-headless"
                }

        results["summary"] = f"Demo mode: {len(tools)} tools checked, 0 available (expected in HF Spaces)"
        results["note"] = "This is a demo environment. In production, install tools locally for full functionality."

        return json.dumps(results, indent=2)

    def validate_nextflow_config(self, pipeline_content: str) -> str:
        """Mock Nextflow validation for demo."""
        if not pipeline_content.strip():
            return json.dumps({"error": "No pipeline content provided"}, indent=2)

        # Basic syntax checks for demo
        validation_results = {
            "validation_status": "demo_mode",
            "pipeline_analysis": {
                "dsl_version": "DSL2" if "nextflow.enable.dsl=2" in pipeline_content or "workflow {" in pipeline_content else "DSL1",
                "processes_found": pipeline_content.count("process "),
                "workflows_found": pipeline_content.count("workflow "),
                "includes_found": pipeline_content.count("include "),
                "line_count": len(pipeline_content.split('\n'))
            },
            "basic_checks": {
                "has_shebang": pipeline_content.startswith("#!/usr/bin/env nextflow"),
                "has_workflow_block": "workflow {" in pipeline_content,
                "has_process_definitions": "process " in pipeline_content,
                "uses_containers": "container " in pipeline_content or "docker" in pipeline_content,
            },
            "recommendations": [],
            "demo_note": "This is a syntax analysis demo. For full validation, use: nextflow config -check pipeline.nf"
        }

        # Add recommendations based on analysis
        if not validation_results["basic_checks"]["has_shebang"]:
            validation_results["recommendations"].append("Add shebang: #!/usr/bin/env nextflow")
        if not validation_results["basic_checks"]["uses_containers"]:
            validation_results["recommendations"].append("Consider using containers for reproducibility")
        if validation_results["pipeline_analysis"]["dsl_version"] == "DSL1":
            validation_results["recommendations"].append("Upgrade to DSL2 for better features")

        return json.dumps(validation_results, indent=2)

    def analyze_nextflow_log(self, log_content: str) -> str:
        """Mock log analysis for demo."""
        if not log_content.strip():
            return json.dumps({"error": "No log content provided"}, indent=2)

        analysis = {
            "log_analysis": {
                "total_lines": len(log_content.split('\n')),
                "timestamp": "Demo analysis",
                "log_size_chars": len(log_content)
            },
            "issues_found": [],
            "patterns_detected": [],
            "performance_indicators": {},
            "recommendations": []
        }

        # Pattern matching for common issues
        lines = log_content.split('\n')

        for line in lines:
            line_lower = line.lower()
            if "error" in line_lower:
                analysis["issues_found"].append({
                    "type": "error",
                    "line": line.strip(),
                    "pattern": "Error detected",
                    "suggestion": "Review error details and check input parameters"
                })
            elif "failed" in line_lower:
                analysis["issues_found"].append({
                    "type": "failure",
                    "line": line.strip(),
                    "pattern": "Process failure",
                    "suggestion": "Check process resource requirements and inputs"
                })
            elif "exit status 137" in line_lower:
                analysis["issues_found"].append({
                    "type": "oom",
                    "line": line.strip(),
                    "pattern": "Out of memory (exit status 137)",
                    "suggestion": "Increase memory allocation or optimize data processing"
                })

        # Detect patterns
        if "nextflow" in log_content.lower():
            analysis["patterns_detected"].append("Nextflow execution log")
        if "docker" in log_content.lower():
            analysis["patterns_detected"].append("Docker container usage")
        if "process >" in log_content:
            analysis["patterns_detected"].append("Process execution details")

        analysis["summary"] = f"Analyzed {len(lines)} lines, found {len(analysis['issues_found'])} potential issues"
        analysis["demo_note"] = "This is a pattern-based analysis demo. Full analysis requires log context."

        return json.dumps(analysis, indent=2)

    def get_documentation(self, doc_type: str) -> str:
        """Get sample documentation for demo."""
        docs = {
            "nextflow": """# Nextflow DSL2 Best Practices

## Overview
Nextflow enables scalable and reproducible scientific workflows using software containers.

## Essential DSL2 Patterns

### Basic Pipeline Structure
```nextflow
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {
    input_ch = Channel.fromPath(params.input)
    PROCESS_NAME(input_ch)
}

process PROCESS_NAME {
    container 'biocontainers/tool:version'

    input:
    path input_file

    output:
    path "output.txt"

    script:
    \"\"\"
    tool --input ${input_file} --output output.txt
    \"\"\"
}
```

## Resource Management
- Always specify memory and CPU requirements
- Use dynamic resource allocation for variable workloads
- Implement retry strategies for robust execution

## OpenProblems Integration
- Follow OpenProblems naming conventions
- Use standardized input/output formats (h5ad)
- Include comprehensive metadata and documentation
""",
            "viash": """# Viash Component Development Guide

## Component Structure
Every Viash component consists of:
- config.vsh.yaml: Component configuration
- script.py/R: Core functionality implementation
- test.py/R: Unit tests

## Best Practices
- Keep components focused on single tasks
- Use descriptive parameter names and types
- Include comprehensive help documentation
- Implement proper error handling
- Follow semantic versioning

## OpenProblems Standards
- Use h5ad format for single-cell data
- Include spatial coordinates in obsm['spatial']
- Validate input data structure
- Generate standardized output formats
""",
            "docker": """# Docker Optimization for Bioinformatics

## Multi-stage Builds
Use multi-stage builds to reduce image size:
```dockerfile
FROM python:3.10-slim as builder
RUN pip install --user package

FROM python:3.10-slim
COPY --from=builder /root/.local /root/.local
```

## Bioinformatics-Specific Tips
- Use biocontainers as base images when available
- Pin specific versions for reproducibility
- Optimize layer caching for iterative development
- Use .dockerignore to exclude large data files
""",
            "spatial-workflows": """# Spatial Transcriptomics Pipeline Templates

## 1. Basic Preprocessing Pipeline
```nextflow
process SPATIAL_QC {
    input: path spatial_data
    output: path "qc_results.h5ad"
    script:
    \"\"\"
    python qc_spatial.py --input ${spatial_data} --output qc_results.h5ad
    \"\"\"
}
```

## 2. Spatially Variable Genes
```nextflow
process FIND_SVG {
    input: path processed_data
    output: path "svg_results.csv"
    script:
    \"\"\"
    python spatial_variable_genes.py --input ${processed_data} --output svg_results.csv
    \"\"\"
}
```

## 3. Label Transfer
```nextflow
process LABEL_TRANSFER {
    input:
        path query_data
        path reference_data
    output: path "annotated_data.h5ad"
    script:
    \"\"\"
    python label_transfer.py --query ${query_data} --reference ${reference_data} --output annotated_data.h5ad
    \"\"\"
}
```
""",
            "server-status": json.dumps({
                "server_name": "OpenProblems Spatial Transcriptomics MCP",
                "version": "0.1.0",
                "status": "demo_mode",
                "environment": "Hugging Face Spaces",
                "capabilities": {
                    "nextflow_execution": "demo_mode",
                    "viash_components": "demo_mode",
                    "docker_builds": False,
                    "automated_testing": True,
                    "log_analysis": True,
                    "web_interface": True
                },
                "supported_formats": ["h5ad", "json", "yaml", "nf", "vsh.yaml"],
                "documentation_available": True,
                "demo_note": "This is a demonstration environment. Full functionality available in local deployment."
            }, indent=2)
        }

        return docs.get(doc_type, f"Documentation for {doc_type} not available in demo mode.")


def create_spatial_mcp_demo():
    """Create the HF Spaces demo interface."""

    mcp = MockMCPServer()

    # Custom CSS for better appearance
    css = """
    .gradio-container {
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    }
    .demo-header {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 20px;
        border-radius: 10px;
        margin-bottom: 20px;
    }
    .tool-section {
        border: 1px solid #e0e0e0;
        border-radius: 8px;
        padding: 20px;
        margin: 10px 0;
        background: #fafafa;
    }
    .success { color: #28a745; }
    .warning { color: #ffc107; }
    .error { color: #dc3545; }
    """

    with gr.Blocks(
        title="OpenProblems Spatial Transcriptomics MCP Server Demo",
        theme=gr.themes.Soft(),
        css=css
    ) as demo:

        gr.HTML("""
        <div class="demo-header">
            <h1>🧬 OpenProblems Spatial Transcriptomics MCP Server</h1>
            <h3>Interactive Demo - Model Context Protocol for AI-Powered Bioinformatics</h3>
            <p>🚀 This demo showcases the MCP server that enables AI agents like Continue.dev to automate spatial transcriptomics workflows</p>
        </div>
        """)

        gr.Markdown("""
        ## 🎯 What is this?

        This is a **Model Context Protocol (MCP) server** designed for spatial transcriptomics research. It provides:
        - **11 specialized tools** for workflow automation
        - **5 knowledge resources** with curated documentation
        - **AI agent integration** for Continue.dev and other MCP-compatible tools
        - **Production deployment** via Docker and local installation

        > **Note**: This is a demo environment. For full functionality with Nextflow, Viash, and Docker, deploy locally.
        """)

        with gr.Tabs():

            # Environment Check Tab
            with gr.Tab("🔧 Environment Validation"):
                gr.Markdown("### Check Bioinformatics Environment")
                gr.Markdown("*Verify that required tools are installed and configured properly.*")

                with gr.Row():
                    tools_input = gr.Textbox(
                        value="nextflow,viash,docker,java",
                        label="Tools to Check",
                        placeholder="Comma-separated list: nextflow,viash,docker,java",
                        info="Enter tools to validate in your environment"
                    )
                    check_btn = gr.Button("🔍 Check Environment", variant="primary")

                env_output = gr.JSON(
                    label="Environment Check Results",
                    show_label=True
                )

                check_btn.click(mcp.check_environment, tools_input, env_output)

                gr.Markdown("""
                **💡 What this tool does:**
                - Validates bioinformatics tool installations
                - Checks version compatibility
                - Provides installation recommendations
                - Assesses environment readiness for spatial workflows
                """)

            # Pipeline Validation Tab
            with gr.Tab("⚡ Pipeline Validation"):
                gr.Markdown("### Nextflow Pipeline Syntax Analysis")
                gr.Markdown("*Analyze Nextflow DSL2 pipelines for syntax and best practices.*")

                pipeline_input = gr.Textbox(
                    label="Nextflow Pipeline Code",
                    value="""#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow {
    input_ch = Channel.fromPath(params.input)
    SPATIAL_QC(input_ch)
}

process SPATIAL_QC {
    container 'biocontainers/scanpy:1.9.1'

    input:
    path spatial_data

    output:
    path "qc_results.h5ad"

    script:
    '''
    python -c "
    import scanpy as sc
    import squidpy as sq
    adata = sc.read_h5ad('${spatial_data}')
    # Quality control analysis
    sc.pp.calculate_qc_metrics(adata)
    adata.write('qc_results.h5ad')
    "
    '''
}""",
                    lines=20,
                    placeholder="Paste your Nextflow pipeline code here..."
                )

                validate_btn = gr.Button("🔍 Validate Pipeline", variant="primary")
                validation_output = gr.JSON(label="Validation Results")

                validate_btn.click(mcp.validate_nextflow_config, pipeline_input, validation_output)

                gr.Markdown("""
                **💡 What this tool does:**
                - Analyzes DSL2 syntax and structure
                - Checks for best practices compliance
                - Identifies potential issues and improvements
                - Validates container usage and resource specifications
                """)

            # Log Analysis Tab
            with gr.Tab("🔍 Log Analysis"):
                gr.Markdown("### Nextflow Execution Log Analysis")
                gr.Markdown("*AI-powered analysis of pipeline execution logs to identify issues and optimization opportunities.*")

                log_input = gr.Textbox(
                    label="Nextflow Log Content",
                    value="""N E X T F L O W  ~  version 23.04.0
Launching `main.nf` [abc123] DSL2 - revision: def456

executor >  local (4)
[12/abc123] process > SPATIAL_QC     [100%] 2 of 2 ✓
[34/def456] process > FIND_SVG       [ 50%] 1 of 2, failed: 1 ✗

ERROR ~ Error executing process > 'FIND_SVG'

Caused by:
  Process `FIND_SVG` terminated with an error exit status (137)

Command executed:
  python spatial_variable_genes.py --input data.h5ad --output svg_results.csv

Command exit status:
  137

Work dir:
  /work/34/def456...

Tip: you can replicate the issue by changing to the process work dir and entering the command shown above""",
                    lines=15,
                    placeholder="Paste Nextflow execution logs here..."
                )

                analyze_btn = gr.Button("🔍 Analyze Log", variant="primary")
                log_output = gr.JSON(label="Log Analysis Results")

                analyze_btn.click(mcp.analyze_nextflow_log, log_input, log_output)

                gr.Markdown("""
                **💡 What this tool does:**
                - Identifies common execution errors and failures
                - Detects out-of-memory issues (exit status 137)
                - Provides specific troubleshooting recommendations
                - Analyzes performance patterns and bottlenecks
                """)

            # Documentation Tab
            with gr.Tab("📚 Knowledge Resources"):
                gr.Markdown("### Access Curated Documentation")
                gr.Markdown("*Browse comprehensive documentation and templates for spatial transcriptomics workflows.*")

                doc_type = gr.Dropdown(
                    choices=[
                        ("Nextflow Best Practices", "nextflow"),
                        ("Viash Component Development", "viash"),
                        ("Docker Optimization", "docker"),
                        ("Spatial Workflow Templates", "spatial-workflows"),
                        ("Server Status", "server-status")
                    ],
                    value="nextflow",
                    label="Documentation Type",
                    info="Select documentation category to explore"
                )

                doc_btn = gr.Button("📖 Get Documentation", variant="primary")
                doc_output = gr.Textbox(
                    label="Documentation Content",
                    lines=20,
                    max_lines=30
                )

                doc_btn.click(mcp.get_documentation, doc_type, doc_output)

                gr.Markdown("""
                **💡 Available Resources:**
                - **Nextflow**: DSL2 patterns, resource management, OpenProblems integration
                - **Viash**: Component structure, best practices, testing guidelines
                - **Docker**: Multi-stage builds, bioinformatics optimization
                - **Spatial Templates**: Ready-to-use pipeline examples
                - **Server Status**: Current capabilities and configuration
                """)

            # MCP Integration Tab
            with gr.Tab("🤖 AI Agent Integration"):
                gr.Markdown("### Connect with Continue.dev and Other AI Agents")

                gr.Markdown("""
                ## 🚀 Local Installation & Integration

                To use this MCP server with AI agents like Continue.dev:

                ### 1. Install the MCP Server
                ```bash
                git clone https://github.com/openproblems-bio/SpatialAI_MCP.git
                cd SpatialAI_MCP
                pip install -e .
                ```

                ### 2. Configure Continue.dev
                Add this to your `~/.continue/config.json`:
                ```json
                {
                  "experimental": {
                    "modelContextProtocolServers": [
                      {
                        "name": "openproblems-spatial",
                        "transport": {
                          "type": "stdio",
                          "command": "python",
                          "args": ["-m", "mcp_server.main"],
                          "cwd": "/path/to/your/SpatialAI_MCP"
                        }
                      }
                    ]
                  }
                }
                ```

                ### 3. Test the Integration
                Ask your AI agent: *"Check my spatial transcriptomics environment and help me create a quality control pipeline"*

                ## 🛠️ Available MCP Tools
                """)

                # Display tools information
                tools_info = []
                for tool, desc in mcp.tools_info.items():
                    tools_info.append(f"• **{tool}**: {desc}")

                gr.Markdown("### Tools (11 available):\n" + "\n".join(tools_info))

                # Display resources information
                resources_info = []
                for resource, desc in mcp.resources_info.items():
                    resources_info.append(f"• **{resource}**: {desc}")

                gr.Markdown("### Resources (5 available):\n" + "\n".join(resources_info))

                gr.Markdown("""
                ## 🎯 Example AI Agent Interactions

                **User**: *"Help me set up spatial transcriptomics quality control"*

                **AI Agent Response**:
                ```
                I'll help you create a comprehensive spatial QC pipeline. Let me first assess your environment.

                [Uses check_environment tool]
                ✅ Docker: Available (version 28.1.1)
                ❌ Nextflow: Not found
                ❌ Viash: Not found

                [Uses list_directory tool]
                Found spatial data in: data/spatial_samples/
                Existing configs: config/

                Based on OpenProblems best practices, I'll:
                1. Install missing dependencies
                2. Create a modular QC pipeline
                3. Generate Viash components
                4. Set up comprehensive testing

                [Creates optimized pipeline with proper error handling and documentation]
                ```

                ## 📖 Additional Resources
                - **[Setup Guide](https://github.com/openproblems-bio/SpatialAI_MCP/blob/main/docs/CONTINUE_DEV_SETUP.md)**: Complete integration instructions
                - **[Agent Rules](https://github.com/openproblems-bio/SpatialAI_MCP/blob/main/docs/AGENT_RULES.md)**: Best practices for AI agents
                - **[Docker Deployment](https://github.com/openproblems-bio/SpatialAI_MCP/blob/main/docker/)**: Production deployment options
                """)

        gr.Markdown("""
        ---
        ## 🎉 Try It Yourself!

        1. **Explore the tools** above to see MCP capabilities in action
        2. **Install locally** for full Nextflow/Viash/Docker integration
        3. **Connect with Continue.dev** for AI-powered spatial transcriptomics workflows

        **🔗 Links**:
        [GitHub Repository](https://github.com/openproblems-bio/SpatialAI_MCP) |
        [OpenProblems Project](https://openproblems.bio) |
        [Model Context Protocol](https://modelcontextprotocol.io)

        *Transforming spatial transcriptomics research through AI-powered workflow automation.* 🧬✨
        """)

    return demo


# For HF Spaces deployment
if __name__ == "__main__":
    demo = create_spatial_mcp_demo()
    demo.launch(
        server_name="0.0.0.0",
        server_port=7860,
        show_error=True,
        share=False  # HF Spaces handles sharing
    )
