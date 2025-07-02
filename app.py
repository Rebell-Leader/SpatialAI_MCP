#!/usr/bin/env python3
"""
Enhanced OpenProblems Spatial Transcriptomics MCP Server Demo
"""

import gradio as gr
import json
import asyncio
from typing import Dict, List, Any
from dotenv import load_dotenv

# Import our existing MCP server tools
from src.mcp_server.main import (
    handle_call_tool,
    handle_read_resource,
    handle_list_tools,
    handle_list_resources
)

# Load environment variables from .env file
load_dotenv()

async def setup_enhanced_services(openai_key: str, mixedbread_key: str, base_url: str, model: str):
    """Setup both LLM and documentation indexer via MCP tools"""
    # This function can be used to pass API keys to the backend if needed.
    # For this demo, we assume the backend is already configured.
    return "✅ Services are managed by the backend and should be ready.", "🚀 Services Ready"

async def enhanced_index_documentation():
    """Index documentation with progress via MCP tool"""
    try:
        result = await handle_call_tool("index_documentation", {})
        return result[0].text, "📚 Indexed Successfully"
    except Exception as e:
        return f"❌ Indexing failed: {str(e)}", ""

async def enhanced_search_docs(query: str):
    """Enhanced documentation search via MCP tool"""
    if not query.strip():
        return "Please provide a search query"

    try:
        results_json = await handle_call_tool("search_documentation", {"query": query})
        results = json.loads(results_json[0].text)

        if not results or results[0]["score"] == 0:
            return "No results found. Make sure documentation is indexed first."

        response = f"🔍 **Search Results for: '{query}'**\n\n"
        for i, result in enumerate(results, 1):
            response += f"**{i}. Relevance Score: {result['score']:.3f}**\n"
            response += f"📄 **Source:** [{result['url']}]({result['url']})\n"
            response += f"📝 **Content:** {result['text']}\n\n---\n\n"

        return response
    except Exception as e:
        return f"Search error: {e}"

async def enhanced_analyze_code(code: str, context: str):
    """Enhanced code analysis with documentation context via MCP tool"""
    if not code.strip():
        return "Please provide code to analyze"

    try:
        result = await handle_call_tool("analyze_code", {"code": code, "context": context})
        return result[0].text
    except Exception as e:
        return f"Analysis error: {e}"

async def enhanced_debug_logs(logs: str):
    """Enhanced log analysis with documentation context via MCP tool"""
    if not logs.strip():
        return "Please provide logs to analyze"

    try:
        result = await handle_call_tool("debug_logs", {"logs": logs})
        return result[0].text
    except Exception as e:
        return f"Analysis error: {e}"

# Create enhanced Gradio interface
def create_enhanced_interface():
    with gr.Blocks(
        title="OpenProblems Spatial Transcriptomics - Enhanced MCP Demo",
        theme=gr.themes.Soft(),
        css='''
        .highlight { background-color: #f0f8ff; padding: 10px; border-radius: 5px; }
        .success { color: #28a745; }
        .error { color: #dc3545; }
        '''
    ) as demo:

        gr.Markdown('''
        # 🧬 OpenProblems Spatial Transcriptomics MCP Server
        ## Enhanced Demo with AI-Powered Tools

        **🚀 Interactive demonstration of the MCP server's capabilities.**
        ''')

        with gr.Tab("🚀 Enhanced Setup"):
            gr.Markdown("### Configure AI Services")
            gr.Markdown("Setup both LLM and MixedBread embeddings for the best experience")

            with gr.Row():
                with gr.Column():
                    openai_key = gr.Textbox(
                        label="OpenAI API Key",
                        type="password",
                        placeholder="sk-..."
                    )
                    base_url = gr.Textbox(
                        label="Base URL (Optional)",
                        placeholder="https://api.openai.com/v1"
                    )
                    model = gr.Textbox(
                        label="Model",
                        value="gpt-3.5-turbo"
                    )

                with gr.Column():
                    mixedbread_key = gr.Textbox(
                        label="MixedBread API Key (Optional)",
                        type="password",
                        placeholder="Enter MixedBread API key for better embeddings"
                    )
                    gr.Markdown('''
                    **🧠 MixedBread Embeddings:**
                    - Superior semantic search capabilities
                    - Optimized for scientific documentation
                    - [Get API key](https://www.mixedbread.ai/)
                    - Without key: Uses local FastEmbed for demo
                    ''')

            setup_btn = gr.Button("🚀 Initialize Enhanced Services", variant="primary", size="lg")
            setup_status = gr.Textbox(label="Setup Status", lines=4, interactive=False)
            setup_indicator = gr.Textbox(label="Status", visible=False)

            setup_btn.click(
                setup_enhanced_services,
                inputs=[openai_key, mixedbread_key, base_url, model],
                outputs=[setup_status, setup_indicator]
            )

        with gr.Tab("📚 Smart Documentation"):
            gr.Markdown("### Live Documentation Indexing and Search")

            index_btn = gr.Button("📥 Index Documentation", variant="primary", size="lg")
            index_status = gr.Textbox(label="Indexing Progress", lines=10, interactive=False)

            gr.Markdown("### 🔍 Semantic Documentation Search")
            with gr.Row():
                search_input = gr.Textbox(
                    label="Search Query",
                    placeholder="e.g., 'nextflow memory allocation', 'viash component testing', 'spatial data formats'",
                    scale=3
                )
                search_btn = gr.Button("🔍 Search", variant="secondary")

            search_output = gr.Markdown(label="Search Results")

            index_btn.click(enhanced_index_documentation, outputs=[index_status, gr.Textbox(visible=False)])
            search_btn.click(enhanced_search_docs, inputs=[search_input], outputs=[search_output])

        with gr.Tab("🔍 AI Code Analysis"):
            gr.Markdown("### Advanced Spatial Transcriptomics Code Review")

            code_input = gr.Code(
                label="Spatial Transcriptomics Code",
                language="python",
                lines=20
            )

            context_input = gr.Textbox(
                label="Analysis Context",
                placeholder="Describe what this code should accomplish (e.g., 'cell segmentation', 'quality control', 'spatial clustering')",
                lines=2
            )

            analyze_btn = gr.Button("🔍 Analyze with AI + Documentation", variant="primary", size="lg")
            analysis_output = gr.Markdown(label="Detailed Analysis")

            analyze_btn.click(
                enhanced_analyze_code,
                inputs=[code_input, context_input],
                outputs=[analysis_output]
            )

        with gr.Tab("🐛 Smart Log Debugging"):
            gr.Markdown("### AI-Powered Execution Log Analysis")

            logs_input = gr.Textbox(
                label="Execution Logs",
                placeholder='''Paste your Nextflow, Viash, or Docker logs here...''',
                lines=20
            )

            debug_btn = gr.Button("🐛 Debug with AI + Documentation", variant="primary", size="lg")
            debug_output = gr.Markdown(label="Debugging Analysis")

            debug_btn.click(enhanced_debug_logs, inputs=[logs_input], outputs=[debug_output])

        with gr.Tab("🛠️ MCP Tools"):
            gr.Markdown("### Core MCP Server Tools")
            gr.Markdown("*Directly interact with the MCP server's built-in tools for validation and analysis.*")

            with gr.Accordion("Environment Validation", open=False):
                gr.Markdown("#### Check Bioinformatics Environment")
                with gr.Row():
                    tools_input = gr.Textbox(
                        value="nextflow,viash,docker,java",
                        label="Tools to Check",
                        placeholder="Comma-separated list: nextflow,viash,docker,java",
                    )
                    check_btn = gr.Button("🔍 Check Environment", variant="secondary")
                env_output = gr.JSON(label="Environment Check Results")

                async def check_environment_mcp(tools_to_check: str) -> str:
                    try:
                        result = await handle_call_tool("check_environment", {"tools_to_check": tools_to_check})
                        return json.loads(result[0].text)
                    except Exception as e:
                        return {"error": str(e)}

                check_btn.click(check_environment_mcp, tools_input, env_output)


            with gr.Accordion("Pipeline Validation", open=False):
                gr.Markdown("#### Nextflow Pipeline Syntax Analysis")
                pipeline_input = gr.Code(
                    label="Nextflow Pipeline Code",
                    value='''#!/usr/bin/env nextflow
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

                        python -c "
                        import scanpy as sc
                        adata = sc.read_h5ad('${spatial_data}')
                        sc.pp.calculate_qc_metrics(adata)
                        adata.write('qc_results.h5ad')
                        "
                    }
                    ''',
                    language="markdown",
                    lines=15
                )
                validate_btn = gr.Button("🔍 Validate Pipeline", variant="secondary")
                validation_output = gr.JSON(label="Validation Results")

                async def validate_nextflow_config_mcp(pipeline_content: str) -> str:
                    try:
                        result = await handle_call_tool("validate_nextflow_config", {"pipeline_content": pipeline_content})
                        return json.loads(result[0].text)
                    except Exception as e:
                        return {"error": str(e)}

                validate_btn.click(validate_nextflow_config_mcp, pipeline_input, validation_output)

            with gr.Accordion("Log Analysis", open=False):
                gr.Markdown("#### Nextflow Execution Log Analysis")
                log_file_path_input = gr.Textbox(
                    label="Path to Nextflow Log File",
                    placeholder="e.g., /path/to/.nextflow.log",
                    value="work/.nextflow.log"
                )
                analyze_log_btn = gr.Button("🔍 Analyze Log", variant="secondary")
                log_analysis_output = gr.JSON(label="Log Analysis Results")

                async def analyze_nextflow_log_mcp(log_file_path: str) -> str:
                    try:
                        result = await handle_call_tool("analyze_nextflow_log", {"log_file_path": log_file_path})
                        return json.loads(result[0].text)
                    except Exception as e:
                        return {"error": str(e)}

                analyze_log_btn.click(analyze_nextflow_log_mcp, log_file_path_input, log_analysis_output)

        with gr.Tab("📚 Knowledge Resources"):
            gr.Markdown("### Access Curated Documentation")
            gr.Markdown("*Browse comprehensive documentation and templates for spatial transcriptomics workflows.*")

            doc_type = gr.Dropdown(
                choices=[
                    ("Nextflow Best Practices", "documentation://nextflow"),
                    ("Viash Component Development", "documentation://viash"),
                    ("Docker Optimization", "documentation://docker"),
                    ("Spatial Workflow Templates", "templates://spatial-workflows"),
                    ("Server Status", "server://status")
                ],
                value="documentation://nextflow",
                label="Documentation Type",
                info="Select documentation category to explore"
            )

            doc_btn = gr.Button("📖 Get Documentation", variant="primary")
            doc_output = gr.Textbox(
                label="Documentation Content",
                lines=20,
                max_lines=30
            )

            async def get_documentation_mcp(doc_uri: str) -> str:
                try:
                    content = await handle_read_resource(doc_uri)
                    try:
                        return json.dumps(json.loads(content), indent=2)
                    except (json.JSONDecodeError, TypeError):
                        return content
                except Exception as e:
                    return f"Error retrieving documentation: {str(e)}"

            doc_btn.click(get_documentation_mcp, doc_type, doc_output)

        with gr.Tab("🤖 AI Agent Integration"):
            gr.Markdown("### Connect with Continue.dev and Other AI Agents")

            gr.Markdown('''
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
                      "args": ["-m", "src.mcp_server.main"],
                      "cwd": "/path/to/your/SpatialAI_MCP"
                    }
                  }
                ]
              }
            }
            ```
            ''')

            gr.Markdown("## 🛠️ Available MCP Tools & Resources")
            with gr.Row():
                tools_info_output = gr.Markdown("Loading tools...")
                resources_info_output = gr.Markdown("Loading resources...")

            async def get_mcp_info():
                tools = await handle_list_tools()
                resources = await handle_list_resources()
                tools_info = [f"• **{tool.name}**: {tool.description}" for tool in tools]
                resources_info = [f"• **{resource.uri}**: {resource.description}" for resource in resources]
                return (f"### Tools ({len(tools)} available):\n" + "\n".join(tools_info),
                        f"### Resources ({len(resources)} available):\n" + "\n".join(resources_info))

            demo.load(get_mcp_info, outputs=[tools_info_output, resources_info_output])

    return demo

# Launch application
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Launch the Gradio interface for OpenProblems Spatial Transcriptomics MCP Server.")
    parser.add_argument("--server-port", type=int, default=7860, help="Port to run the Gradio server on.")
    parser.add_argument("--share", action="store_true", help="Share the Gradio app publicly.")
    args = parser.parse_args()

    demo = create_enhanced_interface()
    demo.launch(
        server_name="0.0.0.0",
        server_port=args.server_port,
        share=args.share
    )
