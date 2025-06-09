#!/usr/bin/env python3
"""
Continue.dev + OpenProblems MCP Server Demo

This demonstrates how a Continue.dev agent would interact with our MCP server
to accomplish common computational biology tasks.

Scenario: AI agent helping computational biologist prepare and validate
spatial transcriptomics pipeline.
"""

import asyncio
import json
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

async def continue_dev_demo():
    """Simulate Continue.dev agent workflow with MCP server."""

    # Connect to MCP server (this would be automatic in Continue.dev)
    server_params = StdioServerParameters(
        command="python",
        args=["-m", "mcp_server.main"],
        env=None
    )

    async with stdio_client(server_params) as (read, write):
        async with ClientSession(read, write) as session:

            print("🤖 Continue.dev Agent: Starting spatial transcriptomics pipeline analysis...")

            # Step 1: Check environment setup
            print("\n📋 STEP 1: Checking computational environment...")
            env_result = await session.call_tool("check_environment", {})
            env_data = json.loads(env_result.content[0].text)

            print(f"   Environment Status: {env_data['overall_status']}")
            if env_data['tools']['docker']['available']:
                print("   ✅ Docker is available")
            else:
                print("   ❌ Docker not found")

            # Step 2: Explore project structure
            print("\n📁 STEP 2: Exploring project structure...")
            dir_result = await session.call_tool("list_directory", {"directory_path": "."})
            files = json.loads(dir_result.content[0].text)

            project_files = [f['name'] for f in files if not f['is_directory']]
            print(f"   Found {len(files)} items in project directory")
            print(f"   Key files: {', '.join(project_files[:5])}")

            # Step 3: Get best practices documentation
            print("\n📚 STEP 3: Retrieving Nextflow best practices...")
            nextflow_docs = await session.read_resource("documentation://nextflow")
            docs_preview = nextflow_docs.contents[0].text[:200] + "..."
            print(f"   Documentation loaded: {len(nextflow_docs.contents[0].text)} characters")
            print(f"   Preview: {docs_preview}")

            # Step 4: Create example pipeline file
            print("\n✏️ STEP 4: Creating example Nextflow pipeline...")
            example_pipeline = '''#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Spatial transcriptomics quality control pipeline
process SPATIAL_QC {
    container 'openproblems/spatial-transcriptomics:latest'

    input:
    path spatial_data

    output:
    path "qc_results.h5ad"
    path "qc_metrics.json"

    script:
    """
    python /app/spatial_qc.py \\
        --input ${spatial_data} \\
        --output qc_results.h5ad \\
        --metrics qc_metrics.json
    """
}

workflow {
    Channel.fromPath(params.input_dir + "/*.h5ad") | SPATIAL_QC
}
'''

            await session.call_tool("write_file", {
                "file_path": "example_spatial_pipeline.nf",
                "content": example_pipeline
            })
            print("   ✅ Created example_spatial_pipeline.nf")

            # Step 5: Validate the pipeline
            print("\n🔍 STEP 5: Validating pipeline syntax...")
            validation_result = await session.call_tool("validate_nextflow_config", {
                "pipeline_path": "example_spatial_pipeline.nf"
            })
            validation_data = json.loads(validation_result.content[0].text)

            print(f"   Validation status: {validation_data['status']}")
            if validation_data.get('warnings'):
                print(f"   Warnings: {len(validation_data['warnings'])}")
                for warning in validation_data['warnings']:
                    print(f"     ⚠️  {warning}")

            # Step 6: Get spatial workflow templates
            print("\n🧬 STEP 6: Loading spatial transcriptomics templates...")
            templates = await session.read_resource("templates://spatial-workflows")
            templates_content = templates.contents[0].text
            print(f"   Templates loaded: {len(templates_content)} characters")
            print("   Available workflow patterns for spatial analysis")

            print("\n🎉 Continue.dev Agent: Pipeline analysis complete!")
            print("   ✅ Environment checked")
            print("   ✅ Project structure mapped")
            print("   ✅ Best practices retrieved")
            print("   ✅ Example pipeline created")
            print("   ✅ Pipeline validated")
            print("   ✅ Templates ready for use")

            return {
                "environment": env_data,
                "validation": validation_data,
                "files_created": ["example_spatial_pipeline.nf"],
                "status": "ready_for_spatial_analysis"
            }

if __name__ == "__main__":
    result = asyncio.run(continue_dev_demo())
    print(f"\n📊 Final Result: {json.dumps(result, indent=2)}")
