#!/usr/bin/env python3
"""
Enhanced OpenProblems Spatial Transcriptomics MCP Server Demo
Featuring MixedBread embeddings and Qdrant vector search
Based on: https://qdrant.tech/documentation/embeddings/mixedbread/#using-mixedbread-with-qdrant
"""

import gradio as gr
import json
import asyncio
import subprocess
import tempfile
import os
from typing import Dict, List, Any, Optional
import requests
from pathlib import Path
import hashlib
import time
import re
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

# Vector search and embeddings
try:
    from qdrant_client import QdrantClient
    from qdrant_client.models import Distance, VectorParams, PointStruct, Batch
    import numpy as np
    VECTOR_SEARCH_AVAILABLE = True
except ImportError:
    VECTOR_SEARCH_AVAILABLE = False
    print("Vector search not available - install qdrant-client")

# LLM integration
try:
    import openai
    LLM_AVAILABLE = True
except ImportError:
    LLM_AVAILABLE = False
    print("OpenAI not available - install openai package")

# MixedBread embeddings
try:
    import requests as embedding_requests
    MIXEDBREAD_AVAILABLE = True
except ImportError:
    MIXEDBREAD_AVAILABLE = False

# FastEmbed for local embeddings
try:
    from fastembed import TextEmbedding
    FASTEMBED_AVAILABLE = True
except ImportError:
    FASTEMBED_AVAILABLE = False
    print("FastEmbed not available - install fastembed package for local embeddings")

class MixedBreadEmbeddings:
    """MixedBread embeddings API integration with FastEmbed fallback"""

    def __init__(self, api_key: str = None, model: str = "mixedbread-ai/mxbai-embed-large-v1"):
        self.api_key = api_key
        self.model = model
        self.base_url = "https://api.mixedbread.ai/v1"
        self.local_embedding_model = None
        if not self.api_key and FASTEMBED_AVAILABLE:
            print("Initializing local FastEmbed model...")
            self.local_embedding_model = TextEmbedding()
            print("Local FastEmbed model initialized.")

    def embed_texts(self, texts: List[str]) -> List[List[float]]:
        """Generate embeddings using MixedBread API or FastEmbed (local fallback)"""
        if self.api_key:
            try:
                headers = {
                    "Authorization": f"Bearer {self.api_key}",
                    "Content-Type": "application/json"
                }

                payload = {
                    "model": self.model,
                    "input": texts,
                    "encoding_format": "float"
                }

                response = embedding_requests.post(
                    f"{self.base_url}/embeddings",
                    headers=headers,
                    json=payload,
                    timeout=30
                )

                if response.status_code == 200:
                    data = response.json()
                    embeddings = [item["embedding"] for item in data["data"]]
                    return embeddings
                else:
                    print(f"MixedBread API error: {response.status_code}. Falling back to local embeddings.")
                    if self.local_embedding_model:
                        return self.local_embedding_model.embed(texts).tolist()
                    else:
                        return [[float(i % 1000) / 1000 for i in range(1024)] for _ in texts]

            except Exception as e:
                print(f"MixedBread embedding error: {e}. Falling back to local embeddings.")
                if self.local_embedding_model:
                    return self.local_embedding_model.embed(texts).tolist()
                else:
                    return [[float(i % 1000) / 1000 for i in range(1024)] for _ in texts]
        else:
            if self.local_embedding_model:
                print("Using local FastEmbed model for embeddings.")
                return self.local_embedding_model.embed(texts).tolist()
            else:
                print("No MixedBread API key and FastEmbed not available - using mock embeddings.")
                return [[float(i % 1000) / 1000 for i in range(1024)] for _ in texts]

class EnhancedDocumentationIndexer:
    """Enhanced documentation indexer with MixedBread embeddings"""

    def __init__(self, mixedbread_api_key: str = None):
        self.qdrant_client = None
        self.collection_name = "spatial_docs_enhanced"
        self.chunk_size = 800
        self.chunk_overlap = 100
        self.embeddings_model = MixedBreadEmbeddings(mixedbread_api_key)

        if VECTOR_SEARCH_AVAILABLE:
            self.setup_qdrant()

    def setup_qdrant(self):
        """Initialize Qdrant client (in-memory for demo)"""
        try:
            self.qdrant_client = QdrantClient(":memory:")
            # Create collection for documentation
            self.qdrant_client.create_collection(
                collection_name=self.collection_name,
                vectors_config=VectorParams(size=1024, distance=Distance.COSINE)
            )
            print("✅ Qdrant vector database initialized")
        except Exception as e:
            print(f"Failed to setup Qdrant: {e}")

    def extract_text_from_html(self, html_content: str) -> str:
        """Extract clean text from HTML content"""
        # Remove scripts and style elements
        html_content = re.sub(r'<script[^>]*>.*?</script>', '', html_content, flags=re.DOTALL | re.IGNORECASE)
        html_content = re.sub(r'<style[^>]*>.*?</style>', '', html_content, flags=re.DOTALL | re.IGNORECASE)

        # Remove HTML tags
        text = re.sub(r'<[^>]+>', ' ', html_content)

        # Clean up whitespace
        text = re.sub(r'\s+', ' ', text).strip()

        # Remove common navigation/footer text
        text = re.sub(r'(Copyright|©|\d{4}|Terms|Privacy|Cookie)', '', text, flags=re.IGNORECASE)

        return text[:100000]  # Limit for demo

    def chunk_text_smart(self, text: str, url: str) -> List[Dict]:
        """Intelligently chunk text based on structure"""
        chunks = []

        # Split by common delimiters (headers, paragraphs)
        sections = re.split(r'\n\s*\n|\. (?=[A-Z])', text)

        current_chunk = ""
        chunk_id = 0

        for section in sections:
            section = section.strip()
            if not section:
                continue

            # If adding this section would exceed chunk size, save current chunk
            if len(current_chunk) + len(section) > self.chunk_size and current_chunk:
                chunks.append({
                    "text": current_chunk.strip(),
                    "url": url,
                    "chunk_id": chunk_id,
                    "char_count": len(current_chunk)
                })
                current_chunk = section
                chunk_id += 1
            else:
                current_chunk += f" {section}" if current_chunk else section

        # Add the last chunk
        if current_chunk.strip():
            chunks.append({
                "text": current_chunk.strip(),
                "url": url,
                "chunk_id": chunk_id,
                "char_count": len(current_chunk)
            })

        return chunks

    def scrape_documentation_enhanced(self, urls: List[str]) -> Dict[str, Any]:
        """Enhanced documentation scraping with better text extraction"""
        docs = {}

        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
        }

        for url in urls:
            try:
                print(f"📄 Scraping: {url}")
                response = requests.get(url, headers=headers, timeout=15)

                if response.status_code == 200:
                    clean_text = self.extract_text_from_html(response.text)

                    # Extract title if possible
                    title_match = re.search(r'<title[^>]*>([^<]+)</title>', response.text, re.IGNORECASE)
                    title = title_match.group(1) if title_match else "Documentation"

                    docs[url] = {
                        "content": clean_text,
                        "title": title.strip(),
                        "timestamp": time.time(),
                        "char_count": len(clean_text),
                        "status": "success"
                    }
                    print(f"✅ {url}: {len(clean_text)} characters extracted")
                else:
                    docs[url] = {
                        "error": f"HTTP {response.status_code}",
                        "timestamp": time.time(),
                        "status": "error"
                    }
                    print(f"❌ {url}: HTTP {response.status_code}")

            except Exception as e:
                docs[url] = {
                    "error": str(e),
                    "timestamp": time.time(),
                    "status": "error"
                }
                print(f"❌ {url}: {str(e)}")

        return docs

    def index_documentation(self, progress_callback=None) -> str:
        """Index documentation with MixedBread embeddings"""
        if not self.qdrant_client:
            return "❌ Vector search not available - install qdrant-client"

        # Official documentation URLs
        urls = [
            "https://viash.io/docs/",
            "https://www.nextflow.io/docs/latest/",
            "https://openproblems.bio/documentation",
            "https://docs.docker.com/"
        ]

        results = []

        if progress_callback:
            progress_callback(0.1, "🌐 Scraping documentation...")

        docs = self.scrape_documentation_enhanced(urls)

        if progress_callback:
            progress_callback(0.3, "📝 Processing and chunking content...")

        all_chunks = []
        for url, doc_data in docs.items():
            if doc_data.get("status") == "success" and "content" in doc_data:
                chunks = self.chunk_text_smart(doc_data["content"], url)
                all_chunks.extend(chunks)
                results.append(f"✅ {doc_data.get('title', url)}: {len(chunks)} chunks ({doc_data['char_count']} chars)")
            else:
                results.append(f"❌ {url}: {doc_data.get('error', 'Failed')}")

        if progress_callback:
            progress_callback(0.6, "🧠 Generating MixedBread embeddings...")

        if all_chunks:
            # Generate embeddings with MixedBread
            texts = [f"{chunk['url']} | {chunk['text']}" for chunk in all_chunks]
            embeddings = self.embeddings_model.embed_texts(texts)

            if progress_callback:
                progress_callback(0.9, "💾 Storing in Qdrant vector database...")

            # Store in Qdrant using batch upsert
            points = []
            for i, (chunk, embedding) in enumerate(zip(all_chunks, embeddings)):
                points.append(PointStruct(
                    id=i,
                    vector=embedding,
                    payload={
                        **chunk,
                        "indexed_at": time.time()
                    }
                ))

            # Batch upsert for efficiency
            batch_size = 100
            for i in range(0, len(points), batch_size):
                batch = points[i:i + batch_size]
                self.qdrant_client.upsert(
                    collection_name=self.collection_name,
                    points=batch
                )

            results.append(f"\n🎉 Successfully indexed {len(all_chunks)} chunks with MixedBread embeddings!")
            results.append(f"📊 Vector database contains {len(embeddings)} vectors")

        if progress_callback:
            progress_callback(1.0, "✅ Documentation indexing complete!")

        return "\n".join(results)

    def search_documentation(self, query: str, limit: int = 5) -> List[Dict]:
        """Search indexed documentation using MixedBread embeddings"""
        if not self.qdrant_client:
            return [{"text": "Vector search not available", "url": "", "score": 0}]

        try:
            # Generate query embedding
            query_embedding = self.embeddings_model.embed_texts([query])[0]

            # Search in Qdrant
            search_results = self.qdrant_client.search(
                collection_name=self.collection_name,
                query_vector=query_embedding,
                limit=limit,
                with_payload=True
            )

            results = []
            for result in search_results:
                results.append({
                    "text": result.payload["text"][:600] + "..." if len(result.payload["text"]) > 600 else result.payload["text"],
                    "url": result.payload["url"],
                    "score": float(result.score),
                    "chunk_id": result.payload.get("chunk_id", 0),
                    "char_count": result.payload.get("char_count", 0)
                })

            return results
        except Exception as e:
            return [{"text": f"Search error: {str(e)}", "url": "", "score": 0}]

class LLMAssistant:
    """Enhanced LLM assistant for spatial transcriptomics"""

    def __init__(self):
        self.client = None
        self.model = "gpt-3.5-turbo"

    def setup_llm(self, api_key: str, base_url: str = "", model: str = "gpt-3.5-turbo"):
        """Setup LLM client"""
        if not LLM_AVAILABLE:
            return False, "OpenAI package not installed"

        try:
            if base_url:
                self.client = openai.OpenAI(api_key=api_key, base_url=base_url)
            else:
                self.client = openai.OpenAI(api_key=api_key)
            self.model = model

            # Test connection
            response = self.client.chat.completions.create(
                model=self.model,
                messages=[{"role": "user", "content": "Hello"}],
                max_tokens=5
            )

            return True, "LLM connected successfully"
        except Exception as e:
            return False, f"Connection failed: {str(e)}"

    def get_spatial_system_prompt(self) -> str:
        """Specialized system prompt for spatial transcriptomics"""
        return """You are an expert AI assistant for spatial transcriptomics and computational biology, integrated with the OpenProblems MCP server.

**Core Expertise:**
- Spatial transcriptomics analysis (SpatialData, AnnData, zarr)
- Nextflow pipeline development and optimization
- Viash component architecture and best practices
- Docker containerization for reproducible bioinformatics
- OpenProblems framework standards and benchmarking
- Python scientific stack (scanpy, spatialdata, anndata, numpy, pandas)

**Response Guidelines:**
1. **Biological Context First**: Always consider the biological meaning and data requirements
2. **Executable Solutions**: Provide complete, runnable code with proper error handling
3. **Best Practices**: Reference established patterns for reproducible research
4. **Troubleshooting Focus**: Help debug common issues (memory, containers, data formats)
5. **Framework Integration**: Suggest appropriate OpenProblems patterns and tools

**Code Standards:**
- Include proper imports and dependencies
- Add meaningful comments explaining biological steps
- Implement error handling for common failures
- Consider memory requirements for large spatial datasets
- Use established file formats (zarr, h5ad) correctly

Be specific, practical, and focus on real-world spatial transcriptomics research scenarios."""

    def analyze_with_context(self, content: str, analysis_type: str, docs_context: List[Dict] = None) -> str:
        """Analyze content with documentation context"""
        if not self.client:
            return "❌ LLM not connected. Please provide API key first."

        # Build documentation context
        doc_context = "\n\n📚 **Relevant Documentation Context:**\n"
        if docs_context:
            for i, doc in enumerate(docs_context, 1):
                doc_context += f"\n{i}. **Source:** {doc['url']}\n"
                doc_context += f"   **Relevance:** {doc['score']:.3f}\n"
                doc_context += f"   **Content:** {doc['text'][:400]}...\n"

        if analysis_type == "code":
            user_prompt = f"""Please analyze this spatial transcriptomics code:

```python
{content}
```

{doc_context}

Provide:
1. **Code Review**: Issues, bugs, improvements
2. **Best Practices**: Spatial transcriptomics patterns
3. **Optimization**: Memory, performance, reproducibility
4. **Integration**: OpenProblems/Nextflow/Viash compatibility
5. **Debugging Steps**: If issues found"""

        else:  # logs
            user_prompt = f"""Please analyze these execution logs:

```
{content}
```

{doc_context}

Provide:
1. **Error Identification**: What went wrong?
2. **Root Cause**: Why did it happen?
3. **Fix Recommendations**: Specific solutions
4. **Prevention**: Avoid future issues
5. **Resource Optimization**: Memory/CPU adjustments"""

        try:
            response = self.client.chat.completions.create(
                model=self.model,
                messages=[
                    {"role": "system", "content": self.get_spatial_system_prompt()},
                    {"role": "user", "content": user_prompt}
                ],
                max_tokens=1200,
                temperature=0.2
            )

            analysis = response.choices[0].message.content

            # Add documentation references
            if docs_context:
                analysis += "\n\n---\n\n📚 **Referenced Documentation:**\n"
                for doc in docs_context:
                    analysis += f"- [{doc['url']}]({doc['url']}) (relevance: {doc['score']:.3f})\n"

            return analysis

        except Exception as e:
            return f"❌ Analysis error: {str(e)}"

# Global instances
doc_indexer = None
llm_assistant = LLMAssistant()

def setup_enhanced_services(openai_key: str, mixedbread_key: str, base_url: str, model: str):
    """Setup both LLM and documentation indexer"""
    global doc_indexer

    results = []

    # Setup LLM
    if openai_key.strip():
        success, message = llm_assistant.setup_llm(openai_key, base_url, model)
        if success:
            results.append(f"✅ LLM: {message}")
        else:
            results.append(f"❌ LLM: {message}")
    else:
        results.append("⚠️ LLM: No API key provided")

    # Setup documentation indexer with MixedBread
    doc_indexer = EnhancedDocumentationIndexer(mixedbread_key.strip() if mixedbread_key.strip() else None)
    if mixedbread_key.strip():
        results.append("✅ MixedBread: API key configured")
    else:
        results.append("⚠️ MixedBread: Using mock embeddings (no API key)")

    return "\n".join(results), "🚀 Services Ready" if "✅ LLM" in results[0] else "⚠️ Partial Setup"

def enhanced_index_documentation():
    """Index documentation with progress"""
    if not doc_indexer:
        return "❌ Documentation indexer not initialized", ""

    try:
        result = doc_indexer.index_documentation()
        return result, "📚 Indexed Successfully"
    except Exception as e:
        return f"❌ Indexing failed: {str(e)}", ""

def enhanced_search_docs(query: str):
    """Enhanced documentation search"""
    if not doc_indexer:
        return "❌ Documentation indexer not initialized. Please setup services first."

    if not query.strip():
        return "Please provide a search query"

    results = doc_indexer.search_documentation(query, limit=5)

    if not results or results[0]["score"] == 0:
        return "No results found. Make sure documentation is indexed first."

    response = f"🔍 **Search Results for: '{query}'**\n\n"
    for i, result in enumerate(results, 1):
        response += f"**{i}. Relevance Score: {result['score']:.3f}**\n"
        response += f"📄 **Source:** [{result['url']}]({result['url']})\n"
        response += f"📝 **Content:** {result['text']}\n\n---\n\n"

    return response

def enhanced_analyze_code(code: str, context: str):
    """Enhanced code analysis with documentation context"""
    if not code.strip():
        return "Please provide code to analyze"

    # Search for relevant documentation
    search_query = f"spatial transcriptomics {context} {code[:150]}"
    docs_context = doc_indexer.search_documentation(search_query, limit=3) if doc_indexer else []

    # Get LLM analysis
    analysis = llm_assistant.analyze_with_context(code, "code", docs_context)

    return analysis

def enhanced_debug_logs(logs: str):
    """Enhanced log analysis with documentation context"""
    if not logs.strip():
        return "Please provide logs to analyze"

    # Search for relevant documentation based on log patterns
    search_query = f"nextflow error debugging {logs[:200]}"
    docs_context = doc_indexer.search_documentation(search_query, limit=3) if doc_indexer else []

    # Get LLM analysis
    analysis = llm_assistant.analyze_with_context(logs, "logs", docs_context)

    return analysis

# Create enhanced Gradio interface
def create_enhanced_interface():
    with gr.Blocks(
        title="OpenProblems Spatial Transcriptomics - Enhanced MCP Demo",
        theme=gr.themes.Soft(),
        css="""
        .highlight { background-color: #f0f8ff; padding: 10px; border-radius: 5px; }
        .success { color: #28a745; }
        .error { color: #dc3545; }
        """
    ) as demo:

        gr.Markdown("""
        # 🧬 OpenProblems Spatial Transcriptomics MCP Server
        ## Enhanced Demo with MixedBread Embeddings & Qdrant Vector Search

        **🚀 Interactive demonstration featuring:**
        - 🤖 **LLM-powered spatial transcriptomics assistance**
        - 🧠 **MixedBread embeddings for semantic search**
        - 📊 **Qdrant vector database for documentation**
        - 🔧 **Live documentation indexing and analysis**

        *Based on: [MixedBread + Qdrant Integration](https://qdrant.tech/documentation/embeddings/mixedbread/#using-mixedbread-with-qdrant)*
        """)

        # Enhanced Setup Tab
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
                    gr.Markdown("""
                    **🧠 MixedBread Embeddings:**
                    - Superior semantic search capabilities
                    - Optimized for scientific documentation
                    - [Get API key](https://www.mixedbread.ai/)
                    - Without key: Uses mock embeddings for demo
                    """)

            setup_btn = gr.Button("🚀 Initialize Enhanced Services", variant="primary", size="lg")
            setup_status = gr.Textbox(label="Setup Status", lines=4, interactive=False)
            setup_indicator = gr.Textbox(label="Status", visible=False)

            setup_btn.click(
                setup_enhanced_services,
                inputs=[openai_key, mixedbread_key, base_url, model],
                outputs=[setup_status, setup_indicator]
            )

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

            async def check_environment_mcp(tools_to_check: str) -> str:
                try:
                    result = await handle_call_tool("check_environment", {"tools_to_check": tools_to_check})
                    return json.dumps(result[0].text, indent=2)
                except Exception as e:
                    return json.dumps({"error": str(e)}, indent=2)

            check_btn.click(check_environment_mcp, tools_input, env_output)

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

            async def validate_nextflow_config_mcp(pipeline_content: str) -> str:
                try:
                    result = await handle_call_tool("validate_nextflow_config", {"pipeline_content": pipeline_content})
                    return json.dumps(result[0].text, indent=2)
                except Exception as e:
                    return json.dumps({"error": str(e)}, indent=2)

            validate_btn.click(validate_nextflow_config_mcp, pipeline_input, validation_output)

            gr.Markdown("""
            **💡 What this tool does:**
            - Analyzes DSL2 syntax and structure
            - Checks for best practices compliance
            - Identifies potential issues and improvements
            - Validates container usage and resource specifications
            """)

        # Log Analysis Tab (using MCP tool)
        with gr.Tab("🔍 Log Analysis"):
            gr.Markdown("### Nextflow Execution Log Analysis")
            gr.Markdown("*Analyze Nextflow execution logs using the MCP tool.*")

            log_file_path_input = gr.Textbox(
                label="Path to Nextflow Log File",
                placeholder="e.g., /path/to/.nextflow.log",
                value="work/.nextflow.log"
            )

            analyze_log_btn = gr.Button("🔍 Analyze Log (MCP Tool)", variant="primary")
            log_analysis_output = gr.JSON(label="Log Analysis Results (MCP Tool)")

            async def analyze_nextflow_log_mcp(log_file_path: str) -> str:
                try:
                    result = await handle_call_tool("analyze_nextflow_log", {"log_file_path": log_file_path})
                    return json.dumps(result[0].text, indent=2)
                except Exception as e:
                    return json.dumps({"error": str(e)}, indent=2)

            analyze_log_btn.click(analyze_nextflow_log_mcp, log_file_path_input, log_analysis_output)

            gr.Markdown("""
            **💡 What this tool does:**
            - Identifies common execution errors and failures
            - Detects out-of-memory issues (exit status 137)
            - Provides specific troubleshooting recommendations
            - Analyzes performance patterns and bottlenecks
            """)

        # Enhanced Documentation Tab
        with gr.Tab("📚 Smart Documentation"):
            gr.Markdown("### Live Documentation Indexing with MixedBread Embeddings")

            gr.Markdown("""
            **📄 Documentation Sources:**
            - 🧬 **[Viash](https://viash.io/docs/)**: Component development framework
            - 🔄 **[Nextflow](https://www.nextflow.io/docs/latest/)**: Workflow orchestration system
            - 🎯 **[OpenProblems](https://openproblems.bio/documentation)**: Benchmarking framework
            - 🐳 **[Docker](https://docs.docker.com/)**: Containerization platform

            **🧠 Enhanced with:**
            - MixedBread semantic embeddings
            - Qdrant vector database
            - Smart text chunking
            - Relevance scoring
            """)

            index_btn = gr.Button("📥 Index Documentation with MixedBread", variant="primary", size="lg")
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

            # Example searches
            gr.Markdown("### 💡 Example Searches")
            with gr.Row():
                gr.Button("Memory Issues").click(
                    lambda: "nextflow memory allocation OOM error debugging",
                    outputs=[search_input]
                )
                gr.Button("Spatial Data").click(
                    lambda: "spatial transcriptomics data formats SpatialData AnnData",
                    outputs=[search_input]
                )
                gr.Button("Component Testing").click(
                    lambda: "viash component testing docker platform",
                    outputs=[search_input]
                )

            index_btn.click(enhanced_index_documentation, outputs=[index_status])
            search_btn.click(enhanced_search_docs, inputs=[search_input], outputs=[search_output])

        # Enhanced Code Analysis
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

            # Code examples
            gr.Markdown("### 📋 Example Code Snippets")
            with gr.Row():
                gr.Button("Load Example: Basic QC").click(
                    lambda: ("""import spatialdata as sd
import scanpy as sc

# Load spatial data
sdata = sd.read_zarr('spatial_data.zarr')

# Basic quality control
adata = sdata.tables['table']
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")""", "Basic quality control for spatial transcriptomics data"),
                    outputs=[code_input, context_input]
                )

                gr.Button("Load Example: Memory Error").click(
                    lambda: ("""import spatialdata as sd
import numpy as np

# This might cause memory issues
sdata = sd.read_zarr('large_dataset.zarr')
huge_array = np.zeros((50000, 50000))  # Too large!
result = sdata.images['image'].values * huge_array

print(result.shape)""", "Code that might cause out of memory errors"),
                    outputs=[code_input, context_input]
                )

            analyze_btn.click(
                enhanced_analyze_code,
                inputs=[code_input, context_input],
                outputs=[analysis_output]
            )

        # Enhanced Log Analysis
        with gr.Tab("🐛 Smart Log Debugging"):
            gr.Markdown("### AI-Powered Execution Log Analysis")

            logs_input = gr.Textbox(
                label="Execution Logs",
                placeholder="""Paste your Nextflow, Viash, or Docker logs here...

Example Nextflow error:
ERROR ~ Error executing process > 'SPATIAL_ANALYSIS'

Caused by:
  Process `SPATIAL_ANALYSIS` terminated with an error exit status (137)

Command executed:
  python spatial_analysis.py --input data.zarr --output results.zarr

Command exit status:
  137

Work dir:
  /tmp/work/12/abc123def456

Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`""",
                lines=20
            )

            debug_btn = gr.Button("🐛 Debug with AI + Documentation", variant="primary", size="lg")
            debug_output = gr.Markdown(label="Debugging Analysis")

            # Log examples
            gr.Markdown("### 📋 Example Log Types")
            with gr.Row():
                gr.Button("Memory Error Log").click(
                    lambda: """ERROR ~ Error executing process > 'SPATIAL_QC'

Caused by:
  Process `SPATIAL_QC` terminated with an error exit status (137)

Command executed:
  python spatial_qc.py --input large_dataset.zarr --min_genes 200

Command exit status:
  137

Command error:
  Killed

Work dir:
  /tmp/nextflow-work/a1/b2c3d4e5f6""",
                    outputs=[logs_input]
                )

                gr.Button("Container Error Log").click(
                    lambda: """ERROR ~ Error executing process > 'BUILD_COMPONENT'

Caused by:
  Failed to pull Docker image 'openproblems/spatial:latest'

Command executed:
  docker run openproblems/spatial:latest viash run config.vsh.yaml

Command exit status:
  125

Command error:
  Unable to find image 'openproblems/spatial:latest' locally
  docker: Error response from daemon: pull access denied""",
                    outputs=[logs_input]
                )

            debug_btn.click(enhanced_debug_logs, inputs=[logs_input], outputs=[debug_output])

        # AI Agent Integration Tab
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

            tools_info_output = gr.Markdown("Loading tools...")
            resources_info_output = gr.Markdown("Loading resources...")

            async def get_mcp_info():
                tools = await handle_list_tools()
                resources = await handle_list_resources()

                tools_info = []
                for tool in tools:
                    tools_info.append(f"• **{tool.name}**: {tool.description}")

                resources_info = []
                for resource in resources:
                    resources_info.append(f"• **{resource.uri}**: {resource.description}")

                return ("### Tools ({} available):\n".format(len(tools)) + "\n".join(tools_info),
                        "### Resources ({} available):\n".format(len(resources)) + "\n".join(resources_info))

            demo.load(get_mcp_info, outputs=[tools_info_output, resources_info_output])

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

        # Knowledge Resources Tab (using MCP tool)
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
                    # Attempt to pretty print JSON if it's valid JSON
                    try:
                        return json.dumps(json.loads(content), indent=2)
                    except json.JSONDecodeError:
                        return content
                except Exception as e:
                    return f"Error retrieving documentation: {str(e)}"

            doc_btn.click(get_documentation_mcp, doc_type, doc_output)

            gr.Markdown("""
            **💡 Available Resources:**
            - **Nextflow**: DSL2 patterns, resource management, OpenProblems integration
            - **Viash**: Component structure, best practices, testing guidelines
            - **Docker**: Multi-stage builds, bioinformatics optimization
            - **Spatial Templates**: Ready-to-use pipeline examples
            - **Server Status**: Current capabilities and configuration
            """)

        # AI Agent Integration Tab
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

            tools_info_output = gr.Markdown("Loading tools...")
            resources_info_output = gr.Markdown("Loading resources...")

            async def get_mcp_info():
                tools = await handle_list_tools()
                resources = await handle_list_resources()

                tools_info = []
                for tool in tools:
                    tools_info.append(f"• **{tool.name}**: {tool.description}")

                resources_info = []
                for resource in resources:
                    resources_info.append(f"• **{resource.uri}**: {resource.description}")

                return ("### Tools ({} available):\n".format(len(tools)) + "\n".join(tools_info),
                        "### Resources ({} available):\n".format(len(resources)) + "\n".join(resources_info))

            demo.load(get_mcp_info, outputs=[tools_info_output, resources_info_output])

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

        # Usage Guide
        with gr.Tab("📖 Complete Guide"):
            gr.Markdown("""
            # 🚀 Complete Usage Guide

            ## 🔧 Quick Start

            ### 1. **Enhanced Setup** (🚀 Enhanced Setup tab)
            - Add your **OpenAI API key** for LLM assistance
            - Optional: Add **MixedBread API key** for superior embeddings
            - Click "Initialize Enhanced Services"

            ### 2. **Index Documentation** (📚 Smart Documentation tab)
            - Click "Index Documentation with MixedBread"
            - Wait for live scraping and vectorization
            - Try semantic searches with natural language

            ### 3. **Analyze Your Code** (🔍 AI Code Analysis tab)
            - Paste spatial transcriptomics code
            - Add context about what it should do
            - Get AI analysis with documentation references

            ### 4. **Debug Execution Issues** (🐛 Smart Log Debugging tab)
            - Paste Nextflow/Viash/Docker logs
            - Get intelligent debugging with fix recommendations

            ## 🧬 Spatial Transcriptomics Focus

            This demo specializes in:

            ### 📊 **Data Formats**
            - **SpatialData**: Modern spatial omics data structure
            - **AnnData**: Annotated data matrices for single cells
            - **Zarr**: Efficient storage for large datasets
            - **H5AD**: HDF5-based AnnData format

            ### 🔄 **Workflow Tools**
            - **Nextflow**: Pipeline orchestration and scalability
            - **Viash**: Component-based workflow development
            - **Docker**: Reproducible containerized environments
            - **OpenProblems**: Standardized benchmarking framework

            ### 🧪 **Analysis Types**
            - Cell segmentation and spatial mapping
            - Quality control and preprocessing
            - Spatial clustering and pattern detection
            - Integration with reference datasets

            ## 🎯 **Key Features**

            ### 🧠 **MixedBread Embeddings**
            - State-of-the-art semantic embeddings
            - Optimized for scientific documentation
            - Better understanding of technical concepts
            - Improved search relevance

            ### 📊 **Qdrant Vector Database**
            - Fast similarity search
            - Scalable vector storage
            - Real-time indexing
            - Precise relevance scoring

            ### 🤖 **AI-Powered Analysis**
            - Context-aware code review
            - Biological relevance checking
            - Performance optimization suggestions
            - Framework-specific best practices

            ## 💡 **Example Workflows**

            ### 🔍 **Debugging Workflow**
            1. Paste error logs → Get root cause analysis
            2. Search docs for "memory allocation" → Find solutions
            3. Apply fixes → Verify with code analysis

            ### 📝 **Code Review Workflow**
            1. Submit spatial analysis code → Get detailed review
            2. Search for "spatial clustering methods" → Learn best practices
            3. Implement suggestions → Re-analyze improved code

            ### 🧬 **Method Development Workflow**
            1. Search "viash component structure" → Understand framework
            2. Analyze example component code → Learn patterns
            3. Debug integration issues → Get specific fixes

            ## 🔗 **Integration Benefits**

            - **Live Documentation**: Always up-to-date information
            - **Semantic Search**: Natural language queries
            - **Context-Aware AI**: Domain-specific assistance
            - **Real-World Focus**: Production-ready solutions

            ## 🚀 **Next Steps**

            After using this demo:
            1. **Deploy the MCP server** in your local environment
            2. **Connect to Continue.dev** for IDE integration
            3. **Use real spatial datasets** with your workflows
            4. **Contribute to OpenProblems** with new methods

            ---

            **🧬 Ready to revolutionize your spatial transcriptomics workflow?**

            Start with the Enhanced Setup tab and explore the power of AI-assisted bioinformatics! 🚀
            """)

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

