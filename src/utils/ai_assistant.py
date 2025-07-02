import openai
from typing import List, Dict

# LLM integration
try:
    import openai
    LLM_AVAILABLE = True
except ImportError:
    LLM_AVAILABLE = False
    print("OpenAI not available - install openai package")

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
        return '''You are an expert AI assistant for spatial transcriptomics and computational biology, integrated with the OpenProblems MCP server.

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

Be specific, practical, and focus on real-world spatial transcriptomics research scenarios.'''

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
            user_prompt = f'''Please analyze this spatial transcriptomics code:

```python
{content}
```

{doc_context}

Provide:
1. **Code Review**: Issues, bugs, improvements
2. **Best Practices**: Spatial transcriptomics patterns
3. **Optimization**: Memory, performance, reproducibility
4. **Integration**: OpenProblems/Nextflow/Viash compatibility
5. **Debugging Steps**: If issues found'''

        else:  # logs
            user_prompt = f'''Please analyze these execution logs:

```
{content}
```

{doc_context}

Provide:
1. **Error Identification**: What went wrong?
2. **Root Cause**: Why did it happen?
3. **Fix Recommendations**: Specific solutions
4. **Prevention**: Avoid future issues
5. **Resource Optimization**: Memory/CPU adjustments'''

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

