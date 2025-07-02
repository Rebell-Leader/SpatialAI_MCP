---

title: OpenProblems Spatial Transcriptomics MCP Demo

emoji: 🧬

colorFrom: blue

colorTo: green

sdk: gradio

app_file: app.py

---

# 🚀 OpenProblems Spatial Transcriptomics MCP Server Demo

This Hugging Face Space hosts an interactive Gradio demo for the OpenProblems Spatial Transcriptomics MCP Server.

## ✨ Features

The demo showcases a powerful, AI-assisted platform for spatial transcriptomics research, including:

-   🤖 **LLM-Powered Assistance**: Integrates with OpenAI-compatible APIs for code analysis and debugging.
-   🧠 **Semantic Documentation Search**: Uses MixedBread embeddings (with a local FastEmbed fallback) and a Qdrant vector database to provide intelligent, real-time search across key bioinformatics documentation (Nextflow, Viash, Docker, OpenProblems).
-   🛠️ **Core MCP Tools**: Provides direct access to the MCP server for validating environments, checking pipeline syntax, and analyzing logs.
-   🧬 **Spatial Transcriptomics Focus**: Tailored prompts, examples, and tools specifically for spatial data analysis.

## 🔑 API Key Configuration

To use the full capabilities of this demo, you will need to configure API keys as **Secrets** in your Hugging Face Space settings:

1.  **`OPENAI_API_KEY`**: (Required) Your API key for an OpenAI-compatible service (e.g., OpenAI, Groq). This is necessary for all AI-powered analysis features.
2.  **`MIXEDBREAD_API_KEY`**: (Optional) Your API key for MixedBread AI. Providing this key enables higher-quality semantic embeddings for the documentation search. If omitted, the application will fall back to a local `FastEmbed` model, which runs on the Space's CPU.

## 🚀 How to Use

1.  **Setup**: Navigate to the **🚀 Enhanced Setup** tab and enter your API keys.
2.  **Index Docs**: Go to the **📚 Smart Documentation** tab and click the "Index Documentation" button to scrape and embed the latest documentation.
3.  **Explore**: Use the various tabs to analyze code, debug logs, and search for information using natural language.