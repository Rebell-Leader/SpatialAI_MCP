# 🚀 Enhanced MCP Demo Deployment Guide

## Overview

The enhanced demo transforms your basic Gradio app into a powerful AI-assisted spatial transcriptomics platform featuring:

- 🤖 **LLM Integration**: OpenAI-compatible API for code analysis and debugging
- 🧠 **MixedBread Embeddings**: State-of-the-art semantic search capabilities
- 📊 **Qdrant Vector Database**: Fast, scalable documentation search
- 🔧 **Live Documentation**: Real-time indexing from official sources
- 🧬 **Spatial Transcriptomics Focus**: Domain-specific AI assistance

## 🎯 Key Improvements Over Basic Demo

| Feature | Basic Demo | Enhanced Demo |
|---------|------------|---------------|
| **LLM Capability** | ❌ None | ✅ Full OpenAI integration with domain expertise |
| **Documentation** | ❌ Static mock data | ✅ Live scraping from official sources |
| **Search** | ❌ Basic text search | ✅ Semantic vector search with MixedBread |
| **Code Analysis** | ❌ Mock responses | ✅ AI-powered review with context |
| **Log Debugging** | ❌ Pattern matching only | ✅ Intelligent analysis with solutions |
| **Vector Database** | ❌ None | ✅ Qdrant with 1024-dim embeddings |
| **Real-time Updates** | ❌ Static | ✅ Live documentation indexing |

## 🔧 Deployment Options

### Option 1: Hugging Face Spaces (Recommended)

#### Step 1: Prepare Files
```bash
# Copy the enhanced demo files
cp app_mixedbread.py app.py
cp mixedbread_requirements.txt requirements.txt
```

#### Step 2: Create HF Spaces Repository
1. Go to [Hugging Face Spaces](https://huggingface.co/spaces)
2. Click "Create new Space"
3. Choose **Gradio** as the SDK
4. Upload these files:
   - `app.py` (the enhanced demo)
   - `requirements.txt` (dependencies)
   - `README.md` (optional description)

#### Step 3: Configure Environment
In your HF Space settings, add these environment variables:
```bash
# Optional: Pre-configure API keys (not recommended for public spaces)
OPENAI_API_KEY=your_openai_key_here
MIXEDBREAD_API_KEY=your_mixedbread_key_here
```

**Note**: For public demos, let users enter API keys through the UI for security.

### Option 2: Local Development

#### Step 1: Setup Environment
```bash
# Create virtual environment
python -m venv enhanced_demo_env
source enhanced_demo_env/bin/activate  # Linux/Mac
# enhanced_demo_env\Scripts\activate  # Windows

# Install dependencies
pip install -r mixedbread_requirements.txt
```

#### Step 2: Run Enhanced Demo
```bash
python app_mixedbread.py
```

#### Step 3: Access Interface
Open browser to: `http://localhost:7860`

### Option 3: Docker Deployment

#### Step 1: Create Dockerfile
```dockerfile
FROM python:3.9-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python packages
COPY mixedbread_requirements.txt requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Copy application
COPY app_mixedbread.py app.py

# Expose port
EXPOSE 7860

# Run application
CMD ["python", "app.py"]
```

#### Step 2: Build and Run
```bash
docker build -t openproblems-enhanced-demo .
docker run -p 7860:7860 openproblems-enhanced-demo
```

## 🔑 API Key Configuration

### Required Keys

#### 1. OpenAI API Key (Required for LLM features)
- **Purpose**: Code analysis, log debugging, AI assistance
- **Get Key**: [OpenAI Platform](https://platform.openai.com/api-keys)
- **Alternative**: Any OpenAI-compatible provider (Groq, Ollama, etc.)
- **Cost**: ~$0.01-0.03 per analysis (GPT-3.5-turbo)

#### 2. MixedBread API Key (Optional but Recommended)
- **Purpose**: Superior semantic embeddings for documentation search
- **Get Key**: [MixedBread AI](https://www.mixedbread.ai/)
- **Fallback**: Uses mock embeddings if no key provided
- **Cost**: Check MixedBread pricing for embedding API

### Alternative LLM Providers

You can use any OpenAI-compatible API:

```python
# Example configurations in the UI:

# Groq (fast inference)
Base URL: https://api.groq.com/openai/v1
Model: llama3-70b-8192

# Local Ollama
Base URL: http://localhost:11434/v1
Model: llama3

# Anthropic Claude (via compatible wrapper)
Base URL: https://api.anthropic.com/v1
Model: claude-3-sonnet-20240229
```

## 📚 Documentation Sources

The enhanced demo indexes live documentation from:

| Source | URL | Content Type |
|--------|-----|--------------|
| **Viash** | https://viash.io/docs/ | Component development |
| **Nextflow** | https://www.nextflow.io/docs/latest/ | Workflow orchestration |
| **OpenProblems** | https://openproblems.bio/documentation | Benchmarking framework |
| **Docker** | https://docs.docker.com/ | Containerization |

### Indexing Process
1. **Scraping**: Extracts clean text from HTML
2. **Chunking**: Splits into 800-character semantic chunks
3. **Embedding**: Generates 1024-dim vectors with MixedBread
4. **Storage**: Indexes in Qdrant vector database
5. **Search**: Enables semantic similarity search

## 🎮 Usage Workflows

### 🔍 Debugging Memory Issues

1. **Setup**: Add OpenAI API key in Enhanced Setup tab
2. **Index**: Click "Index Documentation with MixedBread"
3. **Search**: Query "nextflow memory allocation OOM error"
4. **Debug**: Paste actual error logs in Smart Log Debugging
5. **Analyze**: Get AI-powered root cause analysis + solutions

### 📝 Code Review Workflow

1. **Prepare**: Ensure both LLM and documentation are ready
2. **Submit**: Paste spatial transcriptomics code in AI Code Analysis
3. **Context**: Add description of what code should accomplish
4. **Review**: Get comprehensive analysis with best practices
5. **Improve**: Apply suggestions and re-analyze

### 🧬 Method Development

1. **Research**: Search for "viash component structure"
2. **Template**: Use search results to understand patterns
3. **Code**: Write spatial analysis component
4. **Debug**: Test and fix issues with AI assistance
5. **Validate**: Ensure OpenProblems compatibility

## 🔧 Advanced Configuration

### Custom Embedding Models

To use different embedding models, modify the `MixedBreadEmbeddings` class:

```python
class MixedBreadEmbeddings:
    def __init__(self, api_key: str = None, model: str = "mixedbread-ai/mxbai-embed-large-v1"):
        # Available models:
        # - mixedbread-ai/mxbai-embed-large-v1 (1024 dim, best quality)
        # - mixedbread-ai/mxbai-embed-base-v1 (768 dim, faster)
        # - Custom models via API
```

### Vector Database Customization

For production deployment with persistent storage:

```python
# Replace in-memory Qdrant with persistent storage
self.qdrant_client = QdrantClient(
    host="localhost",  # Your Qdrant server
    port=6333,
    # api_key="your_qdrant_api_key"  # For Qdrant Cloud
)
```

### Performance Tuning

```python
# Adjust chunk sizes for your use case
self.chunk_size = 800  # Smaller = more precise, larger = more context
self.chunk_overlap = 100  # Overlap between chunks

# Batch processing for large documents
batch_size = 100  # Process embeddings in batches
```

## 🚀 Production Deployment

### Performance Considerations

1. **Memory**: ~2-4GB RAM for full functionality
2. **Storage**: ~500MB for cached documentation
3. **Network**: Reliable internet for API calls and doc scraping
4. **Compute**: CPU sufficient, no GPU required

### Security Best Practices

1. **API Keys**: Never hardcode in public repositories
2. **Environment Variables**: Use secure key management
3. **Rate Limiting**: Implement usage limits for public demos
4. **Input Validation**: Sanitize user inputs

### Monitoring

```python
# Add logging for production monitoring
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('demo.log'),
        logging.StreamHandler()
    ]
)
```

## 🎯 Success Metrics

Track these metrics to measure demo effectiveness:

- **User Engagement**: Time spent, features used
- **API Usage**: Successful vs failed LLM/embedding calls
- **Search Quality**: Relevance scores, user feedback
- **Code Analysis**: Successful analyses vs errors
- **Documentation Hits**: Most searched topics

## 🔄 Continuous Improvement

### Regular Updates

1. **Documentation Refresh**: Re-index weekly for latest docs
2. **Model Updates**: Upgrade to newer LLM/embedding models
3. **Feature Additions**: Add new spatial transcriptomics tools
4. **Performance Optimization**: Monitor and improve response times

### User Feedback Integration

1. **Search Relevance**: Track which searches find helpful results
2. **Code Analysis Quality**: Monitor AI analysis accuracy
3. **Feature Requests**: Add most-requested capabilities
4. **Bug Reports**: Fix issues reported by users

## 📞 Support & Resources

### Documentation
- [Viash Documentation](https://viash.io/docs/)
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [OpenProblems Documentation](https://openproblems.bio/documentation)
- [MixedBread API Docs](https://www.mixedbread.ai/api-reference)
- [Qdrant Documentation](https://qdrant.tech/documentation/)

### Community
- [OpenProblems GitHub](https://github.com/openproblems-bio/openproblems)
- [Nextflow Community](https://www.nextflow.io/community.html)
- [Viash Community](https://viash.io/community/)

---

**🧬 Ready to deploy your enhanced spatial transcriptomics demo?**

Start with Option 1 (HF Spaces) for the easiest deployment, or Option 2 (Local) for development and testing! 🚀
