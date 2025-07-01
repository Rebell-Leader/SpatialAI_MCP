# 🎯 4-Minute Demo Presentation Script
**OpenProblems Spatial Transcriptomics MCP Server**

---

## **[0:00-0:45] Project Introduction & Problem Statement**

"Hello! Today I'm presenting the **OpenProblems Spatial Transcriptomics MCP Server** - a solution addressing a critical challenge in computational biology.

**The Problem**: Computational biologists working with spatial transcriptomics face enormous complexity - managing 10-100x larger datasets, dealing with multiple tool ecosystems like Nextflow, Viash, and Docker, while trying to focus on their actual scientific research rather than wrestling with technical infrastructure.

**Our Solution**: We've developed a Model Context Protocol server that bridges AI agents with bioinformatics tools, specifically designed for spatial transcriptomics workflows. This creates an intelligent assistant that understands both the biological context and the technical requirements of spatial data analysis."

---

## **[0:45-1:30] MCP Capabilities Overview**

"Our MCP server implements the full Model Context Protocol specification with **14 specialized tools** and **5 resource endpoints**:

**Core Infrastructure Tools**:
- Environment checking for Nextflow, Viash, Docker, and Java
- Nextflow workflow execution with proper parameter handling
- Viash component building and testing
- Docker image management and containerization

**Spatial Transcriptomics Specialized Tools**:
- Spatial component template generation
- Spatial data validation for SpatialData and AnnData formats
- Environment setup with spatial transcriptomics dependencies

**Resource Endpoints** provide:
- Live server status and capabilities
- Curated documentation for Nextflow, Viash, and Docker
- Spatial workflow templates and best practices

The key innovation is that this isn't just another chatbot - it's a **production-ready tool integration system** that can actually execute commands, validate data, and manage complex bioinformatics workflows."

---

## **[1:30-2:15] Demo UI Overview & Enhanced Setup**

"Let me walk you through our Gradio demonstration interface, which showcases these capabilities:

**Enhanced Setup Tab** - *[Video shows the setup interface]*
- This is where users configure their AI services
- **OpenAI API key input** for LLM-powered analysis
- **MixedBread API key** for superior semantic embeddings
- **Base URL and model selection** supporting alternative providers like Groq, Ollama, or local models
- **Connection testing** with real-time status indicators

The setup supports both required services for full functionality and optional enhancements for improved performance."

---

## **[2:15-3:00] Documentation & Analysis Capabilities**

"**Smart Documentation Tab** - *[Video shows documentation indexing]*
- **Live documentation scraping** from official sources: Viash, Nextflow, OpenProblems, and Docker
- **Real-time indexing** with MixedBread embeddings and Qdrant vector database
- **Semantic search capabilities** - users can query with natural language like 'nextflow memory allocation' or 'spatial data formats'
- **Relevance scoring** shows the most helpful documentation chunks

**AI Code Analysis Tab** - *[Video shows code input and analysis]*
- **Syntax-highlighted code editor** for spatial transcriptomics scripts
- **Context input** where users describe what their code should accomplish
- **Example code snippets** for testing - quality control, memory error scenarios
- **AI-powered analysis** that would provide detailed reviews with biological context and best practices"

---

## **[3:00-3:45] Debugging & Tools Integration**

"**Smart Log Debugging Tab** - *[Video shows log analysis interface]*
- **Log input area** for Nextflow, Viash, or Docker execution logs
- **Example log patterns** users can test with - memory errors, container issues
- **AI analysis output** that would identify root causes and provide specific fix recommendations

**MCP Tools Tab** - *[Video shows tool selection and execution]*
- **Tool dropdown** with all 14 available MCP tools
- **JSON argument input** for tool parameters
- **Quick example buttons** for common operations like environment checking, component creation, log analysis
- **Real-time tool execution** showing actual MCP server responses

The tools tab demonstrates the **bidirectional communication** between the AI agent and the bioinformatics infrastructure."

---

## **[3:45-4:00] Impact & Future Vision**

"This system represents a paradigm shift in computational biology workflows. Instead of researchers spending hours debugging Nextflow pipelines or configuring Docker containers, they can focus on their scientific questions while an AI assistant handles the technical complexity.

**The vision**: Every computational biologist has an intelligent assistant that understands both their biological research context and the technical tools required to execute it. This MCP server is the foundation for that future.

**Next steps**: Integration with Continue.dev for real-time IDE assistance, deployment in research environments, and expansion to other omics workflows.

Thank you!"

---

## **🎬 Video Direction Notes**

### **Timing & Transitions** (for background video)
- **0:45**: Transition to MCP server architecture diagram
- **1:30**: Switch to Gradio interface, focus on Enhanced Setup tab
- **2:15**: Move to Documentation tab, show indexing process
- **2:45**: Switch to Code Analysis tab, show example loading
- **3:00**: Transition to Log Debugging tab
- **3:30**: Move to MCP Tools tab, show tool selection
- **3:45**: Return to overview or architecture view

### **UI Elements to Highlight**
1. **Setup Tab**: API key fields, connection status, provider options
2. **Documentation Tab**: Index button, search bar, example searches
3. **Code Analysis Tab**: Code editor, context field, example buttons
4. **Log Debugging Tab**: Log input area, example log buttons
5. **MCP Tools Tab**: Tool dropdown, JSON input, quick examples

### **Key Visual Moments**
- **Setup**: Show connection success indicators
- **Documentation**: Demonstrate search results appearing
- **Code**: Show syntax highlighting and example code loading
- **Logs**: Display example error logs being pasted
- **Tools**: Show JSON responses from tool execution

---

## **🎯 Presentation Success Metrics**

**What the audience should understand**:
1. **Problem scope**: Spatial transcriptomics complexity
2. **Technical solution**: MCP protocol for tool integration
3. **Practical value**: Real workflow assistance vs simple chatbots
4. **Demo capabilities**: Each tab's purpose and functionality
5. **Future potential**: AI-assisted bioinformatics research

**Key takeaway**: This isn't just a demo - it's a working prototype of intelligent bioinformatics assistance that can be deployed in real research environments.
