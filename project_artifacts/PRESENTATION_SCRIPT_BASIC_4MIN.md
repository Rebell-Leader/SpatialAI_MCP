# 🎯 4-Minute Demo Presentation Script
**OpenProblems Spatial Transcriptomics MCP Server - Basic Demo**

---

## **[0:00-0:45] Project Introduction & Problem Statement**

"Hello! Today I'm presenting the **OpenProblems Spatial Transcriptomics MCP Server** - a solution addressing a critical challenge in computational biology.

**The Problem**: Computational biologists working with spatial transcriptomics face enormous complexity - managing 10-100x larger datasets, dealing with multiple tool ecosystems like Nextflow, Viash, and Docker, while trying to focus on their actual scientific research rather than wrestling with technical infrastructure.

**Our Solution**: We've developed a Model Context Protocol server that provides standardized access to bioinformatics tools through a unified interface. This allows AI agents and developers to interact with complex spatial transcriptomics workflows through simple, consistent commands rather than learning dozens of different tool syntaxes."

---

## **[0:45-1:30] MCP Server Capabilities**

"Our MCP server implements the full Model Context Protocol specification with **14 specialized tools** and **5 resource endpoints**:

**Core Infrastructure Tools**:
- **Environment checking** - Verify Nextflow, Viash, Docker, and Java installations
- **Nextflow workflow execution** - Run complex pipelines with proper parameter handling
- **Viash component management** - Build, test, and deploy workflow components
- **Docker operations** - Build images and manage containerized environments

**File & Configuration Tools**:
- **File operations** - Read, write, and navigate project structures
- **Configuration validation** - Check Nextflow configs and pipeline syntax
- **Log analysis** - Parse execution logs and identify common issues

**Spatial Transcriptomics Specialized Tools**:
- **Component template generation** - Create Viash components for spatial methods
- **Spatial data validation** - Verify SpatialData and AnnData format compliance
- **Environment setup** - Generate conda environments with spatial dependencies

The key innovation is **standardized access** - instead of remembering complex command-line syntaxes, users interact through simple MCP tool calls."

---

## **[1:30-2:15] Demo Interface Overview & Server Status**

"Let me walk you through our demonstration interface that showcases these MCP capabilities:

**Server Status Tab** - *[Video shows server overview]*
- **Real-time server information** showing the MCP server is running
- **Available tools list** - all 14 tools with descriptions
- **Resource endpoints** - 5 different information sources
- **Capability overview** - what the server can execute

**Server Resources Tab** - *[Video shows resource display]*
- **Server status** endpoint providing live system information
- **Documentation resources** for Nextflow, Viash, and Docker best practices
- **Spatial workflow templates** with pre-configured pipeline patterns
- **Resource content** displayed in structured format showing what's available to AI agents

This demonstrates how the MCP protocol exposes both **executable tools** and **information resources** through standardized endpoints."

---

## **[2:15-3:00] Tool Execution & Environment Management**

"**Tool Execution Tab** - *[Video shows tool selection and execution]*
- **Tool dropdown** with all 14 available MCP tools
- **Parameter input** in JSON format for tool arguments
- **Execution results** showing actual MCP server responses
- **Quick examples** for common operations:
  - **Environment Check** - Verify all required tools are installed
  - **Component Creation** - Generate new Viash component templates
  - **Log Analysis** - Parse Nextflow execution logs for errors

**File Operations Demo** - *[Video shows file management]*
- **Read file** tool demonstrating configuration file access
- **Write file** tool for creating new workflow components
- **Directory listing** for project navigation
- **Configuration validation** checking Nextflow pipeline syntax

These tools provide the **foundational operations** that AI agents need to assist with spatial transcriptomics workflows."

---

## **[3:00-3:45] Spatial Transcriptomics Integration**

"**Spatial Tools Demo** - *[Video shows spatial-specific functionality]*
- **Create Spatial Component** tool generating complete Viash component structures
  - Input: Component name and method type (segmentation, assignment, etc.)
  - Output: Complete directory structure with config.vsh.yaml, script.py, and documentation
- **Validate Spatial Data** tool checking SpatialData format compliance
- **Setup Spatial Environment** generating conda environment specifications

**Workflow Integration** - *[Video shows pipeline execution]*
- **Nextflow execution** tool running actual spatial transcriptomics pipelines
- **Parameter management** for complex workflow configurations
- **Docker integration** building and managing containerized environments
- **Log monitoring** tracking execution progress and identifying issues

This demonstrates how the MCP server bridges the gap between **high-level research goals** and **low-level technical implementation** in spatial transcriptomics."

---

## **[3:45-4:00] Production Deployment & Impact**

"This MCP server is **production-ready** and can be deployed in several ways:

**Continue.dev Integration**: Direct connection to IDEs for real-time development assistance
**Standalone Deployment**: Local server for research team collaboration
**API Integration**: RESTful access for custom applications

**The Impact**: Instead of researchers spending hours learning tool-specific syntaxes and debugging configuration issues, they can focus on their scientific questions while the MCP server handles the technical complexity.

**What makes this special**: This isn't just a wrapper around existing tools - it's a **standardized protocol implementation** that provides consistent, reliable access to the entire spatial transcriptomics toolchain.

**Next steps**: Deployment in research environments, integration with additional spatial analysis tools, and expansion to other omics workflows.

Thank you!"

---

## **🎬 Video Direction Notes**

### **Timing & Transitions** (for background video)
- **0:45**: Show MCP server architecture diagram with tool/resource overview
- **1:30**: Switch to Gradio interface, focus on Server Status tab
- **2:00**: Move to Server Resources tab, show resource content
- **2:15**: Switch to Tool Execution tab, show tool dropdown
- **2:45**: Demonstrate file operations tools
- **3:00**: Focus on spatial-specific tools (create component, validate data)
- **3:30**: Show workflow execution tools (Nextflow, Docker)
- **3:45**: Return to architecture overview or deployment diagram

### **UI Elements to Highlight**
1. **Server Status Tab**: Tool list, resource endpoints, capability overview
2. **Server Resources Tab**: Resource content display, documentation structure
3. **Tool Execution Tab**: Tool dropdown, JSON parameter input, execution results
4. **File Operations**: Read/write/list tools in action
5. **Spatial Tools**: Component creation output, data validation results
6. **Workflow Tools**: Nextflow execution, Docker operations

### **Key Visual Moments**
- **Tool List**: Show all 14 tools available
- **Resource Display**: Demonstrate structured information access
- **Tool Execution**: Show JSON input → structured output
- **Component Creation**: Display generated file structure
- **Log Analysis**: Show error detection and parsing
- **Environment Check**: Display system status validation

---

## **🎯 Presentation Success Metrics**

**What the audience should understand**:
1. **MCP Protocol Value**: Standardized access to complex bioinformatics tools
2. **Production Readiness**: Real server that can execute actual workflows
3. **Spatial Focus**: Specialized for spatial transcriptomics research needs
4. **Integration Potential**: Can connect to AI agents, IDEs, and custom applications
5. **Developer Experience**: Simplifies complex tool interactions

**Key takeaway**: This MCP server transforms fragmented bioinformatics toolchains into a unified, programmable interface that AI agents and developers can reliably use to assist spatial transcriptomics research.

---

## **🔧 Technical Demonstration Points**

**What the demo shows working**:
- **Real MCP protocol** implementation with proper tool/resource structure
- **Actual tool execution** with realistic parameter handling
- **Structured outputs** that AI agents can parse and use
- **Error handling** and validation for robust operation
- **Production deployment** readiness with proper configuration

**What makes it impressive**:
- **14 functional tools** covering the complete spatial transcriptomics workflow
- **5 resource endpoints** providing comprehensive documentation
- **Real file operations** and environment management
- **Standardized interface** that works across different platforms
- **Extensible architecture** for adding new tools and capabilities
