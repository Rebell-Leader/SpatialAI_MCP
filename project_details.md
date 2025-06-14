# **Model Context Protocol for Enhanced Spatial Transcriptomics Workflow Management in OpenProblems**

## **1\. Introduction: Bridging the Gap in Computational Biology Research**

Computational biology, particularly in the realm of single-cell and spatial transcriptomics, is experiencing an unprecedented surge in data complexity and analytical challenges. While researchers are primarily focused on developing novel scientific methods, the underlying computational infrastructure and auxiliary tools often present significant bottlenecks, diverting valuable scientific attention away from core biological questions. This report outlines a strategic approach to address this challenge within the OpenProblems project through the implementation of a Model Context Protocol (MCP) server, designed to streamline and standardize AI agent interaction with critical bioinformatics tools and data.

### **1.1 The OpenProblems Project: A Platform for Benchmarking Single-Cell Genomics**

The OpenProblems project stands as a pioneering initiative, characterized as a "living, extensible, community-guided benchmarking platform" dedicated to formalizing and evaluating open problems in single-cell genomics.1 This ambitious endeavor encompasses a wide array of critical tasks, including the preprocessing and rigorous evaluation of spatial transcriptomics simulators, as exemplified by task\_ist\_preprocessing and task\_spatial\_simulators.2 The project's robust benchmarking platform facilitates the contribution and standardized evaluation of containerized methods against well-defined datasets, leveraging the power of AWS Batch and Nextflow to ensure analyses are both scalable and reproducible.4 The underlying codebase of OpenProblems reflects a versatile, polyglot development environment, incorporating Shell, Python, Nextflow, and R, which highlights its adaptability and reliance on a diverse ecosystem of computational tools.1

### **1.2 The Bottleneck: Auxiliary Tools and Frameworks in Spatial Transcriptomics**

A central challenge articulated by the OpenProblems community is the tendency for computational biology researchers to prioritize the development of scientific methods themselves over the intermediate or auxiliary tools and frameworks essential for their practical implementation. This often results in a significant disconnect between innovative methodological advancements and their efficient, widespread application.

Spatial transcriptomics, while offering unparalleled insights into cellular interactions and tissue architecture, introduces formidable technical and computational hurdles.5 These include the management of exceptionally large datasets, which can be 10 to 100 times larger than those from single-cell RNA sequencing, frequently reaching terabytes per experimental run.5 Such data intensity demands substantial memory and processing power, often exceeding 128GB RAM and 32 CPU cores per sample, with processing times extending over several hours, rendering local analysis impractical for most researchers.5

Furthermore, ensuring reproducibility remains a significant challenge due to the diversity of platforms and computational workflows in spatial transcriptomics. Unlike single-cell RNA sequencing, the field currently lacks universally accepted computational pipelines, and the rapid evolution of analytical methods makes reliable replication difficult.5 Many researchers develop custom-built pipelines that often suffer from minimal documentation, severely impeding their reusability and broader adoption.5

Moreover, many existing software tools for transcriptomics analysis were originally designed for less complex data types and are not inherently equipped to handle the scale and intricacy of spatial transcriptomics data.6 This necessitates considerable manual effort for testing, validation, and the development of workarounds to ensure functionality and accuracy.6 The integration of multi-modal data and the need to bridge skills gaps between image processing and computational biology further compound these complexities.5

### **1.3 The Transformative Potential of AI Agents and the Model Context Protocol (MCP)**

The emergence of AI agents represents a pivotal shift in addressing these computational bottlenecks. AI agents are defined as autonomous software programs capable of interacting with their environment, collecting data, and performing self-determined tasks to achieve predefined goals.7 They are designed to execute complex, multi-step actions, learn, and adapt over time, making rational decisions based on their perceptions and available data.7 This capability is foundational to the burgeoning field of "agentic bioinformatics," which specifically deploys autonomous, adaptive, and intelligent AI agents to optimize, automate, and innovate biological data analysis workflows, thereby tackling complex biological challenges.9

The Model Context Protocol (MCP) serves as a critical open standard, developed by Anthropic, to standardize how AI applications—including custom agents—connect with external tools, data sources, and systems.10 It functions as a universal connector, enabling Large Language Models (LLMs) to dynamically interact with various APIs, databases, and business applications.11 This is particularly relevant given the observation that bioinformaticians prioritize methods over auxiliary tools. MCP's fundamental purpose is to standardize how AI agents interact with external tools and data, which extends beyond simple API calls. It establishes a structured, standardized interface that allows AI to "perceive environments, make decisions, and execute actions" within the intricate bioinformatics ecosystem.9

By abstracting the underlying complexities of tool integration, environment management, and workflow orchestration, the MCP server can act as a foundational layer, akin to a "Bioinformatics Operating System" for AI agents. This "OS" provides a standardized interface for AI applications to interact with computational resources and domain-specific software, enabling AI agents to operate at a higher, more conceptual level within bioinformatics. This paradigm suggests a transformative future where AI agents can more readily contribute to complex scientific domains beyond bioinformatics. By providing a universal, computable interface to domain-specific tools and data, it significantly lowers the barrier to entry for AI-driven scientific discovery and accelerates automation across diverse research fields.

## **2\. Foundational Technologies for Reproducible Bioinformatics**

The successful implementation of an MCP server for OpenProblems relies on a robust foundation of existing bioinformatics technologies. Nextflow, Viash, and Docker collectively provide the necessary framework for scalable, reproducible, and modular computational workflows.

### **2.1 Nextflow: A Robust Framework for Scalable Pipeline Orchestration**

Nextflow is recognized as a highly effective workflow framework, specifically engineered to enable bioinformaticians to integrate diverse scripts—including Bash, Python, Perl, and R—into cohesive, portable, reproducible, scalable, and checkpointed pipelines.12 Its inherent support for containerization technologies like Docker and Singularity ensures the consistent reproducibility of analyses across different environments.12 The framework’s ability to execute workflows seamlessly across various computational infrastructures, ranging from local machines to High-Performance Computing (HPC) clusters (e.g., Slurm, SGE, PBS) and cloud platforms (e.g., Google, Kubernetes, AWS), guarantees exceptional portability and scalability.12 The OpenProblems project already leverages Nextflow extensively for its benchmarking efforts on AWS, highlighting its proven utility in large-scale scientific endeavors.4

Key features of Nextflow that contribute to its efficacy include rapid prototyping, which allows for quick development of computational pipelines from smaller tasks. It also offers efficient unified parallelism, achieved by sharding data and submitting each shard as a separate job, particularly beneficial for single-threaded tools.12 Furthermore, its continuous checkpointing mechanism allows for seamless resumption of pipeline execution from the last successfully completed step, even in the event of failures, thereby enhancing robustness and efficiency.12 For optimizing large-scale pipelines, best practices involve minimizing data transfer between steps, enhancing I/O performance by co-locating data with compute resources, and strategically utilizing scalable storage options such as Amazon S3.15 Nextflow also provides robust error handling mechanisms, including errorStrategy directives (ignore, retry) and maxRetries for managing transient conditions, alongside capabilities for dynamic resource allocation based on task attempts, which can prevent out-of-memory errors and other common issues.16

### **2.2 Viash: Modularizing and Standardizing Bioinformatics Components**

Viash is an open-source meta-framework that directly addresses the prevalent challenge of tightly coupled software components in bioinformatics workflows. It actively promotes reusability and significantly reduces maintenance overhead by decoupling component functionality from workflow logic.18 This design principle allows developers to focus on implementing the core functionality of a tool without needing expert knowledge of specific workflow frameworks like Nextflow or cloud environments.18

Viash facilitates a "code-first" prototyping approach: users write a core script and add minimal metadata in a YAML configuration file. From this, Viash automatically generates boilerplate code for modular Nextflow components, standalone executables with auto-generated command-line interfaces (CLIs), and Docker images.18 This automation significantly speeds up development and reduces time spent on repetitive coding tasks.

The transformation of a simple script and metadata into various deployable artifacts—Docker images, Nextflow modules, and standalone executables—positions Viash as a crucial "compiler" for MCP-ready bioinformatics components. It automates the generation of CLIs, documentation, and enforces best practices such as versioning and robust argument validation.19 The MCP specification mandates that servers "wrap external capabilities according to the MCP specification".10 For AI agents to effectively utilize bioinformatics tools, these tools must be standardized and consistently packaged. Viash directly addresses this by acting as a critical factory that translates human-written bioinformatics logic into standardized, containerized, and well-documented components. These components are then inherently ready to be exposed as "Tools" by the MCP server, significantly streamlining the process of creating MCP-compatible bioinformatics operations. This automated generation of standardized, containerized components directly reduces the manual effort and potential for errors in preparing bioinformatics tools for MCP integration, thereby accelerating the development and deployment of AI-driven bioinformatics solutions within the OpenProblems project. This directly addresses the need to abstract away auxiliary tool complexities.

Viash further enhances reproducibility through automated versioning of artifacts, intelligent argument parsing and validation, and seamless integration with containerization technologies (Docker) and Continuous Integration (CI) tools like GitHub Actions and Jenkins.18 Its polyglot support for Bash, Python, R, Docker, and Nextflow makes it exceptionally well-suited for the diverse technological landscape of bioinformatics.18 Data Intuitive, a key contributor to Viash, offers the Viash Catalogue, an extensive collection of over 150 industry-ready, open-source bioinformatics workflows and tools, including specialized solutions for single-cell transcriptomics, further exemplifying the framework's utility.22

### **2.3 Docker: Ensuring Consistent and Portable Computational Environments**

Docker is an indispensable technology for deploying bioinformatics applications and analysis pipelines, providing a consistent and reproducible operating environment by encapsulating software and all its dependencies within isolated containers.24 This containerization approach enables the isolation, capture, reuse, and sharing of computational environments, which is paramount for large-scale analyses that involve numerous tools and diverse programming languages.25

Dockerfiles serve as explicit blueprints, defining the step-by-step instructions for building a Docker image. These instructions include commands such as FROM (specifying the base image), RUN (executing shell commands), COPY (transferring data from host to image), ENTRYPOINT (setting the command to be run when a container is created), and WORKDIR (setting the current working directory).24 Best practices for Dockerfile creation include implementing multi-stage builds for improved caching and combining apt-get update && install commands into a single layer to prevent caching issues and reduce image size.24

Nextflow extensively supports Docker, facilitating the creation of scalable and reproducible scientific workflows that leverage containerization for robust dependency management and environment consistency.4 Fundamental Docker commands such as docker run (to create and start a container), docker ps (to list running containers), docker stop (to stop a container), docker rm (to remove a container), docker images (to list images), and docker rmi (to remove an image) are essential for effective management of containers and images throughout the development and deployment lifecycle.24

## **3\. The Model Context Protocol (MCP): A Standard for AI-Tool Interaction**

The Model Context Protocol (MCP) is central to enabling AI agents to interact effectively and intelligently with complex bioinformatics workflows and data. It provides the necessary standardization and structure for seamless communication.

### **3.1 Core Concepts of MCP: Tools, Resources, and Communication Mechanisms**

The Model Context Protocol (MCP) is an open standard, primarily championed by Anthropic, designed to standardize how AI applications seamlessly connect with external tools, data sources, and systems.10 It operates on a client-server architecture and employs JSON-RPC as its underlying communication protocol.29

Within the MCP architecture, several key roles are defined:

* **Hosts:** These represent the user-facing applications, such as Claude Desktop, Integrated Development Environments (IDEs) like Cursor, or custom AI agents, which manage the overall communication flow with MCP servers.10  
* **Clients:** Embedded within Host applications, clients are responsible for managing connections, discovering available capabilities, forwarding requests, and handling responses from specific MCP servers.10  
* **Servers:** These are the crucial bridge or API components. MCP servers expose the specific functionalities of external systems—such as APIs, databases, or local files—by wrapping them according to the MCP specification.10 Servers can be built in various programming languages, provided they can communicate over the supported transports.

MCP defines fundamental primitives that govern how AI agents interact with external capabilities:

* **Tools (Model-controlled):** These represent functions or actions that Large Language Models (LLMs) can invoke to perform specific operations, akin to function calling mechanisms. An example is a weather API, where the AI decides to call the function to retrieve data.10  
* **Resources (Application-controlled):** These are data sources that LLMs can access to retrieve contextual information. They function similarly to GET endpoints in a REST API, providing data without initiating significant computation or side effects. Resources are considered part of the context or request provided to the AI.10  
* **Prompts (User-controlled):** These are predefined templates or instructions that are triggered by user actions, guiding the AI's initial interaction or task.30

Communication between MCP servers and clients primarily occurs through two robust methods:

* **stdio (Standard Input/Output):** This method is employed when the Client and Server are running on the same machine. It is simple and effective for local integrations, such as accessing local files or running a local script.10  
* **HTTP via SSE (Server-Sent Events):** For persistent connections, the Client connects to the Server using HTTP. After an initial setup, the Server can push messages (events) to the Client over this persistent connection, utilizing the Server-Sent Events standard.10

### **3.2 MCP's Role in Enabling Intelligent Bioinformatics Agents**

MCP refines existing patterns in AI agent development by clearly delineating between "Tools" (actions the AI decides to take) and "Resources" (contextual information provided to the AI), thereby enhancing clarity and control over AI interactions.10 This structured approach provides a standardized pathway for AI models to dynamically interact with APIs, databases, and other applications. Such standardization ensures consistent AI integration, offers flexibility (allowing easy switching between different AI models and vendors), and maintains robust security by keeping data within the user's infrastructure.11

AI agents, powered by sophisticated LLMs, are capable of processing multimodal information, performing complex reasoning, learning, and making informed decisions.8 Their effectiveness is significantly amplified by standardized access to external tools and data through MCP. While current AI models, such as GPT-4o and Claude 3.5 Sonnet, still exhibit limitations in performing complex, iterative bioinformatics tasks—for example, accurately interpreting intricate plots, managing diverse data formats, and achieving only approximately 17% accuracy on open-answer tasks in some benchmarks—MCP provides the essential structured interface to mitigate these challenges by externalizing tool usage and context provision.31 This externalization allows the AI to focus on higher-level reasoning and problem-solving, rather than the intricacies of tool invocation and data formatting.

### **3.3 Synergistic Integration: MCP with Nextflow, Viash, and Docker**

The Model Context Protocol serves as the crucial connector, linking AI agents to their necessary tools and knowledge.29 This is precisely where the robust capabilities of Nextflow, Viash, and Docker become indispensable, creating a powerful synergy for bioinformatics research.

Viash components, inherently designed for modularity, standardization, and containerization (Docker), are ideally suited to be directly exposed as MCP "Tools." This leverages Viash's automated code generation capabilities for CLIs, Docker images, and Nextflow modules, ensuring that bioinformatics operations are readily consumable by AI agents.18 Higher-level Nextflow pipelines, which orchestrate these individual Viash/Docker components, can also be exposed as MCP "Tools." This enables AI agents to initiate, monitor, and manage complex, multi-step bioinformatics workflows with a single, standardized command.12 Docker containers play a critical role in ensuring that the execution environment for any tool or pipeline invoked via MCP is entirely consistent and reproducible, irrespective of the underlying computational infrastructure.24

Spatial transcriptomics data and all associated metadata (e.g., adhering to the CELLxGENE schema 3) can be exposed as MCP "Resources." This provides the essential contextual information that AI agents require to accurately understand, analyze, and interpret complex biological data, directly addressing the need for real-time, structured operational data for AI agents.32

This integration fundamentally transforms the role of AI agents into "Cognitive Accelerators" for spatial transcriptomics. Spatial transcriptomics faces formidable challenges related to data scale, reproducibility, multi-modal integration, and specialized skill requirements.5 AI agents, particularly those driven by LLMs, exhibit strong capabilities in areas like pattern recognition, predictive modeling, data preprocessing, and visualization.34 However, they often struggle with the iterative, exploratory, and subjective aspects inherent in bioinformatics analysis.27 By integrating Nextflow, Viash, and Docker through the MCP, the AI agent is liberated from managing the low-level complexities of tool installation, environment setup, or intricate workflow execution. Instead, it interacts with standardized "Tools" (e.g., a Viash-generated Nextflow module for spatial data normalization) and retrieves structured "Resources" (e.g., an AnnData object with spatial coordinates). This abstraction allows the AI's sophisticated capabilities to be focused on the higher-level scientific problems, such as identifying spatially variable genes or integrating multi-modal data. This approach fundamentally shifts the role of AI agents from mere code generators to powerful augmentations for human bioinformaticians. They handle the computationally intensive, repetitive, and infrastructure-heavy aspects of data analysis, freeing human researchers to concentrate on hypothesis generation, deep biological interpretation, and novel method development. This also underscores the continued necessity for human oversight and refinement, particularly for subjective analytical steps where AI currently lacks nuanced understanding.27

## **4\. MCP for OpenProblems: Revolutionizing Spatial Transcriptomics Workflows**

The MCP server will be established as a central hub within the OpenProblems project, providing a standardized and machine-readable interface for AI agents to interact with the computational environment, with a specific focus on spatial transcriptomics. This strategic implementation will significantly enhance the efficiency and reproducibility of spatial transcriptomics tool development, evaluation, and benchmarking.

### **4.1 Strategic Impact Areas of the MCP Server for Scientists**

The MCP server will address several critical areas to empower bioinformaticians and accelerate scientific discovery within the OpenProblems project:

#### **4.1.1 Centralized and Contextualized Documentation for Key Tools**

**Current Challenge:** Bioinformatics tools, particularly custom-built pipelines, frequently suffer from minimal, outdated, or rapidly changing documentation, which severely hinders reproducibility and comprehension.5 The existence of disparate documentation sources for Docker, Viash, and Nextflow further complicates the learning curve for researchers.

**MCP-Enabled Solution:** The MCP server will expose comprehensive, machine-readable documentation for all integrated tools (Nextflow pipelines, Viash components, Docker images) as structured "Resources".10 This documentation will include detailed parameter schemas, practical usage examples, and adherence to best practices, all directly accessible by AI agents and human users through a standardized interface. This approach transforms static, disparate documentation into a computable, queryable "knowledge graph." AI agents require structured data to make informed decisions.32 By exposing not only raw data but also the metadata and functional specifications of Nextflow pipelines, Viash components, and Docker images as MCP Resources, the MCP server enables AI agents to understand the relationships between tools, their inputs/outputs, and their scientific purpose. This allows for a deeper, more active understanding of the bioinformatics ecosystem that goes beyond simple information retrieval, representing a significant advancement over traditional documentation by enabling dynamic, context-aware interaction. This structured, machine-readable documentation and metadata exposed via MCP Resources enables AI agents to build a richer, more actionable understanding of the bioinformatics domain, which, in turn, leads to more effective tool invocation, precise parameter selection, and overall improved problem-solving, directly addressing the critical need for context for coding agents.

#### **4.1.2 Empowering Context-Aware AI Coding Agents for Workflow Development**

**Current Challenge:** Existing AI models often struggle with the specific nuances of Nextflow, for instance, defaulting to DSL1 instead of DSL2, necessitating substantial debugging and validation efforts from human researchers.36 Furthermore, integrating diverse data formats and accurately interpreting complex plots remain significant hurdles for AI agents.31

**MCP-Enabled Solution:** AI coding agents, interacting directly via the MCP server, will gain privileged access to the latest Nextflow, Viash, and Docker best practices, along with structured schemas, all exposed as MCP Resources. This rich context will enable them to generate DSL2-compliant Nextflow code, precise Viash component configurations, and optimized Dockerfiles that inherently adhere to OpenProblems' stringent standards, with integrated testing capabilities.36 This direct access to structured information and best practices, facilitated by MCP, significantly enhances the AI's ability to generate accurate and functional bioinformatics code.

#### **4.1.3 Enforcing Best Practices and Standardized Guidelines**

**Current Challenge:** The absence of universally accepted computational pipelines in spatial transcriptomics contributes significantly to reproducibility issues, and many custom pipelines lack consistent standardization across research groups.5

**MCP-Enabled Solution:** The MCP server will function as a central gatekeeper and enforcer of best practices within the OpenProblems ecosystem. By defining all tools and resources with strict MCP schemas, it ensures that all interactions and generated components automatically adhere to predefined standards for reproducibility, scalability, and maintainability.19 This encompasses detailed guidelines for Dockerfile optimization 24, Nextflow resource tuning 15, and Viash modularity principles 19, aligning with OpenProblems' core mission of formalizing and benchmarking.1 This enforcement through standardized interfaces ensures that all contributions to the OpenProblems project meet a consistent level of quality and reproducibility.

#### **4.1.4 Providing Curated Examples and Reusable Pipeline Templates**

**Current Challenge:** Researchers often resort to developing in-house workflows with minimal documentation, making it challenging to replicate results or share methods effectively.5 Building complex bioinformatics pipelines from scratch is a time-consuming and error-prone endeavor.

**MCP-Enabled Solution:** The MCP server will expose a meticulously curated library of Nextflow pipeline templates (e.g., for spatial transcriptomics basic processing, identification of spatially variable genes, and label transfer, as seen in SpatialNF 38) and Viash component examples (leveraging the Viash Catalogue 22) as easily discoverable and consumable MCP Resources. AI agents can then leverage these templates to rapidly prototype new workflows, significantly accelerating development cycles and ensuring consistency across projects. This direct access to pre-validated and standardized templates reduces the need for researchers to start from scratch, fostering a more collaborative and efficient development environment.

#### **4.1.5 Facilitating Comprehensive Implementation Checklists**

**Current Challenge:** The inherent complexity of integrating multiple sophisticated tools and frameworks—such as Nextflow, Viash, and Docker—can lead to overlooked steps, configuration errors, and significant delays during implementation.

**MCP-Enabled Solution:** The MCP server can provide AI agents with direct access to structured implementation checklists, exposed as MCP Resources.10 These checklists will guide the AI through the systematic setup, configuration, and deployment of new workflows or components. Critically, these checklists can be dynamically updated and validated by the AI agent itself, ensuring strict adherence to OpenProblems' evolving standards and reducing human oversight requirements. This capability allows the AI agent to perform complex, multi-step actions with greater accuracy and completeness 8, minimizing human error in complex setup procedures.

#### **4.1.6 Streamlining Testing and Advanced Troubleshooting**

**Current Challenge:** Reproducibility remains a significant hurdle in spatial transcriptomics due to platform variability and the rapid evolution of analytical standards.5 Debugging complex Nextflow pipelines is often challenging, requiring laborious manual inspection of work directories and log files.16

**MCP-Enabled Solution:** The MCP server will expose specialized "Tools" for automated testing (e.g., generating and executing nf-test scripts 36; running Viash unit tests 18) and advanced troubleshooting (e.g., analyzing Nextflow logs for actionable insights, identifying common errors like Out-Of-Memory (OOM) issues, and suggesting dynamic resource allocation 16). This enables AI-driven "Proactive Troubleshooting" and "Test-Driven Workflow Development." Nextflow provides detailed error reporting 16, and Seqera AI can analyze these logs to provide actionable insights.37 Furthermore, Seqera AI can generate nf-test scripts and offers "one-click testing in an AI sandbox" with self-correction capabilities.36 By exposing these functionalities as MCP Tools, AI agents can transcend reactive debugging. They can proactively initiate tests (e.g., before deployment or after code changes), continuously monitor pipeline execution for anomalies, diagnose errors by analyzing logs (e.g., OOM errors, missing commands 16), and even suggest or implement dynamic resource adjustments or code fixes. This capability significantly enhances the robustness and reliability of bioinformatics workflows by automating error detection and resolution, thereby accelerating the development and validation cycle.

## **5\. Detailed MCP Project Description for OpenProblems**

The Model Context Protocol (MCP) server for OpenProblems will serve as a central, standardized interface, enabling AI agents to interact intelligently with the complex ecosystem of Nextflow pipelines, Viash components, Dockerized workflows, and spatial transcriptomics data. This server will adhere to the MCP specification, exposing capabilities as "Tools" and contextual information as "Resources."

**Project Name:** OpenProblems Spatial Transcriptomics MCP Server

**Purpose:** To provide a standardized, machine-readable interface for AI agents to interact with Nextflow pipelines, sc/spatial transcriptomics data processing methods, and Viash-managed dockerized workflows within the OpenProblems project, thereby abstracting auxiliary tool complexities and enabling bioinformaticians to focus on scientific innovation.

**Target Audience:** AI agents (e.g., LLM-driven coding assistants, autonomous research agents), bioinformaticians, computational biologists, and developers contributing to the OpenProblems project.

**Core Functionality (Exposed via MCP Primitives):**

**5.1 MCP Tools (Model-controlled actions):**

* **Nextflow Workflow Execution:**  
  * **Tool Name:** run\_nextflow\_workflow  
  * **Description:** Executes a specified Nextflow pipeline from the OpenProblems or OpenPipelines-bio repositories.  
  * **Parameters:**  
    * workflow\_name: (string, required) Name of the Nextflow workflow (e.g., task\_ist\_preprocessing/main.nf, openpipeline/main.nf, SpatialNF/main.nf).  
    * github\_repo\_url: (string, required) GitHub URL of the repository containing the workflow (e.g., https://github.com/openproblems-bio/task\_ist\_preprocessing).  
    * profile: (string, optional) Nextflow profile to use (e.g., docker, singularity, test).  
    * params: (JSON object, optional) Key-value pairs for Nextflow pipeline parameters (e.g., {"input\_file": "data.h5ad", "output\_dir": "results"}).  
    * config\_file: (string, optional) Path to a custom Nextflow configuration file.  
  * **Output:** Execution ID, link to Nextflow log, status (running, completed, failed).  
* **Viash Component Execution:**  
  * **Tool Name:** run\_viash\_component  
  * **Description:** Executes a specific Viash component, either as a standalone executable or within a Docker container.  
  * **Parameters:**  
    * component\_name: (string, required) Name of the Viash component (e.g., process\_dataset, metric).  
    * component\_config\_path: (string, required) Path to the Viash config file (.vsh.yaml).  
    * engine: (string, optional, default: docker) Execution engine (native, docker).  
    * args: (JSON object, optional) Key-value pairs for component-specific arguments (e.g., {"input\_sc": "sc.h5ad", "output\_sp": "sp.h5ad"}).  
  * **Output:** Execution ID, link to component logs, output file paths, status.  
* **Dockerized Workflow Building:**  
  * **Tool Name:** build\_docker\_image  
  * **Description:** Builds a Docker image from a specified Dockerfile path.  
  * **Parameters:**  
    * dockerfile\_path: (string, required) Path to the Dockerfile.  
    * image\_tag: (string, required) Tag for the Docker image (e.g., openproblems/spatial-tool:1.0.0).  
    * context\_path: (string, optional, default: .) Build context directory.  
  * **Output:** Docker image ID, build logs, status.  
* **Automated Testing:**  
  * **Tool Name:** run\_nf\_test  
  * **Description:** Generates and executes nf-test scripts for a given Nextflow pipeline or Viash component.  
  * **Parameters:**  
    * pipeline\_path: (string, required) Path to the Nextflow pipeline or Viash component.  
    * test\_scope: (string, optional, default: all) Scope of tests to run (e.g., unit, integration, all).  
  * **Output:** Test report, pass/fail status, log of test execution.  
* **Log Analysis & Troubleshooting:**  
  * **Tool Name:** analyze\_nextflow\_log  
  * **Description:** Analyzes a Nextflow execution log to identify errors, suggest causes, and provide actionable insights.  
  * **Parameters:**  
    * log\_file\_path: (string, required) Path to the .nextflow.log file.  
  * **Output:** Structured error report (JSON), suggested troubleshooting steps, potential fixes (e.g., memory adjustments, command corrections).

**5.2 MCP Resources (Application-controlled context):**

* **Documentation Context:**  
  * **Resource Name:** documentation\_context://{tool\_name}  
  * **Description:** Provides structured, machine-readable documentation for Nextflow, Viash, Docker, and specific OpenProblems tools/pipelines.  
  * **Content:** Parameter schemas (JSON Schema), usage examples, best practices guidelines (e.g., Dockerfile optimization, Nextflow resource tuning), common errors and their resolutions, versioning information.  
* **Pipeline Templates:**  
  * **Resource Name:** pipeline\_template://{template\_id}  
  * **Description:** Access to curated Nextflow pipeline templates and Viash component examples for spatial transcriptomics.  
  * **Content:** Workflow definition files (.nf), Viash config files (.vsh.yaml), example input data paths, READMEs.  
* **Implementation Checklists:**  
  * **Resource Name:** implementation\_checklist://{checklist\_id}  
  * **Description:** Structured checklists for setting up, configuring, and deploying new workflows or components.  
  * **Content:** Step-by-step instructions, required dependencies, configuration parameters, validation criteria.  
* **Spatial Transcriptomics Data Access:**  
  * **Resource Name:** spatial\_data://{dataset\_id}  
  * **Description:** Provides access to preprocessed spatial transcriptomics datasets and associated metadata.  
  * **Content:** File paths to AnnData objects (.h5ad) containing raw counts and metadata (CELLxGENE schema v4.0.0), spatial coordinates, relevant experimental metadata.

**Communication Methods:**

* **Primary:** stdio for local development and testing environments where AI agents run on the same machine as the MCP server.  
* **Secondary:** HTTP via SSE for remote deployments, allowing persistent connections and event streaming for monitoring long-running tasks.

**Technology Stack:**

* **Server Implementation:** Python (using fastmcp or similar SDK for rapid development).  
* **Orchestration:** Nextflow.  
* **Containerization:** Docker.  
* **Component Framework:** Viash.  
* **Data Formats:** AnnData (.h5ad), JSON, YAML, plain text for logs.

## **6\. Implementation Instructions for DEV AI Agent**

The following detailed list of tasks outlines the implementation roadmap for a Development AI Agent responsible for building and integrating the OpenProblems Spatial Transcriptomics MCP Server.

**Phase 1: Environment Setup and Core MCP Server Development**

1. **Initialize Project Repository:**  
   * Create a new GitHub repository for the MCP server (e.g., openproblems-mcp-server).  
   * Set up basic project structure: src/, config/, docs/, tests/, Docker/.  
2. **Set Up Python Environment:**  
   * Create a Python virtual environment.  
   * Install fastmcp (or chosen MCP SDK) and other core dependencies (e.g., pyyaml, requests, nextflow).  
3. **Develop Core MCP Server Application:**  
   * Implement the main MCP server application in src/main.py.  
   * Define the FastMCP instance with a descriptive name (e.g., OpenProblemsBioMCP).  
4. **Implement Basic MCP Tools:**  
   * **echo\_test Tool:** Create a simple @mcp.tool() function that echoes input, to verify basic MCP communication.  
   * **list\_available\_tools Tool:** Implement a tool that dynamically lists all registered MCP tools and their descriptions.  
5. **Implement Basic MCP Resources:**  
   * **server\_status Resource:** Create an @mcp.resource() that returns the server's current status and version.  
   * **read\_file Resource:** Implement a resource that can read and return the content of a specified local file (e.g., README.md).  
6. **Containerize the MCP Server:**  
   * Create a Dockerfile for the MCP server, including Python, fastmcp, and other dependencies.  
   * Ensure the Dockerfile is optimized for size and build time (e.g., multi-stage build, apt-get update && install in single layer).  
   * Build and test the Docker image locally.

**Phase 2: Integrating Foundational Bioinformatics Technologies**

1. **Integrate Nextflow Execution Tool:**  
   * **Tool Name:** run\_nextflow\_workflow  
   * **Implementation:**  
     * The tool will accept workflow\_name, github\_repo\_url, profile, params, and config\_file.  
     * Use subprocess to execute nextflow run {github\_repo\_url}/{workflow\_name} \-profile {profile} \--{params} \-c {config\_file}.  
     * Capture stdout, stderr, and exit code.  
     * Return a unique execution ID and paths to generated log files.  
   * **Error Handling:** Implement Nextflow's errorStrategy and maxRetries logic within the tool's execution for robustness.  
2. **Integrate Viash Component Execution Tool:**  
   * **Tool Name:** run\_viash\_component  
   * **Implementation:**  
     * The tool will accept component\_name, component\_config\_path, engine, and args.  
     * Execute viash run {component\_config\_path} \-p {engine} \-- {args} via subprocess.  
     * Parse Viash's output to identify output file paths and execution status.  
   * **Dependency Management:** Ensure the Docker image for the MCP server includes Viash or that Viash is run within its own container via the tool.  
3. **Integrate Docker Image Building Tool:**  
   * **Tool Name:** build\_docker\_image  
   * **Implementation:**  
     * The tool will accept dockerfile\_path, image\_tag, and context\_path.  
     * Execute docker build \-t {image\_tag} {context\_path} via subprocess.  
     * Capture build logs and return the resulting Docker image ID.  
   * **Best Practices Enforcement:** Automatically check for common Dockerfile best practices (e.g., apt-get update && install in one layer, multi-stage builds) and provide warnings or suggestions as part of the output.  
4. **Develop Data Access Resources:**  
   * **Resource Name:** spatial\_data://{dataset\_id}  
   * **Implementation:**  
     * The resource will map dataset\_id to predefined paths for h5ad files.  
     * Return the file path and relevant metadata (e.g., CELLxGENE schema version, organism, assay type) for the specified spatial transcriptomics dataset.  
     * Ensure secure access control if sensitive data is involved.

**Phase 3: Advanced Features and Documentation**

1. **Implement Automated Testing Tool (run\_nf\_test):**  
   * **Tool Name:** run\_nf\_test  
   * **Implementation:**  
     * The tool will accept pipeline\_path and test\_scope.  
     * Execute nf-test test {pipeline\_path} \--profile {test\_scope}.  
     * Parse nf-test output to generate a structured test report (JSON) indicating pass/fail status and details of failed tests.  
2. **Implement Log Analysis and Troubleshooting Tool (analyze\_nextflow\_log):**  
   * **Tool Name:** analyze\_nextflow\_log  
   * **Implementation:**  
     * The tool will accept log\_file\_path.  
     * Parse the Nextflow log file (.nextflow.log, .command.err, .command.out) to identify error patterns (e.g., exit status 137 for OOM, "command not found").  
     * Use rule-based logic or a small, fine-tuned LLM (if available and feasible) to suggest specific troubleshooting steps (e.g., increase memory, install missing software, check file paths).  
     * Return a structured report of identified issues and suggested actions.  
3. **Develop Comprehensive Documentation Resources:**  
   * **Resource Name:** documentation\_context://{tool\_name}  
   * **Implementation:**  
     * For each implemented MCP Tool and Resource, create a corresponding structured documentation entry.  
     * Define JSON schemas for all tool parameters and resource outputs.  
     * Provide markdown-formatted usage examples for each tool and resource.  
     * Include sections on best practices for Nextflow, Viash, and Docker relevant to OpenProblems.  
     * Ensure this documentation is dynamically loadable by the MCP server.  
4. **Curate Pipeline Templates and Examples Resource:**  
   * **Resource Name:** pipeline\_template://{template\_id}  
   * **Implementation:**  
     * Identify key Nextflow pipelines from openproblems-bio/task\_ist\_preprocessing, openpipelines-bio/openpipeline, and SpatialNF that serve as valuable templates.  
     * Create structured metadata for each template (description, inputs, outputs, relevant use cases).  
     * Expose the raw .nf and .vsh.yaml files, along with example input data paths, as part of this resource.  
5. **Develop Implementation Checklists Resource:**  
   * **Resource Name:** implementation\_checklist://{checklist\_id}  
   * **Implementation:**  
     * Create structured checklists for common tasks:  
       * "New Nextflow Pipeline Integration Checklist"  
       * "New Viash Component Development Checklist"  
       * "Docker Image Optimization Checklist"  
     * Each checklist item should include a description, a pass/fail criterion, and suggested actions.  
     * Expose these checklists as MCP Resources.

**Phase 4: Testing, Deployment, and Maintenance**

1. **Unit and Integration Testing:**  
   * Write unit tests for each MCP Tool and Resource function.  
   * Develop integration tests to verify the end-to-end functionality of AI agents interacting with the MCP server and underlying bioinformatics tools.  
   * Automate testing using GitHub Actions or a similar CI/CD pipeline.  
2. **Deployment Strategy:**  
   * Define deployment procedures for the MCP server (e.g., Docker Compose for local/on-prem, Kubernetes for cloud).  
   * Ensure the server can be deployed securely and with appropriate access controls.  
3. **Monitoring and Logging:**  
   * Implement robust logging for all MCP server interactions and tool executions.  
   * Integrate with monitoring tools to track server health, performance, and error rates.  
4. **Continuous Improvement:**  
   * Establish a feedback loop for AI agent performance and user experience.  
   * Regularly update MCP Tools and Resources to reflect new versions of Nextflow, Viash, Docker, and evolving best practices in spatial transcriptomics.  
   * Expand the library of pipeline templates and documentation based on community needs.

## **7\. Conclusions and Recommendations**

The implementation of a Model Context Protocol (MCP) server within the OpenProblems project represents a pivotal step towards revolutionizing spatial transcriptomics workflows. By providing a standardized, machine-readable interface, the MCP server will abstract away the complexities of auxiliary tools and frameworks, allowing bioinformaticians to dedicate their focus to scientific innovation. This approach transforms the current landscape by enabling AI agents to act as "Bioinformatics Operating Systems," providing a universal, computable interface to domain-specific tools and data, thereby lowering the barrier to entry for AI-driven scientific discovery.

The synergistic integration of MCP with Nextflow, Viash, and Docker facilitates the creation of "Cognitive Accelerators" in the form of AI agents. These agents, liberated from low-level computational complexities, can concentrate on higher-level scientific problems such as identifying spatially variable genes, integrating multi-modal data, and performing complex analyses with unprecedented efficiency. Furthermore, the MCP server will function as a "Knowledge Graph Interface" for bioinformatics, converting disparate documentation into computable resources that AI agents can actively query and understand. This will also enable AI-driven "Proactive Troubleshooting" and "Test-Driven Workflow Development," where AI agents can automatically initiate tests, diagnose issues, and even suggest or implement fixes, significantly enhancing the robustness and reliability of bioinformatics pipelines.

**Recommendations for OpenProblems Project:**

1. **Prioritize MCP Server Development:** Allocate dedicated resources to the development and maintenance of the OpenProblems Spatial Transcriptomics MCP Server as outlined in this report. This server is foundational to integrating AI agents effectively.  
2. **Standardize Tool Exposure:** Ensure all existing and new bioinformatics tools and pipelines within OpenProblems are wrapped as Viash components, making them inherently compatible for exposure as MCP "Tools." This will maximize reusability and standardization.  
3. **Invest in Structured Documentation:** Develop and maintain comprehensive, machine-readable documentation (e.g., JSON schemas, usage examples) for all tools and datasets, accessible as MCP "Resources." This is critical for enabling AI agents to understand and effectively utilize the bioinformatics ecosystem.  
4. **Foster AI Agent Integration:** Actively encourage the development and integration of AI agents (e.g., LLM-driven coding assistants, automated analysis agents) that leverage the MCP server. Provide clear guidelines and examples for agent developers.  
5. **Establish Continuous Feedback and Improvement:** Implement mechanisms for collecting feedback on the MCP server's performance and utility from both human users and AI agents. Continuously refine the MCP implementation, tools, and resources based on evolving research needs and technological advancements.  
6. **Promote Community Contribution:** Leverage the open-source nature of MCP, Nextflow, Viash, and Docker to foster community contributions to the MCP server, its tools, and associated documentation, aligning with the OpenProblems project's community-guided mission.

By embracing the Model Context Protocol, OpenProblems can significantly enhance the efficiency, reproducibility, and accessibility of spatial transcriptomics research, empowering bioinformaticians to push the boundaries of biological discovery.

#### **Источники**

1. openproblems-bio/openproblems: Formalizing and ... \- GitHub, дата последнего обращения: мая 28, 2025, [https://github.com/openproblems-bio/openproblems](https://github.com/openproblems-bio/openproblems)  
2. дата последнего обращения: января 1, 1970, [https://github.com/openproblems-bio/task\_ist\_preprocessing](https://github.com/openproblems-bio/task_ist_preprocessing)  
3. openproblems-bio/task\_spatial\_simulators: Benchmarking ... \- GitHub, дата последнего обращения: мая 28, 2025, [https://github.com/openproblems-bio/task\_spatial\_simulators](https://github.com/openproblems-bio/task_spatial_simulators)  
4. Driving innovation in single-cell analysis on AWS | AWS Public Sector Blog, дата последнего обращения: мая 28, 2025, [https://aws.amazon.com/blogs/publicsector/driving-innovation-single-cell-analysis-aws/](https://aws.amazon.com/blogs/publicsector/driving-innovation-single-cell-analysis-aws/)  
5. Spatial Transcriptomics at Scale: How to Overcome the Top 5 Data ..., дата последнего обращения: мая 28, 2025, [https://www.viascientific.com/blogs/spatial-transcriptomics-at-scale-how-to-overcome-the-top-5-data-hurdles](https://www.viascientific.com/blogs/spatial-transcriptomics-at-scale-how-to-overcome-the-top-5-data-hurdles)  
6. From bulk to spatial: How transcriptomics is changing the way we see biology \- Ardigen, дата последнего обращения: мая 28, 2025, [https://ardigen.com/from-bulk-to-spatial-how-transcriptomics-is-changing-the-way-we-see-biology/](https://ardigen.com/from-bulk-to-spatial-how-transcriptomics-is-changing-the-way-we-see-biology/)  
7. What are AI Agents?- Agents in Artificial Intelligence Explained \- AWS, дата последнего обращения: мая 28, 2025, [https://aws.amazon.com/what-is/ai-agents/](https://aws.amazon.com/what-is/ai-agents/)  
8. What are AI agents? Definition, examples, and types | Google Cloud, дата последнего обращения: мая 28, 2025, [https://cloud.google.com/discover/what-are-ai-agents](https://cloud.google.com/discover/what-are-ai-agents)  
9. (PDF) Agentic Bioinformatics \- ResearchGate, дата последнего обращения: мая 28, 2025, [https://www.researchgate.net/publication/389284860\_Agentic\_Bioinformatics](https://www.researchgate.net/publication/389284860_Agentic_Bioinformatics)  
10. Model Context Protocol (MCP) an overview \- Philschmid, дата последнего обращения: мая 28, 2025, [https://www.philschmid.de/mcp-introduction](https://www.philschmid.de/mcp-introduction)  
11. Model Context Protocol (MCP): A Guide With Demo Project \- DataCamp, дата последнего обращения: мая 28, 2025, [https://www.datacamp.com/tutorial/mcp-model-context-protocol](https://www.datacamp.com/tutorial/mcp-model-context-protocol)  
12. Introduction to NextFlow \- Bioinformatics Workbook, дата последнего обращения: мая 28, 2025, [https://bioinformaticsworkbook.org/dataAnalysis/nextflow/01\_introductionToNextFlow.html](https://bioinformaticsworkbook.org/dataAnalysis/nextflow/01_introductionToNextFlow.html)  
13. Nextflow | Core Bioinformatics group \- University of Cambridge, дата последнего обращения: мая 28, 2025, [https://www.corebioinf.stemcells.cam.ac.uk/pipelines-tools/pipelines/nextflow](https://www.corebioinf.stemcells.cam.ac.uk/pipelines-tools/pipelines/nextflow)  
14. Introduction to Bioinformatics workflows with Nextflow and nf-core: All in One View, дата последнего обращения: мая 28, 2025, [https://carpentries-incubator.github.io/workflows-nextflow/aio.html](https://carpentries-incubator.github.io/workflows-nextflow/aio.html)  
15. Help with Optimizing Nextflow Pipeline for Large Datasets \- Seqera Community, дата последнего обращения: мая 28, 2025, [https://community.seqera.io/t/help-with-optimizing-nextflow-pipeline-for-large-datasets/1761](https://community.seqera.io/t/help-with-optimizing-nextflow-pipeline-for-large-datasets/1761)  
16. Troubleshooting \- training.nextflow.io, дата последнего обращения: мая 28, 2025, [https://training.nextflow.io/2.1/basic\_training/debugging/](https://training.nextflow.io/2.1/basic_training/debugging/)  
17. Troubleshooting Guide \- Documentation \- EPI2ME, дата последнего обращения: мая 28, 2025, [https://epi2me.nanoporetech.com/epi2me-docs/help/troubleshooting/](https://epi2me.nanoporetech.com/epi2me-docs/help/troubleshooting/)  
18. (PDF) Viash: A meta-framework for building reusable workflow modules \- ResearchGate, дата последнего обращения: мая 28, 2025, [https://www.researchgate.net/publication/377671642\_Viash\_A\_meta-framework\_for\_building\_reusable\_workflow\_modules](https://www.researchgate.net/publication/377671642_Viash_A_meta-framework_for_building_reusable_workflow_modules)  
19. www.theoj.org, дата последнего обращения: мая 28, 2025, [https://www.theoj.org/joss-papers/joss.06089/10.21105.joss.06089.pdf](https://www.theoj.org/joss-papers/joss.06089/10.21105.joss.06089.pdf)  
20. Create a new component \- Viash, дата последнего обращения: мая 28, 2025, [https://viash.io/guide/component/create-component.html](https://viash.io/guide/component/create-component.html)  
21. Viash \- Data Intuitive, дата последнего обращения: мая 28, 2025, [https://www.data-intuitive.com/products/viash.html](https://www.data-intuitive.com/products/viash.html)  
22. Data Intuitive: Where Data Meets Intuition, дата последнего обращения: мая 28, 2025, [https://www.data-intuitive.com/](https://www.data-intuitive.com/)  
23. Intuitive Data Workflow Approach, дата последнего обращения: мая 28, 2025, [https://www.data-intuitive.com/approach/approach.html](https://www.data-intuitive.com/approach/approach.html)  
24. Containerised Bioinformatics, дата последнего обращения: мая 28, 2025, [https://www.melbournebioinformatics.org.au/tutorials/tutorials/docker/media/](https://www.melbournebioinformatics.org.au/tutorials/tutorials/docker/media/)  
25. A Robust Method for Constructing Docker Images for Reproducible Research. \- OSF, дата последнего обращения: мая 28, 2025, [https://osf.io/preprints/osf/8pgd7\_v1](https://osf.io/preprints/osf/8pgd7_v1)  
26. Docker for Bioinformatics Analysis \- Omics tutorials, дата последнего обращения: мая 28, 2025, [https://omicstutorials.com/docker-for-bioinformatics-analysis/](https://omicstutorials.com/docker-for-bioinformatics-analysis/)  
27. Looking for good examples of reproducible scRNA-seq pipeline with Nextflow, Docker, renv, дата последнего обращения: мая 28, 2025, [https://www.reddit.com/r/bioinformatics/comments/1ig3spm/looking\_for\_good\_examples\_of\_reproducible/](https://www.reddit.com/r/bioinformatics/comments/1ig3spm/looking_for_good_examples_of_reproducible/)  
28. Dependencies and containers \- training.nextflow.io, дата последнего обращения: мая 28, 2025, [https://training.nextflow.io/2.1/es/basic\_training/containers/](https://training.nextflow.io/2.1/es/basic_training/containers/)  
29. An open-source protocol for AI agents to interact \- IBM Research, дата последнего обращения: мая 28, 2025, [https://research.ibm.com/blog/agent-communication-protocol-ai](https://research.ibm.com/blog/agent-communication-protocol-ai)  
30. A beginners Guide on Model Context Protocol (MCP) \- OpenCV, дата последнего обращения: мая 28, 2025, [https://opencv.org/blog/model-context-protocol/](https://opencv.org/blog/model-context-protocol/)  
31. Researchers from FutureHouse and ScienceMachine Introduce BixBench: A Benchmark Designed to Evaluate AI Agents on Real-World Bioinformatics Task \- MarkTechPost, дата последнего обращения: мая 28, 2025, [https://www.marktechpost.com/2025/03/04/researchers-from-futurehouse-and-sciencemachine-introduce-bixbench-a-benchmark-designed-to-evaluate-ai-agents-on-real-world-bioinformatics-task/](https://www.marktechpost.com/2025/03/04/researchers-from-futurehouse-and-sciencemachine-introduce-bixbench-a-benchmark-designed-to-evaluate-ai-agents-on-real-world-bioinformatics-task/)  
32. Structured retrieval AI agent tools \- Databricks Documentation, дата последнего обращения: мая 28, 2025, [https://docs.databricks.com/aws/en/generative-ai/agent-framework/structured-retrieval-tools](https://docs.databricks.com/aws/en/generative-ai/agent-framework/structured-retrieval-tools)  
33. With AI Agents on the Scene, Structured Data is Back in Vogue \- RTInsights, дата последнего обращения: мая 28, 2025, [https://www.rtinsights.com/with-ai-agents-on-the-scene-structured-data-is-back-in-vogue/](https://www.rtinsights.com/with-ai-agents-on-the-scene-structured-data-is-back-in-vogue/)  
34. www.akira.ai, дата последнего обращения: мая 28, 2025, [https://www.akira.ai/blog/ai-agents-for-genomic-data-analysis\#:\~:text=Pattern%20Recognition%20Agent%3A%20AI%20Agents,lead%20to%20more%20accurate%20diagnoses.](https://www.akira.ai/blog/ai-agents-for-genomic-data-analysis#:~:text=Pattern%20Recognition%20Agent%3A%20AI%20Agents,lead%20to%20more%20accurate%20diagnoses.)  
35. How AI Agents Enhances Genomic Data Analysis for Precision Healthcare \- Akira AI, дата последнего обращения: мая 28, 2025, [https://www.akira.ai/blog/ai-agents-for-genomic-data-analysis](https://www.akira.ai/blog/ai-agents-for-genomic-data-analysis)  
36. From legacy scripts to ready-to-run Nextflow pipelines with Seqera AI, дата последнего обращения: мая 28, 2025, [https://seqera.io/blog/legacy-scripts-to-nextflow-seqera-ai/](https://seqera.io/blog/legacy-scripts-to-nextflow-seqera-ai/)  
37. Bringing Seqera AI to the Nextflow VS Code extension, дата последнего обращения: мая 28, 2025, [https://seqera.io/blog/seqera-ai--nextflow-vs-code/](https://seqera.io/blog/seqera-ai--nextflow-vs-code/)  
38. aertslab/SpatialNF: Spatial transcriptomics NextFlow pipelines \- GitHub, дата последнего обращения: мая 28, 2025, [https://github.com/aertslab/SpatialNF](https://github.com/aertslab/SpatialNF)  
39. Docs: Troubleshooting basics \- nf-core, дата последнего обращения: мая 28, 2025, [https://nf-co.re/docs/usage/troubleshooting/basics](https://nf-co.re/docs/usage/troubleshooting/basics)
