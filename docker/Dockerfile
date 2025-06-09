# Multi-stage build for optimized Docker image
FROM python:3.11-slim as python-base

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    curl \
    wget \
    ca-certificates \
    openjdk-17-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# Install Docker CLI (for building images)
RUN curl -fsSL https://download.docker.com/linux/debian/gpg | gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg \
    && echo "deb [arch=amd64 signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/debian bullseye stable" > /etc/apt/sources.list.d/docker.list \
    && apt-get update && apt-get install -y --no-install-recommends docker-ce-cli \
    && rm -rf /var/lib/apt/lists/*

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/ \
    && chmod +x /usr/local/bin/nextflow

# Install Viash
RUN curl -fsSL get.viash.io | bash -s -- --bin /usr/local/bin

# Create non-root user
RUN useradd --create-home --shell /bin/bash openproblems

# Set working directory
WORKDIR /app

# Copy requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the application
COPY src/ ./src/
COPY pyproject.toml ./

# Install the package
RUN pip install -e .

# Create necessary directories
RUN mkdir -p /app/logs /app/data /app/work \
    && chown -R openproblems:openproblems /app

# Switch to non-root user
USER openproblems

# Set environment variables for the user
ENV PATH="/home/openproblems/.local/bin:$PATH"

# Expose the default MCP port (not required for stdio but useful for HTTP transport)
EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=30s --start-period=5s --retries=3 \
    CMD python -c "import mcp; print('MCP SDK available')" || exit 1

# Default command
CMD ["python", "-m", "mcp_server.main"]
