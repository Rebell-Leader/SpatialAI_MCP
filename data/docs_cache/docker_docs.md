# Docker Best Practices for Bioinformatics

## Multi-stage Builds

### Optimized Python Environment
```dockerfile
# Build stage
FROM python:3.9-slim as builder
WORKDIR /build
COPY requirements.txt .
RUN pip install --no-cache-dir --user -r requirements.txt

# Production stage
FROM python:3.9-slim
COPY --from=builder /root/.local /root/.local
RUN apt-get update && apt-get install -y procps
WORKDIR /app
```

### Bioinformatics Stack
```dockerfile
FROM python:3.9-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    libblas-dev \
    liblapack-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    scanpy>=1.9.0 \
    anndata>=0.8.0 \
    pandas>=1.5.0 \
    numpy>=1.21.0

WORKDIR /app
```

### OpenProblems Compatible Container
```dockerfile
FROM python:3.9-slim

RUN apt-get update && apt-get install -y procps
RUN pip install --no-cache-dir scanpy anndata pandas numpy

# Create non-root user for Nextflow
RUN groupadd -g 1000 nextflow && \
    useradd -u 1000 -g nextflow nextflow

USER nextflow
WORKDIR /app
ENTRYPOINT ["python"]
```

## Best Practices
- Use specific versions for reproducibility
- Use minimal base images
- Create non-root users
- Combine RUN commands to reduce layers
- Use health checks for services
- Set appropriate resource limits
