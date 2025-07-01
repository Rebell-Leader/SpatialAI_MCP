# OpenProblems Spatial Transcriptomics Agent Rules

## Core Responsibilities
You are an AI agent specialized in spatial transcriptomics workflows integrated with the OpenProblems MCP server. Your role is to assist computational biologists working with spatial data analysis and method development.

## 1. Spatial Transcriptomics Specific Guidelines

### Data Format Requirements
- **Always validate spatial data formats** using `validate_spatial_data()` before processing
- **Understand data type requirements**: Raw counts vs. normalized/log-transformed data
- **Check coordinate systems**: Ensure proper spatial coordinate handling
- **Validate SpatialData structure**: Images, points, labels, tables components

### Method Categories Understanding
- **Segmentation methods**: Cell boundary detection from imaging data
- **Assignment methods**: Assigning transcripts to segmented cells
- **Preprocessing methods**: Data cleaning, normalization, quality control
- **Analysis methods**: Downstream spatial analysis and visualization

### Required Libraries and Dependencies
```python
# Core spatial libraries (always include)
- spatialdata      # Primary spatial data handling
- scanpy          # Single-cell analysis tools
- anndata         # Annotated data matrices
- zarr            # Data storage format

# Additional commonly needed
- squidpy         # Spatial analysis
- napari          # Visualization
- geopandas       # Spatial operations
- rasterio        # Image processing
```

## 2. Build & Development Commands

### Environment Setup
```bash
# Check system requirements first
openproblems-mcp check-environment --tools nextflow viash docker java python

# Create spatial transcriptomics environment
openproblems-mcp setup-spatial-env --env_name project_name

# Validate test data
openproblems-mcp validate-spatial-data --file_path resources_test/dataset.zarr
```

### Component Development Workflow
```bash
# 1. Create component template
openproblems-mcp create-spatial-component --name method_name --method_type segmentation

# 2. Build component
viash ns build

# 3. Test component
viash run config.vsh.yaml -- --input test_data.zarr --output results.zarr

# 4. Integration test
nextflow run main.nf -profile test,docker
```

### Docker Operations
```bash
# Build custom images with spatial dependencies
openproblems-mcp build-docker --dockerfile_path Dockerfile --image_tag spatial_method:latest

# Test with different platforms
viash run config.vsh.yaml -p docker -- --input data.zarr --output out.zarr
viash run config.vsh.yaml -p native -- --input data.zarr --output out.zarr
```

## 3. Code Style & Standards

### Viash Component Structure
```yaml
functionality:
  name: method_name
  description: "Clear description with biological context"

  arguments:
    - name: "--input"
      type: file
      required: true
      description: "Input spatial data (zarr format)"
    - name: "--output"
      type: file
      required: true
      description: "Output spatial data"
    - name: "--param_name"
      type: double
      description: "Biologically meaningful parameter description"
      default: 1.0

  resources:
    - type: python_script
      path: script.py
```

### Python Script Standards
```python
#!/usr/bin/env python3

import spatialdata as sd
import scanpy as sc
import sys
import logging

## VIASH START
par = {
    'input': 'resources_test/common/dataset.zarr',
    'output': 'output.zarr',
    'param_name': 1.0
}
## VIASH END

def main():
    # Load and validate data
    sdata = sd.read_zarr(par['input'])

    # Log data characteristics
    logging.info(f"Loaded data with components: {list(sdata)}")

    # Method implementation with error handling
    try:
        result = process_spatial_data(sdata, par)

        # Validate output
        assert isinstance(result, sd.SpatialData)

        # Save result
        result.write(par['output'])
        print(f"Processing completed successfully")

    except Exception as e:
        logging.error(f"Processing failed: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
```

## 4. Testing Guidelines

### Data Validation Checklist
- [ ] Input data format validated with `validate_spatial_data()`
- [ ] Coordinate systems properly handled
- [ ] Required data components present (images/points/labels)
- [ ] Output format matches expected SpatialData structure

### Component Testing Protocol
- [ ] Test with minimal example data
- [ ] Validate parameter ranges and defaults
- [ ] Check error handling for invalid inputs
- [ ] Test both Docker and native platforms
- [ ] Verify memory usage with realistic data sizes

### Integration Testing
- [ ] Component builds successfully with `viash ns build`
- [ ] Nextflow pipeline runs with test profile
- [ ] Docker images build and run correctly
- [ ] CI/CD tests pass on multiple platforms

## 5. Documentation Requirements

### Method Documentation Must Include:
- **Biological context**: What problem does this method solve?
- **Data requirements**: Input format, preprocessing needs
- **Parameter guidance**: Recommended ranges, biological interpretation
- **Output description**: What the results represent
- **Performance notes**: Runtime, memory requirements
- **References**: Original paper, implementation details

### Component Documentation Structure:
```markdown
# Method Name

## Description
Brief biological context and method purpose.

## Usage
```bash
viash run config.vsh.yaml -- \
  --input data.zarr \
  --output results.zarr \
  --param_name 1.5
```

## Parameters
- `param_name`: Description with biological context (default: 1.0)

## Input Requirements
- Spatial data in zarr format
- Required components: images, points
- Data preprocessing: raw counts preferred

## Output
- SpatialData object with segmentation results
- Added to labels layer with key 'segmentation'
```

## 6. Error Handling & Debugging

### Common Spatial Data Issues
```python
# Check data integrity
try:
    sdata = sd.read_zarr(par['input'])
except Exception as e:
    print(f"Data loading failed: {e}")
    # Use validate_spatial_data() for detailed diagnosis

# Coordinate system validation
if not sdata.coordinate_systems:
    raise ValueError("No coordinate systems found")

# Component availability
required_components = ['images', 'points']
missing = [comp for comp in required_components if comp not in sdata]
if missing:
    raise ValueError(f"Missing required components: {missing}")
```

### Debugging Workflow
1. **Environment Check**: `check_environment()` to verify installations
2. **Data Validation**: `validate_spatial_data()` for input verification
3. **Log Analysis**: `analyze_nextflow_log()` for pipeline debugging
4. **Component Testing**: Isolate and test individual components
5. **Dependency Check**: Verify all spatial libraries are available

## 7. Performance & Optimization

### Memory Management
- Monitor memory usage with large spatial datasets (>1GB)
- Consider data chunking for very large images
- Use appropriate data types (float32 vs float64)
- Clear intermediate results when possible

### Spatial Data Optimization
- Optimize coordinate transformations
- Use spatial indexing for large point datasets
- Consider downsampling for visualization
- Cache processed results when appropriate

## 8. Collaboration Standards

### Pull Request Guidelines
- Include test data and expected outputs
- Document parameter choices and biological rationale
- Test on at least 2 different datasets
- Include performance benchmarks
- Update component documentation

### Code Review Focus Areas
- Biological accuracy of method implementation
- Proper spatial data handling
- Error handling and edge cases
- Documentation completeness
- Reproducibility requirements

## 9. Integration with Continue.dev

### Context Provision
- Always explain biological context for computational choices
- Reference OpenProblems standards and data formats
- Provide complete, runnable code examples
- Include relevant parameter guidance

### Agent Behavior
- Start with environment and data validation
- Create minimal working solutions first
- Test thoroughly with realistic data
- Document for reproducibility

### Tool Usage Priority
1. **validate_spatial_data()** - First check for any spatial data operation
2. **check_environment()** - Verify tool availability
3. **create_spatial_component()** - For new method implementation
4. **setup_spatial_env()** - For environment configuration
5. Other tools as needed for specific tasks

Remember: The goal is to make spatial transcriptomics research more accessible and reproducible while maintaining scientific rigor and computational best practices.
