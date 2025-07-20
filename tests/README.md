# Spatial Transcriptomics Tools Test Suite

This directory contains comprehensive unit tests for the spatial transcriptomics data validation and metadata analysis tools implemented in the OpenProblems MCP server.

## Test Structure

### Core Test Files

1. **`test_spatial_validation.py`** - Tests for the `SpatialDataValidator` class
   - Format detection (SpatialData, zarr, AnnData)
   - Validation levels (basic, structure, integrity, domain)
   - File validation workflows
   - Multiple file validation
   - Validation result summarization

2. **`test_metadata_analysis.py`** - Tests for the `BioinformaticsMetadataExtractor` class
   - File type detection (Nextflow, Viash, spatial data)
   - Metadata extraction from various formats
   - Workflow dependency analysis
   - Quality assessment

3. **`test_spatial_tools.py`** - Tests for the `SpatialMCPTools` class
   - MCP tool integration
   - Result formatting (summary, detailed, JSON)
   - Error handling
   - Parameter validation

### Supporting Files

- **`conftest.py`** - Pytest configuration and fixtures
- **`test_runner.py`** - Test runner with different modes
- **`README.md`** - This documentation

## Test Categories

### Unit Tests
- Individual component testing
- Mock-based testing for external dependencies
- Input validation and error handling

### Integration Tests
- End-to-end workflow testing
- File system interaction
- Cross-component integration

### Fixtures and Test Data

The test suite includes fixtures for creating:
- Mock SpatialData directory structures
- Sample Nextflow workflows and configurations
- Sample Viash component configurations
- Temporary test files and directories

## Running Tests

### Prerequisites
```bash
pip install pytest pytest-mock
```

### Run All Tests
```bash
python tests/test_runner.py --all
```

### Run Quick Tests (Development)
```bash
python tests/test_runner.py --quick
```

### Run Integration Tests Only
```bash
python tests/test_runner.py --integration
```

### Run Specific Test Files
```bash
pytest tests/test_spatial_validation.py -v
pytest tests/test_metadata_analysis.py -v
pytest tests/test_spatial_tools.py -v
```

### Run Tests with Coverage
```bash
pytest --cov=src/openproblems_mcp tests/ --cov-report=html
```

## Test Coverage

The test suite covers:

### SpatialDataValidator
- ✅ Format detection for all supported formats
- ✅ Validation levels and their behaviors
- ✅ File existence and accessibility checks
- ✅ Domain-specific validation rules
- ✅ Error handling and reporting
- ✅ Multiple file validation workflows
- ✅ Validation result summarization

### BioinformaticsMetadataExtractor
- ✅ File type detection patterns
- ✅ Nextflow config/workflow parsing
- ✅ Viash component configuration parsing
- ✅ Python/R/Bash script analysis
- ✅ Dependency extraction and analysis
- ✅ Quality metrics calculation
- ✅ Suggestion generation

### SpatialMCPTools
- ✅ All MCP tool methods
- ✅ Parameter validation and JSON handling
- ✅ Result formatting in multiple formats
- ✅ Error handling and graceful degradation
- ✅ File size formatting utilities
- ✅ Integration with validation and metadata modules

## Mock Strategy

The tests use extensive mocking to:
- Avoid dependencies on external libraries (spatialdata, zarr, anndata, h5py)
- Test error conditions and edge cases
- Ensure consistent test execution across environments
- Speed up test execution

## Test Data

Test fixtures create realistic but minimal test data:
- Mock SpatialData directory structures with proper metadata
- Sample Nextflow workflows with processes and configurations
- Viash component configurations with platforms and arguments
- Python scripts with bioinformatics library imports

## Continuous Integration

The test suite is designed to run in CI environments without requiring:
- Heavy bioinformatics dependencies
- Large test data files
- External services or databases

## Adding New Tests

When adding new functionality:

1. Add unit tests for individual methods
2. Add integration tests for workflows
3. Update fixtures if new test data patterns are needed
4. Ensure mocking covers external dependencies
5. Add documentation for new test categories

## Test Debugging

For debugging test failures:

```bash
# Run with detailed output
pytest tests/ -v --tb=long

# Run specific test with debugging
pytest tests/test_spatial_validation.py::TestSpatialDataValidator::test_validate_file_nonexistent -v -s

# Run with pdb on failure
pytest tests/ --pdb
```

## Performance Considerations

The test suite is optimized for speed:
- Uses temporary directories that are automatically cleaned up
- Mocks expensive operations (file I/O, library imports)
- Focuses on logic testing rather than data processing
- Parallel test execution supported

## Known Limitations

- Tests don't validate actual spatial transcriptomics data processing
- Some integration tests require file system access
- Mock objects may not perfectly replicate library behavior
- Performance tests are not included (focused on correctness)

## Future Enhancements

Potential test suite improvements:
- Property-based testing with hypothesis
- Performance benchmarking tests
- Real data integration tests (optional)
- Cross-platform compatibility tests
- Memory usage validation
