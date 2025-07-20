# Testing Guide

This document explains how to run the test suite for the OpenProblems MCP Spatial Tools.

## Quick Start

### 1. Verify Test Setup
```bash
python verify_tests.py
```

### 2. Run All Tests
```bash
python run_tests.py
```

### 3. Run Specific Test Categories
```bash
# Quick development tests
python tests/test_runner.py --quick

# Integration tests only
python tests/test_runner.py --integration

# All tests
python tests/test_runner.py --all
```

## Alternative Test Methods

### Using pytest directly
```bash
# Install in development mode first
pip install -e .

# Run all tests
pytest tests/ -v

# Run specific test file
pytest tests/test_spatial_validation.py -v

# Run with coverage
pytest tests/ --cov=openproblems_mcp --cov-report=html
```

### Using the standalone test script
```bash
python test_spatial_tools.py
```

## Test Structure

- **`tests/test_spatial_validation.py`** - Tests for SpatialDataValidator
- **`tests/test_metadata_analysis.py`** - Tests for BioinformaticsMetadataExtractor
- **`tests/test_spatial_tools.py`** - Tests for SpatialMCPTools (MCP interface)
- **`tests/conftest.py`** - Pytest fixtures and configuration
- **`tests/test_runner.py`** - Custom test runner with different modes

## Troubleshooting

### Import Errors
If you get `ModuleNotFoundError: No module named 'openproblems_mcp'`:

1. Try installing in development mode:
   ```bash
   pip install -e .
   ```

2. Or use the provided test runners that handle paths automatically:
   ```bash
   python run_tests.py
   python verify_tests.py
   ```

### Missing Dependencies
The tests are designed to work without heavy bioinformatics dependencies (spatialdata, zarr, anndata). They use mocking to simulate these libraries.

Required for testing:
```bash
pip install pytest pytest-mock
```

### Path Issues
The test runners automatically set up the Python path. If you're having issues:

1. Run from the project root directory
2. Use the provided test runners instead of calling pytest directly
3. Check that `src/openproblems_mcp/__init__.py` exists

## Test Coverage

The test suite covers:

- ✅ All validation methods and levels
- ✅ All metadata extraction functionality
- ✅ All MCP tool interfaces
- ✅ Error handling and edge cases
- ✅ File format detection
- ✅ Result formatting
- ✅ Integration workflows

## Development Workflow

1. **Before making changes**: Run `python verify_tests.py`
2. **During development**: Run `python tests/test_runner.py --quick`
3. **Before committing**: Run `python run_tests.py`
4. **For CI/CD**: Use `pytest tests/ -v --tb=short`

## Adding New Tests

When adding new functionality:

1. Add unit tests to the appropriate test file
2. Add integration tests for complete workflows
3. Update fixtures in `conftest.py` if needed
4. Run the full test suite to ensure nothing breaks

## Performance

The test suite is optimized for speed:
- Uses mocking to avoid expensive operations
- Creates minimal test data
- Focuses on logic rather than data processing
- Typical run time: < 30 seconds for full suite
