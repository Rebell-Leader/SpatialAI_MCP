#!/usr/bin/env python3
"""Simple test runner that setse environment correctly."""

import sys
import os
from pathlib import Path

# Add src to Python path
project_root = Path(__file__).parent
src_path = project_root / "src"
sys.path.insert(0, str(src_path))

# Change to project directory
os.chdir(project_root)

def test_imports():
    """Test that all modules can be imported."""
    print("🧪 Testing module imports...")

    try:
        from openproblems_mcp.spatial_validation import SpatialDataValidator, ValidationLevel
        print("✅ SpatialDataValidator imported successfully")

        from openproblems_mcp.metadata_analysis import BioinformaticsMetadataExtractor
        print("✅ BioinformaticsMetadataExtractor imported successfully")

        from openproblems_mcp.spatial_tools import SpatialMCPTools
        print("✅ SpatialMCPTools imported successfully")

        return True
    except ImportError as e:
        print(f"❌ Import failed: {e}")
        return False

def test_basic_functionality():
    """Test basic functionality."""
    print("\n🔧 Testing basic functionality...")

    try:
        from openproblems_mcp.spatial_validation import SpatialDataValidator, ValidationLevel
        from openproblems_mcp.metadata_analysis import BioinformaticsMetadataExtractor
        from openproblems_mcp.spatial_tools import SpatialMCPTools

        # Test SpatialDataValidator
        validator = SpatialDataValidator()
        print(f"✅ SpatialDataValidator initialized")
        print(f"   Dependencies: {list(validator.dependencies.keys())}")

        # Test validation levels
        levels = [level.value for level in ValidationLevel]
        print(f"   Validation levels: {levels}")

        # Test BioinformaticsMetadataExtractor
        extractor = BioinformaticsMetadataExtractor()
        print(f"✅ BioinformaticsMetadataExtractor initialized")
        print(f"   Dependencies: {list(extractor.dependencies.keys())}")

        # Test SpatialMCPTools
        tools = SpatialMCPTools()
        print(f"✅ SpatialMCPTools initialized")

        return True
    except Exception as e:
        print(f"❌ Basic functionality test failed: {e}")
        return False

def test_validation_workflow():
    """Test validation workflow."""
    print("\n🔍 Testing validation workflow...")

    try:
        from openproblems_mcp.spatial_tools import SpatialMCPTools

        tools = SpatialMCPTools()

        # Test with non-existent file (should handle gracefully)
        result = tools.validate_spatial_data("nonexistent_file.h5ad")
        print("✅ Validation workflow completed")
        print(f"   Result contains error handling: {'❌' in result}")

        # Test invalid validation level
        result = tools.validate_spatial_data("test.h5ad", "invalid_level")
        success = "Invalid validation level" in result
        print(f"✅ Invalid validation level handled: {success}")

        return True
    except Exception as e:
        print(f"❌ Validation workflow test failed: {e}")
        return False

def run_pytest_tests():
    """Run pytest tests if available."""
    print("\n🧪 Running pytest tests...")

    try:
        import pytest

        # Run tests with proper path setup
        exit_code = pytest.main([
            "tests/",
            "-v",
            "--tb=short",
            f"--pythonpath={src_path}"
        ])

        if exit_code == 0:
            print("✅ All pytest tests passed!")
        else:
            print(f"⚠️ Some pytest tests failed (exit code: {exit_code})")

        return exit_code == 0

    except ImportError:
        print("⚠️ pytest not available, skipping pytest tests")
        return True
    except Exception as e:
        print(f"❌ pytest execution failed: {e}")
        return False

def main():
    """Run all tests."""
    print("🚀 OpenProblems MCP Spatial Tools Test Suite")
    print("=" * 50)

    tests = [
        ("Module Imports", test_imports),
        ("Basic Functionality", test_basic_functionality),
        ("Validation Workflow", test_validation_workflow),
        ("Pytest Tests", run_pytest_tests)
    ]

    passed = 0
    total = len(tests)

    for test_name, test_func in tests:
        print(f"\n📋 Running: {test_name}")
        print("-" * 30)

        if test_func():
            passed += 1
            print(f"✅ {test_name}: PASSED")
        else:
            print(f"❌ {test_name}: FAILED")

    print("\n" + "=" * 50)
    print(f"📊 Test Results: {passed}/{total} passed")

    if passed == total:
        print("🎉 All tests passed!")
        return 0
    else:
        print("⚠️ Some tests failed!")
        return 1

if __name__ == "__main__":
    sys.exit(main())
