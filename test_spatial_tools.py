#!/usr/bin/env python3
"""Quick test script to verify spatial tools are working correctly."""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

def test_imports():
    """Test that all modules can be imported."""
    print("Testing imports...")

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
    """Test basic functionality of the tools."""
    print("\nTesting basic functionality...")

    try:
        from openproblems_mcp.spatial_validation import SpatialDataValidator, ValidationLevel
        from openproblems_mcp.metadata_analysis import BioinformaticsMetadataExtractor
        from openproblems_mcp.spatial_tools import SpatialMCPTools

        # Test SpatialDataValidator
        validator = SpatialDataValidator()
        print(f"✅ SpatialDataValidator initialized")
        print(f"   Dependencies: {validator.dependencies}")

        # Test validation levels
levels = [level.value for level in ValidationLevel]
        print(f"   Validation levels: {levels}")

        # Test BioinformaticsMetadataExtractor
        extractor = BioinformaticsMetadataExtractor()
        print(f"✅ BioinformaticsMetadataExtractor initialized")
        print(f"   Dependencies: {extractor.dependencies}")

        # Test SpatialMCPTools
        tools = SpatialMCPTools()
        print(f"✅ SpatialMCPTools initialized")

        return True
    except Exception as e:
        print(f"❌ Basic functionality test failed: {e}")
        return False

def test_validation_workflow():
    """Test validation workflow with non-existent file."""
    print("\nTesting validation workflow...")

    try:
        from openproblems_mcp.spatial_tools import SpatialMCPTools

        tools = SpatialMCPTools()

        # Test with non-existent file (should handle gracefully)
        result = tools.validate_spatial_data("nonexistent_file.h5ad")
        print("✅ Validation workflow completed")
        print(f"   Result type: {type(result)}")
        print(f"   Result length: {len(result)} characters")

        # Test invalid validation level
        result = tools.validate_spatial_data("test.h5ad", "invalid_level")
        assert "Invalid validation level" in result
        print("✅ Invalid validation level handled correctly")

        return True
    except Exception as e:
        print(f"❌ Validation workflow test failed: {e}")
        return False

def test_metadata_extraction():
    """Test metadata extraction workflow."""
    print("\nTesting metadata extraction...")

    try:
        from openproblems_mcp.spatial_tools import SpatialMCPTools

        tools = SpatialMCPTools()

        # Test with non-existent file (should handle gracefully)
        result = tools.extract_bioinformatics_metadata("nonexistent_file.nf")
        print("✅ Metadata extraction workflow completed")
        print(f"   Result type: {type(result)}")
        print(f"   Result length: {len(result)} characters")

        return True
    except Exception as e:
        print(f"❌ Metadata extraction test failed: {e}")
        return False

def main():
    """Run all tests."""
    print("🧪 Testing Spatial Transcriptomics Tools")
    print("=" * 40)

    tests = [
        test_imports,
        test_basic_functionality,
        test_validation_workflow,
        test_metadata_extraction
    ]

    passed = 0
    total = len(tests)

    for test in tests:
        if test():
            passed += 1
        print()

    print("=" * 40)
    print(f"Tests passed: {passed}/{total}")

    if passed == total:
        print("🎉 All tests passed!")
        return 0
    else:
        print("❌ Some tests failed!")
        return 1

if __name__ == "__main__":
    sys.exit(main())
