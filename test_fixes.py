#!/usr/bin/env python3
"""Test script to verify the fixes work correctly."""

import sys
from pathlib import Path

# Add src to Python path
project_root = Path(__file__).parent
src_path = project_root / "src"
sys.path.insert(0, str(src_path))

def test_file_type_detection():
    """Test that file type detection works correctly."""
    print("🔍 Testing file type detection...")

    try:
        from openproblems_mcp.metadata_analysis import BioinformaticsMetadataExtractor
        from openproblems_mcp.spatial_validation import SpatialDataValidator

        # Test metadata extractor
        extractor = BioinformaticsMetadataExtractor()

        # Test script.py detection
        script_path = Path("script.py")
        result = extractor._detect_file_type(script_path)
        print(f"   script.py detected as: {result}")
        assert result == "viash_script", f"Expected 'viash_script', got '{result}'"

        # Test nextflow.config detection
        config_path = Path("nextflow.config")
        result = extractor._detect_file_type(config_path)
        print(f"   nextflow.config detected as: {result}")
        assert result == "nextflow_config", f"Expected 'nextflow_config', got '{result}'"

        # Test spatial validator
        validator = SpatialDataValidator()

        # Test .spatialdata detection
        spatialdata_path = Path("test.spatialdata")
        result = validator._detect_format(spatialdata_path)
        print(f"   test.spatialdata detected as: {result}")
        assert result == "spatialdata", f"Expected 'spatialdata', got '{result}'"

        print("✅ File type detection works correctly")
        return True

    except Exception as e:
        print(f"❌ File type detection test failed: {e}")
        return False

def test_validation_formatting():
    """Test that validation result formatting works correctly."""
    print("\n📝 Testing validation result formatting...")

    try:
        from openproblems_mcp.spatial_tools import SpatialMCPTools
        from openproblems_mcp.spatial_validation import ValidationResult, ValidationLevel, ValidationStatus, ValidationIssue

        tools = SpatialMCPTools()

        # Create a mock result with issues
        result = ValidationResult(
            file_path="test.h5ad",
            format_type="anndata",
            is_valid=False,
            validation_level=ValidationLevel.STRUCTURE,
            issues=[
                ValidationIssue(ValidationStatus.ERROR, "structure", "Missing required component"),
                ValidationIssue(ValidationStatus.WARNING, "domain", "Low gene count")
            ],
            metadata={},
            file_size=1024
        )

        # Test formatting
        formatted = tools._format_result_summary(result)
        print(f"   Formatted result contains '**Issues Found:** 2': {'**Issues Found:** 2' in formatted}")
        assert "**Issues Found:** 2" in formatted, "Expected '**Issues Found:** 2' in formatted result"

        print("✅ Validation result formatting works correctly")
        return True

    except Exception as e:
        print(f"❌ Validation formatting test failed: {e}")
        return False

def main():
    """Run fix verification tests."""
    print("🔧 Testing Bug Fixes")
    print("=" * 25)

    tests = [
        ("File Type Detection", test_file_type_detection),
        ("Validation Formatting", test_validation_formatting)
    ]

    passed = 0
    total = len(tests)

    for test_name, test_func in tests:
        if test_func():
            passed += 1

    print("\n" + "=" * 25)
    print(f"Fix Tests: {passed}/{total} passed")

    if passed == total:
        print("🎉 All fixes verified!")
        return 0
    else:
        print("❌ Some fixes need attention!")
        return 1

if __name__ == "__main__":
    sys.exit(main())
