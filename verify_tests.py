#!/usr/bin/env python3
"""Verify that the test setup is working correctly."""

import sys
from pathlib import Path

# Add src to Python path
project_root = Path(__file__).parent
src_path = project_root / "src"
sys.path.insert(0, str(src_path))

def verify_imports():
    """Verify that all modules can be imported."""
    print("🔍 Verifying module imports...")

    try:
        # Test core imports
        from openproblems_mcp.spatial_validation import (
            SpatialDataValidator, ValidationLevel, ValidationStatus, ValidationIssue, ValidationResult
        )
        print("✅ spatial_validation imports work")

        from openproblems_mcp.metadata_analysis import (
            BioinformaticsMetadataExtractor, MetadataType, MetadataField, MetadataExtractionResult
        )
        print("✅ metadata_analysis imports work")

        from openproblems_mcp.spatial_tools import SpatialMCPTools
        print("✅ spatial_tools imports work")

        return True

    except ImportError as e:
        print(f"❌ Import failed: {e}")
        return False

def verify_basic_functionality():
    """Verify basic functionality works."""
    print("\n🔧 Verifying basic functionality...")

    try:
        from openproblems_mcp.spatial_validation import SpatialDataValidator, ValidationLevel
        from openproblems_mcp.metadata_analysis import BioinformaticsMetadataExtractor
        from openproblems_mcp.spatial_tools import SpatialMCPTools

        # Test initialization
        validator = SpatialDataValidator()
        extractor = BioinformaticsMetadataExtractor()
        tools = SpatialMCPTools()

        print("✅ All classes initialize correctly")

        # Test enum values
        levels = [level.value for level in ValidationLevel]
        print(f"✅ ValidationLevel enum works: {levels}")

        # Test a simple method call
        result = tools.validate_spatial_data("nonexistent.h5ad")
        print(f"✅ Method calls work (result length: {len(result)})")

        return True

    except Exception as e:
        print(f"❌ Functionality test failed: {e}")
        return False

def verify_test_structure():
    """Verify test file structure."""
    print("\n📁 Verifying test file structure...")

    test_files = [
        "tests/test_spatial_validation.py",
        "tests/test_metadata_analysis.py",
        "tests/test_spatial_tools.py",
        "tests/conftest.py",
        "tests/test_runner.py"
    ]

    missing_files = []
    for test_file in test_files:
        if not Path(test_file).exists():
            missing_files.append(test_file)

    if missing_files:
        print(f"❌ Missing test files: {missing_files}")
        return False
    else:
        print("✅ All test files exist")
        return True

def main():
    """Run verification checks."""
    print("🧪 Test Setup Verification")
    print("=" * 30)

    checks = [
        ("Module Imports", verify_imports),
        ("Basic Functionality", verify_basic_functionality),
        ("Test File Structure", verify_test_structure)
    ]

    passed = 0
    total = len(checks)

    for check_name, check_func in checks:
        if check_func():
            passed += 1
        print()

    print("=" * 30)
    print(f"Verification Results: {passed}/{total} passed")

    if passed == total:
        print("🎉 Test setup is working correctly!")
        print("\nYou can now run tests with:")
        print("  python run_tests.py")
        print("  python tests/test_runner.py --all")
        print("  pytest tests/ -v")
        return 0
    else:
        print("❌ Test setup has issues that need to be fixed")
        return 1

if __name__ == "__main__":
    sys.exit(main())
