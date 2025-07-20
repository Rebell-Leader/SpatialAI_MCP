"""Test runner for spatial transcriptomics validation and metadata analysis tools."""

import sys
import os
import pytest
from pathlib import Path

# Set up paths
project_root = Path(__file__).parent.parent
src_path = project_root / "src"
sys.path.insert(0, str(src_path))

# Change to project root for consistent test execution
os.chdir(project_root)


def run_spatial_tests():
    """Run all spatial transcriptomics related tests."""
    test_files = [
        "tests/test_spatial_validation.py",
        "tests/test_metadata_analysis.py",
        "tests/test_spatial_tools.py"
    ]

    print("🧪 Running Spatial Transcriptomics Tool Tests")
    print("=" * 50)

    # Run tests with verbose output
    exit_code = pytest.main([
        "-v",
        "--tb=short",
        "--color=yes",
        f"--pythonpath={src_path}",
        *test_files
    ])

    if exit_code == 0:
        print("\n✅ All spatial transcriptomics tests passed!")
    else:
        print(f"\n❌ Some tests failed (exit code: {exit_code})")

    return exit_code


def run_quick_tests():
    """Run a quick subset of tests for development."""
    print("🚀 Running Quick Test Suite")
    print("=" * 30)

    # Run basic import and initialization tests
    exit_code = pytest.main([
        "-v",
        "--tb=short",
        "--color=yes",
        f"--pythonpath={src_path}",
        "-k", "test_init or test_detect_format or test_metadata_field_creation",
        "tests/test_spatial_validation.py",
        "tests/test_metadata_analysis.py",
        "tests/test_spatial_tools.py"
    ])

    return exit_code


def run_integration_tests():
    """Run integration tests."""
    print("🔗 Running Integration Tests")
    print("=" * 30)

    exit_code = pytest.main([
        "-v",
        "--tb=short",
        "--color=yes",
        f"--pythonpath={src_path}",
        "-k", "integration or end_to_end",
        "tests/test_spatial_validation.py",
        "tests/test_metadata_analysis.py",
        "tests/test_spatial_tools.py"
    ])

    return exit_code


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run spatial transcriptomics tool tests")
    parser.add_argument("--quick", action="store_true", help="Run quick test suite")
    parser.add_argument("--integration", action="store_true", help="Run integration tests only")
    parser.add_argument("--all", action="store_true", help="Run all tests (default)")

    args = parser.parse_args()

    if args.quick:
        exit_code = run_quick_tests()
    elif args.integration:
        exit_code = run_integration_tests()
    else:
        exit_code = run_spatial_tests()

    sys.exit(exit_code)
