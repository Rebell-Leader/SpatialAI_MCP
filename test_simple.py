#!/usr/bin/env python3
"""Simple test to verify pytest configuration works."""

def test_pytest_setup():
    """Test that pytest can run with the current configuration."""
    try:
        from openproblems_mcp.spatial_validation import SpatialDataValidator
        from openproblems_mcp.metadata_analysis import BioinformaticsMetadataExtractor
        from openproblems_mcp.spatial_tools import SpatialMCPTools

        # Basic initialization test
        validator = SpatialDataValidator()
        extractor = BioinformaticsMetadataExtractor()
        tools = SpatialMCPTools()

        assert validator is not None
        assert extractor is not None
        assert tools is not None

        print("✅ Simple pytest test passed!")

    except ImportError as e:
        print(f"❌ Import failed: {e}")
        raise
    except Exception as e:
        print(f"❌ Test failed: {e}")
        raise

if __name__ == "__main__":
    test_pytest_setup()
    print("✅ Direct execution test passed!")
