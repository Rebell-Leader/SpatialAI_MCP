"""Unit tests for spatial MCP tools module."""

import pytest
import json
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

from src.openproblems_mcp.spatial_tools import SpatialMCPTools
from src.openproblems_mcp.spatial_validation import ValidationResult, ValidationLevel, ValidationStatus, ValidationIssue
from src.openproblems_mcp.metadata_analysis import MetadataExtractionResult, MetadataField


class TestSpatialMCPTools:
    """Test cases for SpatialMCPTools."""

    def setup_method(self):
        """Set up test fixtures."""
        self.tools = SpatialMCPTools()

    def test_init(self):
        """Test SpatialMCPTools initialization."""
        assert hasattr(self.tools, 'validator')
        assert hasattr(self.tools, 'metadata_extractor')

    @patch('src.openproblems_mcp.spatial_tools.SpatialDataValidator')
    def test_validate_spatial_data_success(self, mock_validator_class):
        """Test successful spatial data validation."""
        # Mock validator and result
        mock_validator = Mock()
        mock_validator_class.return_value = mock_validator

        mock_result = ValidationResult(
            file_path="test.h5ad",
            format_type="anndata",
            is_valid=True,
            validation_level=ValidationLevel.STRUCTURE,
            issues=[],
            metadata={"n_obs": 1000, "n_vars": 2000},
            file_size=1024000
        )
        mock_validator.validate_file.return_value = mock_result

        # Create new tools instance to use mocked validator
        tools = SpatialMCPTools()

        result = tools.validate_spatial_data("test.h5ad", "structure", "summary")

        assert "✅" in result
        assert "test.h5ad" in result
        assert "anndata" in result
        assert "Valid" in result

    @patch('src.openproblems_mcp.spatial_tools.SpatialDataValidator')
    def test_validate_spatial_data_with_issues(self, mock_validator_class):
        """Test spatial data validation with issues."""
        # Mock validator and result with issues
        mock_validator = Mock()
        mock_validator_class.return_value = mock_validator

        mock_result = ValidationResult(
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
        mock_validator.validate_file.return_value = mock_result

        # Create new tools instance to use mocked validator
        tools = SpatialMCPTools()

        result = tools.validate_spatial_data("test.h5ad", "structure", "summary")

        assert "❌" in result
        assert "Invalid" in result
        assert "Issues Found: 2" in result
        assert "Errors: 1" in result
        assert "Warnings: 1" in result

    def test_validate_spatial_data_invalid_level(self):
        """Test validation with invalid level."""
        result = self.tools.validate_spatial_data("test.h5ad", "invalid_level")

        assert "❌ Invalid validation level" in result
        assert "invalid_level" in result

    def test_validate_multiple_spatial_files_invalid_json(self):
        """Test multiple file validation with invalid JSON."""
        result = self.tools.validate_multiple_spatial_files("invalid json")

        assert "❌ Invalid JSON format" in result

    def test_validate_multiple_spatial_files_not_list(self):
        """Test multiple file validation with non-list JSON."""
        result = self.tools.validate_multiple_spatial_files('{"not": "a list"}')

        assert "❌ file_paths must be a JSON list" in result

    @patch('src.openproblems_mcp.spatial_tools.SpatialDataValidator')
    def test_validate_multiple_spatial_files_success(self, mock_validator_class):
        """Test successful multiple file validation."""
        # Mock validator
        mock_validator = Mock()
        mock_validator_class.return_value = mock_validator

        # Mock results
        mock_results = [
            ValidationResult("test1.h5ad", "anndata", True, ValidationLevel.STRUCTURE, [], {}, 1000),
            ValidationResult("test2.zarr", "spatialdata", False, ValidationLevel.STRUCTURE,
                           [ValidationIssue(ValidationStatus.ERROR, "test", "Test error")], {}, 2000)
        ]
        mock_validator.validate_multiple_files.return_value = mock_results

        # Mock summary
        mock_summary = {
            'total_files': 2,
            'valid_files': 1,
            'files_with_warnings': 0,
            'files_with_errors': 1,
            'success_rate': 0.5,
            'format_distribution': {'anndata': 1, 'spatialdata': 1}
        }
        mock_validator.get_validation_summary.return_value = mock_summary

        # Create new tools instance
        tools = SpatialMCPTools()

        result = tools.validate_multiple_spatial_files('["test1.h5ad", "test2.zarr"]')

        assert "📊 Multiple File Validation Summary" in result
        assert "Total Files: 2" in result
        assert "Valid Files: 1" in result
        assert "Success Rate: 50.0%" in result

    @patch('src.openproblems_mcp.spatial_tools.SpatialDataValidator')
    def test_analyze_spatial_metadata_success(self, mock_validator_class):
        """Test successful spatial metadata analysis."""
        # Mock validator
        mock_validator = Mock()
        mock_validator_class.return_value = mock_validator

        mock_result = ValidationResult(
            file_path="test.h5ad",
            format_type="anndata",
            is_valid=True,
            validation_level=ValidationLevel.DOMAIN,
            issues=[],
            metadata={
                "n_obs": 1000,
                "n_vars": 2000,
                "spatial_keys": ["X_spatial"],
                "components": ["images", "table"]
            },
            file_size=1024000
        )
        mock_validator.validate_file.return_value = mock_result

        # Create new tools instance
        tools = SpatialMCPTools()

        result = tools.analyze_spatial_metadata("test.h5ad")

        assert "📊 Spatial Data Metadata Analysis" in result
        assert "anndata" in result
        assert "Number of Observations: 1000" in result
        assert "Number of Variables: 2000" in result
        assert "Spatial Coordinate Keys: X_spatial" in result

    @patch('src.openproblems_mcp.spatial_tools.SpatialDataValidator')
    def test_analyze_spatial_metadata_validation_failed(self, mock_validator_class):
        """Test spatial metadata analysis with validation failure."""
        # Mock validator
        mock_validator = Mock()
        mock_validator_class.return_value = mock_validator

        mock_result = ValidationResult(
            file_path="test.h5ad",
            format_type="anndata",
            is_valid=False,
            validation_level=ValidationLevel.DOMAIN,
            issues=[ValidationIssue(ValidationStatus.ERROR, "critical", "File corrupted")],
            metadata={},
            file_size=1024
        )
        mock_validator.validate_file.return_value = mock_result

        # Create new tools instance
        tools = SpatialMCPTools()

        result = tools.analyze_spatial_metadata("test.h5ad")

        assert "❌ Cannot analyze metadata - file validation failed" in result
        assert "File corrupted" in result

    def test_check_spatial_data_compatibility_invalid_json(self):
        """Test compatibility check with invalid JSON."""
        result = self.tools.check_spatial_data_compatibility("invalid json")

        assert "❌ Invalid JSON format" in result

    def test_check_spatial_data_compatibility_insufficient_files(self):
        """Test compatibility check with insufficient files."""
        result = self.tools.check_spatial_data_compatibility('["single_file.h5ad"]')

        assert "❌ Need at least 2 file paths" in result

    @patch('src.openproblems_mcp.spatial_tools.SpatialDataValidator')
    def test_check_spatial_data_compatibility_success(self, mock_validator_class):
        """Test successful compatibility check."""
        # Mock validator
        mock_validator = Mock()
        mock_validator_class.return_value = mock_validator

        # Mock validation results
        mock_results = [
            ValidationResult("test1.h5ad", "anndata", True, ValidationLevel.DOMAIN, [],
                           {"n_obs": 1000, "n_vars": 2000, "spatial_keys": ["X_spatial"]}, 1000),
            ValidationResult("test2.h5ad", "anndata", True, ValidationLevel.DOMAIN, [],
                           {"n_obs": 1200, "n_vars": 2000, "spatial_keys": ["X_spatial"]}, 1200)
        ]
        mock_validator.validate_multiple_files.return_value = mock_results

        # Create new tools instance
        tools = SpatialMCPTools()

        result = tools.check_spatial_data_compatibility('["test1.h5ad", "test2.h5ad"]')

        assert "🔗 Spatial Data Compatibility Analysis" in result
        assert "Analyzing 2 files" in result
        assert "✅ All files use the same format: anndata" in result
        assert "✅ Consistent gene count: 2000 genes" in result

    @patch('src.openproblems_mcp.spatial_tools.BioinformaticsMetadataExtractor')
    def test_extract_bioinformatics_metadata_success(self, mock_extractor_class):
        """Test successful bioinformatics metadata extraction."""
        # Mock extractor
        mock_extractor = Mock()
        mock_extractor_class.return_value = mock_extractor

        mock_result = MetadataExtractionResult(
            file_path="test.nf",
            file_type="nextflow_workflow",
            extraction_success=True,
            metadata_fields=[
                MetadataField("component_name", "test_workflow", "str", "workflow", "Workflow name", "critical"),
                MetadataField("process_count", 3, "int", "workflow", "Number of processes", "high")
            ],
            quality_metrics={"completeness": 0.8},
            issues=[],
            suggestions=["Add more documentation"]
        )
        mock_extractor.extract_metadata.return_value = mock_result

        # Create new tools instance
        tools = SpatialMCPTools()

        result = tools.extract_bioinformatics_metadata("test.nf")

        assert "🔍 Bioinformatics Metadata Extraction" in result
        assert "nextflow_workflow" in result
        assert "✅ Success" in result
        assert "component_name" in result
        assert "process_count" in result

    def test_format_file_size(self):
        """Test file size formatting."""
        # Test various file sizes
        assert self.tools._format_file_size(0) == "0 B"
        assert self.tools._format_file_size(512) == "512.0 B"
        assert self.tools._format_file_size(1024) == "1.0 KB"
        assert self.tools._format_file_size(1024 * 1024) == "1.0 MB"
        assert self.tools._format_file_size(1024 * 1024 * 1024) == "1.0 GB"

    def test_result_to_dict(self):
        """Test ValidationResult to dictionary conversion."""
        result = ValidationResult(
            file_path="test.h5ad",
            format_type="anndata",
            is_valid=True,
            validation_level=ValidationLevel.STRUCTURE,
            issues=[ValidationIssue(ValidationStatus.WARNING, "test", "Test warning")],
            metadata={"test": "value"},
            file_size=1024
        )

        result_dict = self.tools._result_to_dict(result)

        assert result_dict['file_path'] == "test.h5ad"
        assert result_dict['format_type'] == "anndata"
        assert result_dict['is_valid'] is True
        assert result_dict['validation_level'] == "structure"
        assert len(result_dict['issues']) == 1
        assert result_dict['issues'][0]['level'] == "warning"
        assert result_dict['metadata'] == {"test": "value"}
        assert result_dict['file_size'] == 1024


class TestSpatialMCPToolsIntegration:
    """Integration tests for SpatialMCPTools."""

    def setup_method(self):
        """Set up test fixtures."""
        self.tools = SpatialMCPTools()

    def test_end_to_end_validation_workflow(self):
        """Test complete validation workflow."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create a mock data file
            test_file = Path(temp_dir) / "test.h5ad"
            test_file.touch()

            # Test validation (will fail due to invalid format, but tests the workflow)
            result = self.tools.validate_spatial_data(str(test_file))

            # Should contain validation attempt
            assert isinstance(result, str)
            assert len(result) > 0

    def test_json_parameter_handling(self):
        """Test JSON parameter handling in various methods."""
        # Test valid JSON
        valid_json = '["file1.h5ad", "file2.zarr"]'

        # These should not crash with JSON parsing errors
        result1 = self.tools.validate_multiple_spatial_files(valid_json)
        result2 = self.tools.check_spatial_data_compatibility(valid_json)
        result3 = self.tools.assess_data_quality(valid_json)

        # All should be strings (even if they contain error messages)
        assert isinstance(result1, str)
        assert isinstance(result2, str)
        assert isinstance(result3, str)

    def test_error_handling(self):
        """Test error handling in various scenarios."""
        # Test with non-existent file
        result = self.tools.validate_spatial_data("nonexistent_file.h5ad")
        assert "❌" in result or "Validation failed" in result

        # Test with invalid JSON
        result = self.tools.validate_multiple_spatial_files("invalid json")
        assert "❌" in result

        # Test with empty file list
        result = self.tools.check_spatial_data_compatibility("[]")
        assert "❌" in result


class TestSpatialMCPToolsFormatting:
    """Test cases for formatting methods in SpatialMCPTools."""

    def setup_method(self):
        """Set up test fixtures."""
        self.tools = SpatialMCPTools()

    def test_format_result_summary(self):
        """Test summary formatting."""
        result = ValidationResult(
            file_path="test.h5ad",
            format_type="anndata",
            is_valid=True,
            validation_level=ValidationLevel.STRUCTURE,
            issues=[],
            metadata={},
            file_size=1024
        )

        summary = self.tools._format_result_summary(result)

        assert "✅" in summary
        assert "test.h5ad" in summary
        assert "anndata" in summary
        assert "Valid" in summary
        assert "1.0 KB" in summary

    def test_format_result_detailed(self):
        """Test detailed formatting."""
        result = ValidationResult(
            file_path="test.h5ad",
            format_type="anndata",
            is_valid=False,
            validation_level=ValidationLevel.STRUCTURE,
            issues=[
                ValidationIssue(ValidationStatus.ERROR, "structure", "Test error", "Error details", "Fix suggestion")
            ],
            metadata={"test_key": "test_value"},
            file_size=2048
        )

        detailed = self.tools._format_result_detailed(result)

        assert "❌" in detailed
        assert "Detailed Spatial Data Validation Report" in detailed
        assert "test.h5ad" in detailed
        assert "anndata" in detailed
        assert "Test error" in detailed
        assert "Error details" in detailed
        assert "Fix suggestion" in detailed
        assert "test_key" in detailed

    def test_format_result_json(self):
        """Test JSON formatting."""
        result = ValidationResult(
            file_path="test.h5ad",
            format_type="anndata",
            is_valid=True,
            validation_level=ValidationLevel.STRUCTURE,
            issues=[],
            metadata={},
            file_size=1024
        )

        json_result = self.tools._format_result_json(result)

        # Should be valid JSON
        parsed = json.loads(json_result)
        assert parsed['file_path'] == "test.h5ad"
        assert parsed['format_type'] == "anndata"
        assert parsed['is_valid'] is True

    def test_format_issues(self):
        """Test issue formatting."""
        issues = [
            ValidationIssue(ValidationStatus.ERROR, "test", "Error message", suggested_fix="Fix it"),
            ValidationIssue(ValidationStatus.WARNING, "test", "Warning message")
        ]

        formatted = self.tools._format_issues(issues)

        assert "❌ Error message" in formatted
        assert "💡 Fix it" in formatted
        assert "⚠️ Warning message" in formatted


if __name__ == "__main__":
    pytest.main([__file__])
