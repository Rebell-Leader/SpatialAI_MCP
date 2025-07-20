"""Unit tests for spatial transcriptomics data validation module."""

import pytest
import tempfile
import json
from pathlib import Path
from unittest.mock import Mock, patch, mock_open

from openproblems_mcp.spatial_validation import (
    SpatialDataValidator,
    ValidationLevel,
    ValidationStatus,
    ValidationIssue,
    ValidationResult
)


class TestSpatialDataValidator:
    """Test cases for SpatialDataValidator."""

    def setup_method(self):
        """Set up test fixtures."""
        self.validator = SpatialDataValidator()

    def test_init(self):
        """Test validator initialization."""
        assert isinstance(self.validator, SpatialDataValidator)
        assert hasattr(self.validator, 'dependencies')
        assert isinstance(self.validator.dependencies, dict)

    def test_detect_format_spatialdata(self):
        """Test format detection for SpatialData files."""
        # Test .spatialdata extension
        path = Path("test.spatialdata")
        result = self.validator._detect_format(path)
        assert result == "spatialdata", f"Expected 'spatialdata', got '{result}' for {path}"

        # Test .zarr extension (should check if it's SpatialData)
        path = Path("test.zarr")
        with patch.object(self.validator, '_is_spatialdata_zarr', return_value=True):
            result = self.validator._detect_format(path)
            assert result == "spatialdata", f"Expected 'spatialdata', got '{result}' for zarr with SpatialData"

        with patch.object(self.validator, '_is_spatialdata_zarr', return_value=False):
            result = self.validator._detect_format(path)
            assert result == "zarr", f"Expected 'zarr', got '{result}' for zarr without SpatialData"

    def test_detect_format_anndata(self):
        """Test format detection for AnnData files."""
        for ext in ['.h5ad', '.h5', '.hdf5']:
            path = Path(f"test{ext}")
            result = self.validator._detect_format(path)
            assert result == "anndata"

    def test_detect_format_unknown(self):
        """Test format detection for unknown files."""
        path = Path("test.txt")
        result = self.validator._detect_format(path)
        assert result == "unknown"

    def test_is_spatialdata_zarr_with_metadata(self):
        """Test SpatialData zarr detection with metadata."""
        with tempfile.TemporaryDirectory() as temp_dir:
            zarr_path = Path(temp_dir) / "test.zarr"
            zarr_path.mkdir()

            # Create .zattrs file with spatialdata metadata
            zattrs_file = zarr_path / ".zattrs"
            with open(zattrs_file, 'w') as f:
                json.dump({"spatialdata": "0.1.0"}, f)

            result = self.validator._is_spatialdata_zarr(zarr_path)
            assert result is True

    def test_is_spatialdata_zarr_with_structure(self):
        """Test SpatialData zarr detection with directory structure."""
        with tempfile.TemporaryDirectory() as temp_dir:
            zarr_path = Path(temp_dir) / "test.zarr"
            zarr_path.mkdir()

            # Create SpatialData-like structure
            (zarr_path / "images").mkdir()
            (zarr_path / "table").mkdir()

            result = self.validator._is_spatialdata_zarr(zarr_path)
            assert result is True

    def test_validate_file_nonexistent(self):
        """Test validation of non-existent file."""
        result = self.validator.validate_file("nonexistent.zarr")

        assert result.is_valid is False
        assert len(result.issues) == 1
        assert result.issues[0].level == ValidationStatus.CRITICAL
        assert "does not exist" in result.issues[0].message

    def test_validate_file_unknown_format(self):
        """Test validation of unknown format file."""
        with tempfile.NamedTemporaryFile(suffix=".txt") as temp_file:
            result = self.validator.validate_file(temp_file.name)

            assert result.is_valid is False
            assert result.format_type == "unknown"
            assert any("Unknown or unsupported format" in issue.message for issue in result.issues)

    @patch('pathlib.Path.exists', return_value=True)
    @patch('pathlib.Path.is_dir', return_value=True)
    @patch('pathlib.Path.stat')
    def test_validate_spatialdata_basic(self, mock_stat, mock_is_dir, mock_exists):
        """Test basic SpatialData validation."""
        mock_stat.return_value.st_size = 1000

        with patch.object(self.validator, '_detect_format', return_value='spatialdata'):
            with patch('pathlib.Path.iterdir', return_value=[]):
                result = self.validator.validate_file("test.spatialdata")

                assert result.format_type == "spatialdata"
                assert result.metadata['format_details'] == 'SpatialData'

    @patch('pathlib.Path.exists', return_value=True)
    @patch('pathlib.Path.is_file', return_value=True)
    @patch('pathlib.Path.stat')
    def test_validate_anndata_basic(self, mock_stat, mock_is_file, mock_exists):
        """Test basic AnnData validation."""
        mock_stat.return_value.st_size = 1000

        with patch.object(self.validator, '_detect_format', return_value='anndata'):
            result = self.validator.validate_file("test.h5ad")

            assert result.format_type == "anndata"
            assert result.metadata['format_details'] == 'AnnData'

    def test_validate_multiple_files(self):
        """Test validation of multiple files."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create test files
            file1 = Path(temp_dir) / "test1.h5ad"
            file2 = Path(temp_dir) / "test2.zarr"
            file1.touch()
            file2.mkdir()

            results = self.validator.validate_multiple_files([str(file1), str(file2)])

            assert len(results) == 2
            assert all(isinstance(r, ValidationResult) for r in results)

    def test_get_validation_summary(self):
        """Test validation summary generation."""
        # Create mock results
        result1 = ValidationResult(
            file_path="test1.h5ad",
            format_type="anndata",
            is_valid=True,
            validation_level=ValidationLevel.BASIC,
            issues=[],
            metadata={},
            file_size=1000
        )

        result2 = ValidationResult(
            file_path="test2.zarr",
            format_type="spatialdata",
            is_valid=False,
            validation_level=ValidationLevel.BASIC,
            issues=[ValidationIssue(ValidationStatus.ERROR, "test", "Test error")],
            metadata={},
            file_size=2000
        )

        summary = self.validator.get_validation_summary([result1, result2])

        assert summary['total_files'] == 2
        assert summary['valid_files'] == 1
        assert summary['files_with_errors'] == 1
        assert summary['success_rate'] == 0.5
        assert 'anndata' in summary['format_distribution']
        assert 'spatialdata' in summary['format_distribution']


class TestValidationIssue:
    """Test cases for ValidationIssue dataclass."""

    def test_validation_issue_creation(self):
        """Test ValidationIssue creation."""
        issue = ValidationIssue(
            level=ValidationStatus.WARNING,
            category="test",
            message="Test message",
            details="Test details",
            suggested_fix="Test fix"
        )

        assert issue.level == ValidationStatus.WARNING
        assert issue.category == "test"
        assert issue.message == "Test message"
        assert issue.details == "Test details"
        assert issue.suggested_fix == "Test fix"


class TestValidationResult:
    """Test cases for ValidationResult dataclass."""

    def test_validation_result_properties(self):
        """Test ValidationResult properties."""
        # Test with no issues
        result = ValidationResult(
            file_path="test.h5ad",
            format_type="anndata",
            is_valid=True,
            validation_level=ValidationLevel.BASIC,
            issues=[],
            metadata={},
            file_size=1000
        )

        assert result.has_errors is False
        assert result.has_warnings is False

        # Test with warnings
        result.issues = [ValidationIssue(ValidationStatus.WARNING, "test", "Warning")]
        assert result.has_warnings is True
        assert result.has_errors is False

        # Test with errors
        result.issues.append(ValidationIssue(ValidationStatus.ERROR, "test", "Error"))
        assert result.has_errors is True
        assert result.has_warnings is True


class TestValidationLevels:
    """Test cases for different validation levels."""

    def setup_method(self):
        """Set up test fixtures."""
        self.validator = SpatialDataValidator()

    def test_validation_levels_enum(self):
        """Test ValidationLevel enum values."""
        assert ValidationLevel.BASIC.value == "basic"
        assert ValidationLevel.STRUCTURE.value == "structure"
        assert ValidationLevel.INTEGRITY.value == "integrity"
        assert ValidationLevel.DOMAIN.value == "domain"

    @patch('pathlib.Path.exists', return_value=True)
    @patch('pathlib.Path.is_dir', return_value=True)
    @patch('pathlib.Path.stat')
    def test_different_validation_levels(self, mock_stat, mock_is_dir, mock_exists):
        """Test validation with different levels."""
        mock_stat.return_value.st_size = 1000

        with patch.object(self.validator, '_detect_format', return_value='spatialdata'):
            with patch('pathlib.Path.iterdir', return_value=[]):
                # Test each validation level
                for level in ValidationLevel:
                    result = self.validator.validate_file("test.spatialdata", level)
                    assert result.validation_level == level


class TestSpatialDataValidatorIntegration:
    """Integration tests for SpatialDataValidator."""

    def setup_method(self):
        """Set up test fixtures."""
        self.validator = SpatialDataValidator()

    def test_end_to_end_validation_workflow(self):
        """Test complete validation workflow."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create a mock SpatialData structure
            spatialdata_path = Path(temp_dir) / "test.spatialdata"
            spatialdata_path.mkdir()

            # Add basic structure
            (spatialdata_path / "images").mkdir()
            (spatialdata_path / "table").mkdir()

            # Add metadata file
            zattrs_file = spatialdata_path / ".zattrs"
            with open(zattrs_file, 'w') as f:
                json.dump({"spatialdata": "0.1.0"}, f)

            # Validate with different levels
            for level in ValidationLevel:
                result = self.validator.validate_file(str(spatialdata_path), level)
                assert result.format_type == "spatialdata"
                assert result.validation_level == level

    def test_batch_validation_workflow(self):
        """Test batch validation workflow."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create multiple test files
            files = []
            for i, ext in enumerate(['.h5ad', '.zarr', '.spatialdata']):
                file_path = Path(temp_dir) / f"test{i}{ext}"
                if ext == '.h5ad':
                    file_path.touch()
                else:
                    file_path.mkdir()
                files.append(str(file_path))

            # Validate all files
            results = self.validator.validate_multiple_files(files)
            assert len(results) == 3

            # Generate summary
            summary = self.validator.get_validation_summary(results)
            assert summary['total_files'] == 3
            assert len(summary['format_distribution']) > 0


if __name__ == "__main__":
    pytest.main([__file__])
