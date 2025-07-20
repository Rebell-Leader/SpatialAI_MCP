"""Spatial transcriptomics data validation module.

This module provides validation for spatial transcriptomics data formats including
SpatialData, zarr, and AnnData formats with domain-specific integrity checking.
"""

import logging
import os
import json
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Any, Union, Tuple
from dataclasses import dataclass
from enum import Enum

logger = logging.getLogger(__name__)


class ValidationLevel(Enum):
    """Validation levels for spatial data."""
    BASIC = "basic"          # File existence and basic structure
    STRUCTURE = "structure"  # Format-specific structure validation
    INTEGRITY = "integrity"  # Data integrity and consistency
    DOMAIN = "domain"        # Spatial biology domain-specific validation


class ValidationStatus(Enum):
    """Validation result status."""
    VALID = "valid"
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"


@dataclass
class ValidationIssue:
    """Represents a validation issue found during data validation."""
    level: ValidationStatus
    category: str
    message: str
    details: Optional[str] = None
    suggested_fix: Optional[str] = None
    file_path: Optional[str] = None
    line_number: Optional[int] = None


@dataclass
class ValidationResult:
    """Result of spatial data validation."""
    file_path: str
    format_type: str
    is_valid: bool
    validation_level: ValidationLevel
    issues: List[ValidationIssue]
    metadata: Dict[str, Any]
    file_size: int

    @property
    def has_errors(self) -> bool:
        """Check if validation found any errors."""
        return any(issue.level in [ValidationStatus.ERROR, ValidationStatus.CRITICAL]
                  for issue in self.issues)

    @property
    def has_warnings(self) -> bool:
        """Check if validation found any warnings."""
        return any(issue.level == ValidationStatus.WARNING for issue in self.issues)


class SpatialDataValidator:
    """Validator for spatial transcriptomics data formats."""

    # File extensions for different formats
    SPATIALDATA_EXTENSIONS = {'.zarr', '.spatialdata'}
    ZARR_EXTENSIONS = {'.zarr'}
    ANNDATA_EXTENSIONS = {'.h5ad', '.h5', '.hdf5'}

    # Expected spatial data components
    SPATIALDATA_COMPONENTS = {
        'images', 'labels', 'points', 'shapes', 'table'
    }

    # Common AnnData attributes
    ANNDATA_ATTRIBUTES = {
        'X', 'obs', 'var', 'uns', 'obsm', 'varm', 'obsp', 'varp'
    }

    def __init__(self):
        """Initialize the spatial data validator."""
        self._check_optional_dependencies()

    def _check_optional_dependencies(self) -> Dict[str, bool]:
        """Check availability of optional spatial transcriptomics libraries."""
        dependencies = {}

        try:
            import spatialdata
            dependencies['spatialdata'] = True
        except ImportError:
            dependencies['spatialdata'] = False

        try:
            import zarr
            dependencies['zarr'] = True
        except ImportError:
            dependencies['zarr'] = False

        try:
            import anndata
            dependencies['anndata'] = True
        except ImportError:
            dependencies['anndata'] = False

        try:
            import h5py
            dependencies['h5py'] = True
        except ImportError:
            dependencies['h5py'] = False

        self.dependencies = dependencies
        return dependencies

    def validate_file(
        self,
        file_path: Union[str, Path],
        validation_level: ValidationLevel = ValidationLevel.STRUCTURE
    ) -> ValidationResult:
        """Validate a spatial transcriptomics data file.

        Args:
            file_path: Path to the file to validate
            validation_level: Level of validation to perform

        Returns:
            ValidationResult with validation details
        """
        file_path = Path(file_path)

        # Initialize result
        result = ValidationResult(
            file_path=str(file_path),
            format_type="unknown",
            is_valid=False,
            validation_level=validation_level,
            issues=[],
            metadata={},
            file_size=0
        )

        # Basic file existence check
        if not file_path.exists():
            result.issues.append(ValidationIssue(
                level=ValidationStatus.CRITICAL,
                category="file_access",
                message=f"File does not exist: {file_path}",
                suggested_fix="Check file path and ensure file exists"
            ))
            return result

        # Get file size
        try:
            if file_path.is_file():
                result.file_size = file_path.stat().st_size
            elif file_path.is_dir():
                result.file_size = sum(f.stat().st_size for f in file_path.rglob('*') if f.is_file())
        except Exception as e:
            result.issues.append(ValidationIssue(
                level=ValidationStatus.WARNING,
                category="file_access",
                message=f"Could not determine file size: {e}"
            ))

        # Determine format type
        format_type = self._detect_format(file_path)
        result.format_type = format_type

        # Perform format-specific validation
        if format_type == "spatialdata":
            self._validate_spatialdata(file_path, result, validation_level)
        elif format_type == "zarr":
            self._validate_zarr(file_path, result, validation_level)
        elif format_type == "anndata":
            self._validate_anndata(file_path, result, validation_level)
        else:
            result.issues.append(ValidationIssue(
                level=ValidationStatus.ERROR,
                category="format_detection",
                message=f"Unknown or unsupported format for file: {file_path}",
                details=f"File extension: {file_path.suffix}",
                suggested_fix="Ensure file has a supported extension (.zarr, .spatialdata, .h5ad, .h5, .hdf5)"
            ))
            return result

        # Set overall validity
        result.is_valid = not result.has_errors

        return result

    def _detect_format(self, file_path: Path) -> str:
        """Detect the format of a spatial data file."""
        suffix = file_path.suffix.lower()

        # Check for SpatialData
        if suffix in self.SPATIALDATA_EXTENSIONS or file_path.name.endswith('.spatialdata'):
            return "spatialdata"

        # Check for zarr (but not SpatialData)
        if suffix in self.ZARR_EXTENSIONS:
            # Check if it's actually a SpatialData object
            if self._is_spatialdata_zarr(file_path):
                return "spatialdata"
            return "zarr"

        # Check for AnnData
        if suffix in self.ANNDATA_EXTENSIONS:
            return "anndata"

        return "unknown"

    def _is_spatialdata_zarr(self, file_path: Path) -> bool:
        """Check if a zarr file is actually a SpatialData object."""
        try:
            # Look for SpatialData-specific metadata
            metadata_file = file_path / '.zattrs'
            if metadata_file.exists():
                with open(metadata_file, 'r') as f:
                    metadata = json.load(f)
                    return 'spatialdata' in str(metadata).lower()

            # Check for typical SpatialData structure
            expected_dirs = {'images', 'labels', 'points', 'shapes', 'table'}
            actual_dirs = {d.name for d in file_path.iterdir() if d.is_dir()}
            return len(expected_dirs.intersection(actual_dirs)) >= 2

        except Exception:
            return False

    def _validate_spatialdata(
        self,
        file_path: Path,
        result: ValidationResult,
        validation_level: ValidationLevel
    ) -> None:
        """Validate SpatialData format."""
        result.metadata['format_details'] = 'SpatialData'

        # Basic structure validation
        if not file_path.is_dir():
            result.issues.append(ValidationIssue(
                level=ValidationStatus.ERROR,
                category="structure",
                message="SpatialData should be a directory (zarr format)",
                suggested_fix="Ensure the SpatialData object is saved as a zarr directory"
            ))
            return

        # Check for zarr metadata
        zattrs_file = file_path / '.zattrs'
        if not zattrs_file.exists():
            result.issues.append(ValidationIssue(
                level=ValidationStatus.ERROR,
                category="structure",
                message="Missing .zattrs file in SpatialData directory",
                suggested_fix="Ensure the SpatialData object was saved properly"
            ))

        # Check for expected components
        components = {d.name for d in file_path.iterdir() if d.is_dir()}
        result.metadata['components'] = list(components)

        expected_components = self.SPATIALDATA_COMPONENTS
        missing_components = expected_components - components

        if len(components) == 0:
            result.issues.append(ValidationIssue(
                level=ValidationStatus.CRITICAL,
                category="structure",
                message="No SpatialData components found",
                suggested_fix="Ensure the SpatialData object contains at least one component (images, labels, points, shapes, or table)"
            ))
        elif len(missing_components) == len(expected_components):
            result.issues.append(ValidationIssue(
                level=ValidationStatus.WARNING,
                category="structure",
                message="No standard SpatialData components found",
                details=f"Found components: {components}",
                suggested_fix="Verify this is a valid SpatialData object"
            ))

        # Domain-specific validation
        if validation_level in [ValidationLevel.INTEGRITY, ValidationLevel.DOMAIN]:
            self._validate_spatialdata_domain(file_path, result, components)

    def _validate_spatialdata_domain(
        self,
        file_path: Path,
        result: ValidationResult,
        components: set
    ) -> None:
        """Perform domain-specific validation for SpatialData."""
        # Check for spatial transcriptomics specific patterns
        if 'table' in components:
            table_path = file_path / 'table'
            if table_path.exists():
                # Check if table looks like gene expression data
                self._check_expression_table_structure(table_path, result)
                self._validate_table_spatial_annotations(table_path, result)

        if 'images' in components:
            # Validate image data for spatial context
            images_path = file_path / 'images'
            if images_path.exists():
                self._check_spatial_images(images_path, result)
                self._validate_image_spatial_properties(images_path, result)

        if 'points' in components:
            # Validate points data
            points_path = file_path / 'points'
            if points_path.exists():
                self._validate_spatial_points(points_path, result)

        if 'shapes' in components:
            # Validate shapes data
            shapes_path = file_path / 'shapes'
            if shapes_path.exists():
                self._validate_spatial_shapes(shapes_path, result)

        # Check for coordinate systems
        if not any(comp in components for comp in ['images', 'points', 'shapes']):
            result.issues.append(ValidationIssue(
                level=ValidationStatus.WARNING,
                category="domain",
                message="No spatial coordinate data found",
                details="SpatialData object lacks images, points, or shapes components",
                suggested_fix="Ensure spatial coordinates are included for spatial transcriptomics analysis"
            ))

        # Check component relationships
        self._validate_component_relationships(file_path, result, components)

    def _validate_zarr(
        self,
        file_path: Path,
        result: ValidationResult,
        validation_level: ValidationLevel
    ) -> None:
        """Validate zarr format."""
        result.metadata['format_details'] = 'Zarr'

        if not file_path.is_dir():
            result.issues.append(ValidationIssue(
                level=ValidationStatus.ERROR,
                category="structure",
                message="Zarr should be a directory",
                suggested_fix="Ensure the zarr object is saved as a directory"
            ))
            return

        # Check for zarr metadata files
        zarray_files = list(file_path.rglob('.zarray'))
        zattrs_files = list(file_path.rglob('.zattrs'))

        if not zarray_files:
            result.issues.append(ValidationIssue(
                level=ValidationStatus.ERROR,
                category="structure",
                message="No .zarray files found in zarr directory",
                suggested_fix="Ensure this is a valid zarr array"
            ))

        result.metadata['zarray_count'] = len(zarray_files)
        result.metadata['zattrs_count'] = len(zattrs_files)

        # Validate zarr structure if libraries available
        if self.dependencies.get('zarr', False) and validation_level != ValidationLevel.BASIC:
            self._validate_zarr_with_library(file_path, result)

    def _validate_zarr_with_library(self, file_path: Path, result: ValidationResult) -> None:
        """Validate zarr using the zarr library if available."""
        try:
            import zarr

            # Try to open the zarr array
            z = zarr.open(str(file_path), mode='r')

            result.metadata['zarr_shape'] = z.shape if hasattr(z, 'shape') else None
            result.metadata['zarr_dtype'] = str(z.dtype) if hasattr(z, 'dtype') else None
            result.metadata['zarr_chunks'] = z.chunks if hasattr(z, 'chunks') else None

            # Check for reasonable array properties
            if hasattr(z, 'shape') and len(z.shape) == 0:
                result.issues.append(ValidationIssue(
                    level=ValidationStatus.WARNING,
                    category="structure",
                    message="Zarr array has no dimensions",
                    suggested_fix="Verify this is the expected array structure"
                ))

            # Domain-specific validation for spatial transcriptomics zarr arrays
            if hasattr(z, 'shape') and z.shape:
                # Check if this looks like expression data (2D matrix)
                if len(z.shape) == 2:
                    n_obs, n_vars = z.shape
                    result.metadata['expression_matrix_shape'] = [n_obs, n_vars]

                    # Check for reasonable dimensions
                    if n_vars < 10:
                        result.issues.append(ValidationIssue(
                            level=ValidationStatus.WARNING,
                            category="domain",
                            message=f"Very few features/genes: {n_vars}",
                            suggested_fix="Verify this contains gene expression data"
                        ))
                    elif n_vars > 100000:
                        result.issues.append(ValidationIssue(
                            level=ValidationStatus.WARNING,
                            category="domain",
                            message=f"Very many features: {n_vars}",
                            details="This may impact memory usage",
                            suggested_fix="Consider feature filtering if appropriate"
                        ))

                    if n_obs < 10:
                        result.issues.append(ValidationIssue(
                            level=ValidationStatus.WARNING,
                            category="domain",
                            message=f"Very few observations/cells: {n_obs}",
                            suggested_fix="Verify sufficient data for analysis"
                        ))

                # Check data type appropriateness
                if hasattr(z, 'dtype'):
                    dtype_str = str(z.dtype)
                    if 'float' in dtype_str and len(z.shape) == 2:
                        # Check if data looks normalized (values between 0-1 or log-scale)
                        try:
                            sample_data = z[:min(100, z.shape[0]), :min(100, z.shape[1])]
                            data_min, data_max = float(sample_data.min()), float(sample_data.max())

                            if data_max <= 1.0 and data_min >= 0.0:
                                result.metadata['data_appears_normalized'] = True
                                result.issues.append(ValidationIssue(
                                    level=ValidationStatus.WARNING,
                                    category="domain",
                                    message="Data appears to be normalized (0-1 range)",
                                    details="Raw counts may be preferred for some analyses",
                                    suggested_fix="Verify if raw or normalized data is expected"
                                ))
                            elif data_min < 0:
                                result.metadata['data_has_negative_values'] = True
                                result.issues.append(ValidationIssue(
                                    level=ValidationStatus.WARNING,
                                    category="domain",
                                    message="Expression data contains negative values",
                                    details="This suggests log-transformed or centered data",
                                    suggested_fix="Verify data preprocessing is appropriate"
                                ))
                        except Exception:
                            pass  # Skip data sampling if it fails

            # Check chunking strategy for large arrays
            if hasattr(z, 'chunks') and hasattr(z, 'shape') and z.shape:
                total_elements = 1
                for dim in z.shape:
                    total_elements *= dim

                if total_elements > 1e8:  # Large array
                    chunk_size = 1
                    for chunk_dim in z.chunks:
                        chunk_size *= chunk_dim

                    if chunk_size > 1e6:  # Large chunks
                        result.issues.append(ValidationIssue(
                            level=ValidationStatus.WARNING,
                            category="performance",
                            message="Large zarr chunks may impact performance",
                            details=f"Chunk size: {chunk_size:,} elements",
                            suggested_fix="Consider smaller chunk sizes for better I/O performance"
                        ))

        except Exception as e:
            result.issues.append(ValidationIssue(
                level=ValidationStatus.ERROR,
                category="integrity",
                message=f"Failed to read zarr array: {e}",
                suggested_fix="Check zarr array integrity and format"
            ))

    def _validate_anndata(
        self,
        file_path: Path,
        result: ValidationResult,
        validation_level: ValidationLevel
    ) -> None:
        """Validate AnnData format."""
        result.metadata['format_details'] = 'AnnData'

        if not file_path.is_file():
            result.issues.append(ValidationIssue(
                level=ValidationStatus.ERROR,
                category="structure",
                message="AnnData should be a file",
                suggested_fix="Ensure the AnnData object is saved as an .h5ad file"
            ))
            return

        # Basic HDF5 structure check
        if self.dependencies.get('h5py', False):
            self._validate_anndata_hdf5_structure(file_path, result)

        # Validate with AnnData library if available
        if self.dependencies.get('anndata', False) and validation_level != ValidationLevel.BASIC:
            self._validate_anndata_with_library(file_path, result, validation_level)

    def _validate_anndata_hdf5_structure(self, file_path: Path, result: ValidationResult) -> None:
        """Validate AnnData HDF5 structure."""
        try:
            import h5py

            with h5py.File(file_path, 'r') as f:
                # Check for expected AnnData structure
                expected_keys = {'X', 'obs', 'var'}
                found_keys = set(f.keys())
                result.metadata['hdf5_keys'] = list(found_keys)

                missing_keys = expected_keys - found_keys
                if missing_keys:
                    result.issues.append(ValidationIssue(
                        level=ValidationStatus.WARNING,
                        category="structure",
                        message=f"Missing expected AnnData components: {missing_keys}",
                        details=f"Found keys: {found_keys}",
                        suggested_fix="Verify this is a complete AnnData object"
                    ))

                # Check X matrix
                if 'X' in found_keys:
                    x_shape = f['X'].shape if hasattr(f['X'], 'shape') else None
                    result.metadata['expression_matrix_shape'] = x_shape

                    if x_shape and len(x_shape) != 2:
                        result.issues.append(ValidationIssue(
                            level=ValidationStatus.ERROR,
                            category="structure",
                            message=f"Expression matrix X should be 2D, got shape {x_shape}",
                            suggested_fix="Ensure X contains a 2D expression matrix"
                        ))

        except Exception as e:
            result.issues.append(ValidationIssue(
                level=ValidationStatus.ERROR,
                category="integrity",
                message=f"Failed to read HDF5 file: {e}",
                suggested_fix="Check file integrity and HDF5 format"
            ))

    def _validate_anndata_with_library(
        self,
        file_path: Path,
        result: ValidationResult,
        validation_level: ValidationLevel
    ) -> None:
        """Validate AnnData using the anndata library."""
        try:
            import anndata as ad

            # Try to read the AnnData object
            adata = ad.read_h5ad(file_path)

            result.metadata['n_obs'] = adata.n_obs
            result.metadata['n_vars'] = adata.n_vars
            result.metadata['obs_columns'] = list(adata.obs.columns)
            result.metadata['var_columns'] = list(adata.var.columns)

            # Domain-specific validation for spatial transcriptomics
            if validation_level == ValidationLevel.DOMAIN:
                self._validate_spatial_anndata(adata, result)

        except Exception as e:
            result.issues.append(ValidationIssue(
                level=ValidationStatus.ERROR,
                category="integrity",
                message=f"Failed to read AnnData object: {e}",
                suggested_fix="Check AnnData file integrity and format"
            ))

    def _validate_spatial_anndata(self, adata, result: ValidationResult) -> None:
        """Validate AnnData for spatial transcriptomics specific content."""
        # Check for spatial coordinates
        spatial_keys = [key for key in adata.obsm.keys() if 'spatial' in key.lower()]
        if not spatial_keys:
            result.issues.append(ValidationIssue(
                level=ValidationStatus.WARNING,
                category="domain",
                message="No spatial coordinates found in obsm",
                details="Expected keys containing 'spatial' in adata.obsm",
                suggested_fix="Add spatial coordinates to adata.obsm for spatial analysis"
            ))
        else:
            result.metadata['spatial_keys'] = spatial_keys

            # Validate spatial coordinate structure
            for key in spatial_keys:
                coords = adata.obsm[key]
                if coords.shape[1] not in [2, 3]:
                    result.issues.append(ValidationIssue(
                        level=ValidationStatus.ERROR,
                        category="domain",
                        message=f"Spatial coordinates '{key}' should have 2 or 3 dimensions, got {coords.shape[1]}",
                        suggested_fix="Ensure spatial coordinates are 2D (x,y) or 3D (x,y,z)"
                    ))

                # Check for reasonable coordinate ranges
                if coords.min() < -1e6 or coords.max() > 1e6:
                    result.issues.append(ValidationIssue(
                        level=ValidationStatus.WARNING,
                        category="domain",
                        message=f"Spatial coordinates '{key}' have extreme values (min: {coords.min():.2f}, max: {coords.max():.2f})",
                        suggested_fix="Verify coordinate system and scaling"
                    ))

        # Check for reasonable gene expression data
        if adata.n_vars < 100:
            result.issues.append(ValidationIssue(
                level=ValidationStatus.WARNING,
                category="domain",
                message=f"Low number of genes/features: {adata.n_vars}",
                suggested_fix="Verify this contains the expected gene expression data"
            ))
        elif adata.n_vars > 50000:
            result.issues.append(ValidationIssue(
                level=ValidationStatus.WARNING,
                category="domain",
                message=f"Very high number of features: {adata.n_vars}",
                details="This may include non-gene features or need filtering",
                suggested_fix="Consider filtering to protein-coding genes if appropriate"
            ))

        # Check cell count reasonableness
        if adata.n_obs < 10:
            result.issues.append(ValidationIssue(
                level=ValidationStatus.WARNING,
                category="domain",
                message=f"Very low number of cells/spots: {adata.n_obs}",
                suggested_fix="Verify this contains sufficient spatial data for analysis"
            ))
        elif adata.n_obs > 1000000:
            result.issues.append(ValidationIssue(
                level=ValidationStatus.WARNING,
                category="domain",
                message=f"Very high number of observations: {adata.n_obs}",
                details="This may impact memory usage and processing time",
                suggested_fix="Consider data subsampling or chunked processing"
            ))

        # Check for spatial metadata in uns
        spatial_uns_keys = [key for key in adata.uns.keys() if 'spatial' in key.lower()]
        if spatial_uns_keys:
            result.metadata['spatial_uns_keys'] = spatial_uns_keys

            # Validate spatial metadata structure
            for key in spatial_uns_keys:
                spatial_info = adata.uns[key]
                if isinstance(spatial_info, dict):
                    # Check for common spatial metadata fields
                    expected_fields = ['images', 'scalefactors', 'use_quality']
                    missing_fields = [f for f in expected_fields if f not in spatial_info]
                    if missing_fields:
                        result.issues.append(ValidationIssue(
                            level=ValidationStatus.WARNING,
                            category="domain",
                            message=f"Spatial metadata '{key}' missing fields: {missing_fields}",
                            suggested_fix="Add standard spatial metadata fields for compatibility"
                        ))

        # Check expression matrix properties
        if hasattr(adata.X, 'data'):  # Sparse matrix
            # Check sparsity
            sparsity = 1.0 - (adata.X.nnz / (adata.n_obs * adata.n_vars))
            result.metadata['expression_sparsity'] = sparsity

            if sparsity > 0.99:
                result.issues.append(ValidationIssue(
                    level=ValidationStatus.WARNING,
                    category="domain",
                    message=f"Very sparse expression matrix ({sparsity:.1%} zeros)",
                    details="This may indicate over-filtering or technical issues",
                    suggested_fix="Verify filtering parameters and data quality"
                ))
            elif sparsity < 0.5:
                result.issues.append(ValidationIssue(
                    level=ValidationStatus.WARNING,
                    category="domain",
                    message=f"Unusually dense expression matrix ({sparsity:.1%} zeros)",
                    details="This may indicate normalized/imputed data",
                    suggested_fix="Verify if raw counts are expected for downstream analysis"
                ))

        # Check for gene symbols vs IDs
        if hasattr(adata.var, 'index'):
            gene_names = adata.var.index.tolist()[:100]  # Sample first 100

            # Check if using Ensembl IDs
            ensembl_pattern = any(name.startswith('ENS') for name in gene_names)
            if ensembl_pattern:
                result.metadata['gene_id_type'] = 'ensembl'
            else:
                # Check if using gene symbols
                symbol_pattern = any(name.replace('-', '').replace('_', '').isalpha() for name in gene_names)
                if symbol_pattern:
                    result.metadata['gene_id_type'] = 'symbol'
                else:
                    result.metadata['gene_id_type'] = 'unknown'
                    result.issues.append(ValidationIssue(
                        level=ValidationStatus.WARNING,
                        category="domain",
                        message="Gene identifier format not recognized",
                        details="Expected gene symbols or Ensembl IDs",
                        suggested_fix="Ensure gene identifiers follow standard conventions"
                    ))

    def _check_expression_table_structure(self, table_path: Path, result: ValidationResult) -> None:
        """Check if a table component looks like gene expression data."""
        try:
            # Look for zarr array structure
            zarray_file = table_path / '.zarray'
            if zarray_file.exists():
                with open(zarray_file, 'r') as f:
                    zarray_info = json.load(f)
                    shape = zarray_info.get('shape', [])

                    if len(shape) == 2:
                        result.metadata['table_shape'] = shape
                        if shape[1] < 100:  # Assuming genes should be > 100
                            result.issues.append(ValidationIssue(
                                level=ValidationStatus.WARNING,
                                category="domain",
                                message=f"Low number of features in expression table: {shape[1]}",
                                suggested_fix="Verify this contains gene expression data"
                            ))
        except Exception as e:
            logger.debug(f"Could not analyze table structure: {e}")

    def _check_spatial_images(self, images_path: Path, result: ValidationResult) -> None:
        """Check spatial images for reasonable properties."""
        try:
            # Count image arrays
            image_dirs = [d for d in images_path.iterdir() if d.is_dir()]
            result.metadata['image_count'] = len(image_dirs)

            if len(image_dirs) == 0:
                result.issues.append(ValidationIssue(
                    level=ValidationStatus.WARNING,
                    category="domain",
                    message="No images found in images component",
                    suggested_fix="Ensure spatial images are included for spatial analysis"
                ))
        except Exception as e:
            logger.debug(f"Could not analyze images: {e}")

    def _validate_table_spatial_annotations(self, table_path: Path, result: ValidationResult) -> None:
        """Validate spatial annotations in table component."""
        try:
            # Check for region/instance_id columns that link to spatial elements
            zattrs_file = table_path / '.zattrs'
            if zattrs_file.exists():
                with open(zattrs_file, 'r') as f:
                    attrs = json.load(f)

                    # Check for spatial annotation keys
                    if 'spatialdata_attrs' in attrs:
                        spatial_attrs = attrs['spatialdata_attrs']
                        if 'region' in spatial_attrs:
                            result.metadata['has_region_annotations'] = True
                        if 'instance_id' in spatial_attrs:
                            result.metadata['has_instance_annotations'] = True

                    # Warn if no spatial linking found
                    if not any(key in attrs.get('spatialdata_attrs', {}) for key in ['region', 'instance_id']):
                        result.issues.append(ValidationIssue(
                            level=ValidationStatus.WARNING,
                            category="domain",
                            message="Table lacks spatial region annotations",
                            details="No 'region' or 'instance_id' found in spatial attributes",
                            suggested_fix="Add spatial region annotations to link expression data to spatial elements"
                        ))
        except Exception as e:
            logger.debug(f"Could not analyze table spatial annotations: {e}")

    def _validate_image_spatial_properties(self, images_path: Path, result: ValidationResult) -> None:
        """Validate spatial properties of images."""
        try:
            for image_dir in images_path.iterdir():
                if image_dir.is_dir():
                    # Check for coordinate transformation information
                    zattrs_file = image_dir / '.zattrs'
                    if zattrs_file.exists():
                        with open(zattrs_file, 'r') as f:
                            attrs = json.load(f)

                            # Check for coordinate system information
                            if 'spatialdata_attrs' in attrs:
                                spatial_attrs = attrs['spatialdata_attrs']
                                if 'transform' not in spatial_attrs:
                                    result.issues.append(ValidationIssue(
                                        level=ValidationStatus.WARNING,
                                        category="domain",
                                        message=f"Image '{image_dir.name}' lacks coordinate transformation",
                                        suggested_fix="Add coordinate transformation for proper spatial alignment"
                                    ))
                                else:
                                    result.metadata['has_coordinate_transforms'] = True
        except Exception as e:
            logger.debug(f"Could not analyze image spatial properties: {e}")

    def _validate_spatial_points(self, points_path: Path, result: ValidationResult) -> None:
        """Validate spatial points data."""
        try:
            # Check points structure
            for points_dir in points_path.iterdir():
                if points_dir.is_dir():
                    # Check for coordinate arrays
                    coord_files = list(points_dir.glob('**/0'))  # Zarr array chunks
                    if not coord_files:
                        result.issues.append(ValidationIssue(
                            level=ValidationStatus.ERROR,
                            category="structure",
                            message=f"Points '{points_dir.name}' missing coordinate data",
                            suggested_fix="Ensure points contain coordinate arrays"
                        ))

                    # Check for reasonable coordinate dimensions
                    zarray_file = points_dir / '.zarray'
                    if zarray_file.exists():
                        with open(zarray_file, 'r') as f:
                            zarray_info = json.load(f)
                            shape = zarray_info.get('shape', [])

                            if len(shape) >= 2 and shape[1] not in [2, 3]:
                                result.issues.append(ValidationIssue(
                                    level=ValidationStatus.WARNING,
                                    category="domain",
                                    message=f"Points '{points_dir.name}' have unusual coordinate dimensions: {shape[1]}",
                                    details="Expected 2D (x,y) or 3D (x,y,z) coordinates",
                                    suggested_fix="Verify coordinate system and dimensions"
                                ))
        except Exception as e:
            logger.debug(f"Could not analyze spatial points: {e}")

    def _validate_spatial_shapes(self, shapes_path: Path, result: ValidationResult) -> None:
        """Validate spatial shapes data."""
        try:
            # Check shapes structure
            for shapes_dir in shapes_path.iterdir():
                if shapes_dir.is_dir():
                    # Check for geometry data
                    geom_files = list(shapes_dir.rglob('*.parquet'))
                    if not geom_files:
                        # Check for zarr-based geometry
                        zarray_files = list(shapes_dir.rglob('.zarray'))
                        if not zarray_files:
                            result.issues.append(ValidationIssue(
                                level=ValidationStatus.ERROR,
                                category="structure",
                                message=f"Shapes '{shapes_dir.name}' missing geometry data",
                                suggested_fix="Ensure shapes contain geometry information"
                            ))

                    # Check for spatial attributes
                    zattrs_file = shapes_dir / '.zattrs'
                    if zattrs_file.exists():
                        with open(zattrs_file, 'r') as f:
                            attrs = json.load(f)
                            if 'spatialdata_attrs' not in attrs:
                                result.issues.append(ValidationIssue(
                                    level=ValidationStatus.WARNING,
                                    category="domain",
                                    message=f"Shapes '{shapes_dir.name}' lack spatial metadata",
                                    suggested_fix="Add spatial attributes for proper coordinate handling"
                                ))
        except Exception as e:
            logger.debug(f"Could not analyze spatial shapes: {e}")

    def _validate_component_relationships(self, file_path: Path, result: ValidationResult, components: set) -> None:
        """Validate relationships between SpatialData components."""
        try:
            # Check if table has corresponding spatial elements
            if 'table' in components:
                table_path = file_path / 'table'
                zattrs_file = table_path / '.zattrs'

                if zattrs_file.exists():
                    with open(zattrs_file, 'r') as f:
                        attrs = json.load(f)
                        spatial_attrs = attrs.get('spatialdata_attrs', {})

                        # Check if referenced regions exist
                        if 'region' in spatial_attrs:
                            region_refs = spatial_attrs['region']
                            if isinstance(region_refs, list):
                                for region_ref in region_refs:
                                    # Check if referenced component exists
                                    ref_component = region_ref.split('/')[0] if '/' in region_ref else region_ref
                                    if ref_component not in components:
                                        result.issues.append(ValidationIssue(
                                            level=ValidationStatus.ERROR,
                                            category="integrity",
                                            message=f"Table references missing component: {ref_component}",
                                            suggested_fix=f"Ensure component '{ref_component}' exists or update table references"
                                        ))

            # Check coordinate system consistency
            coord_systems = set()
            for component in ['images', 'points', 'shapes']:
                if component in components:
                    comp_path = file_path / component
                    for item_dir in comp_path.iterdir():
                        if item_dir.is_dir():
                            zattrs_file = item_dir / '.zattrs'
                            if zattrs_file.exists():
                                with open(zattrs_file, 'r') as f:
                                    attrs = json.load(f)
                                    spatial_attrs = attrs.get('spatialdata_attrs', {})
                                    if 'transform' in spatial_attrs:
                                        # Extract coordinate system info
                                        transform = spatial_attrs['transform']
                                        if isinstance(transform, dict):
                                            coord_systems.add(str(transform.keys()))

            if len(coord_systems) > 1:
                result.issues.append(ValidationIssue(
                    level=ValidationStatus.WARNING,
                    category="domain",
                    message="Multiple coordinate systems detected",
                    details=f"Found coordinate systems: {coord_systems}",
                    suggested_fix="Ensure consistent coordinate systems across spatial components"
                ))
            elif len(coord_systems) == 1:
                result.metadata['coordinate_systems'] = list(coord_systems)

        except Exception as e:
            logger.debug(f"Could not analyze component relationships: {e}")

    def validate_multiple_files(
        self,
        file_paths: List[Union[str, Path]],
        validation_level: ValidationLevel = ValidationLevel.STRUCTURE
    ) -> List[ValidationResult]:
        """Validate multiple spatial data files.

        Args:
            file_paths: List of file paths to validate
            validation_level: Level of validation to perform

        Returns:
            List of ValidationResult objects
        """
        results = []
        for file_path in file_paths:
            try:
                result = self.validate_file(file_path, validation_level)
                results.append(result)
            except Exception as e:
                # Create error result for failed validation
                error_result = ValidationResult(
                    file_path=str(file_path),
                    format_type="unknown",
                    is_valid=False,
                    validation_level=validation_level,
                    issues=[ValidationIssue(
                        level=ValidationStatus.CRITICAL,
                        category="validation_error",
                        message=f"Validation failed with exception: {e}",
                        suggested_fix="Check file accessibility and format"
                    )],
                    metadata={},
                    file_size=0
                )
                results.append(error_result)

        return results

    def get_validation_summary(self, results: List[ValidationResult]) -> Dict[str, Any]:
        """Generate a summary of validation results.

        Args:
            results: List of validation results

        Returns:
            Summary dictionary with statistics and issues
        """
        total_files = len(results)
        valid_files = sum(1 for r in results if r.is_valid)
        files_with_warnings = sum(1 for r in results if r.has_warnings)
        files_with_errors = sum(1 for r in results if r.has_errors)

        format_counts = {}
        for result in results:
            format_counts[result.format_type] = format_counts.get(result.format_type, 0) + 1

        all_issues = []
        for result in results:
            all_issues.extend(result.issues)

        issue_categories = {}
        for issue in all_issues:
            issue_categories[issue.category] = issue_categories.get(issue.category, 0) + 1

        return {
            'total_files': total_files,
            'valid_files': valid_files,
            'files_with_warnings': files_with_warnings,
            'files_with_errors': files_with_errors,
            'success_rate': valid_files / total_files if total_files > 0 else 0,
            'format_distribution': format_counts,
            'issue_categories': issue_categories,
            'total_issues': len(all_issues),
            'dependencies_available': self.dependencies
        }

    def assess_spatial_data_quality(
        self,
        file_path: Union[str, Path],
        validation_level: ValidationLevel = ValidationLevel.DOMAIN
    ) -> Dict[str, Any]:
        """Perform comprehensive spatial data quality assessment.

        Args:
            file_path: Path to the spatial data file
            validation_level: Level of validation to perform

        Returns:
            Dictionary with comprehensive quality assessment
        """
        result = self.validate_file(file_path, validation_level)

        quality_assessment = {
            'file_path': str(file_path),
            'format_type': result.format_type,
            'overall_quality': 'good' if result.is_valid and not result.has_warnings else
                             'acceptable' if result.is_valid else 'poor',
            'validation_passed': result.is_valid,
            'issues_found': len(result.issues),
            'critical_issues': len([i for i in result.issues if i.level == ValidationStatus.CRITICAL]),
            'error_issues': len([i for i in result.issues if i.level == ValidationStatus.ERROR]),
            'warning_issues': len([i for i in result.issues if i.level == ValidationStatus.WARNING]),
            'file_size_mb': result.file_size / (1024 * 1024) if result.file_size > 0 else 0,
            'metadata': result.metadata,
            'recommendations': []
        }

        # Generate quality-based recommendations
        if result.has_errors:
            quality_assessment['recommendations'].append(
                "Address critical errors before using this data for analysis"
            )

        if result.format_type == 'anndata':
            if 'n_obs' in result.metadata and result.metadata['n_obs'] < 100:
                quality_assessment['recommendations'].append(
                    "Consider combining with additional samples for robust analysis"
                )

            if 'expression_sparsity' in result.metadata:
                sparsity = result.metadata['expression_sparsity']
                if sparsity > 0.95:
                    quality_assessment['recommendations'].append(
                        "High sparsity detected - consider quality filtering parameters"
                    )

        elif result.format_type == 'spatialdata':
            if 'components' in result.metadata:
                components = result.metadata['components']
                if 'table' not in components:
                    quality_assessment['recommendations'].append(
                        "Add expression table component for transcriptomics analysis"
                    )
                if not any(comp in components for comp in ['images', 'points', 'shapes']):
                    quality_assessment['recommendations'].append(
                        "Add spatial coordinate components for spatial analysis"
                    )

        # Performance recommendations
        if quality_assessment['file_size_mb'] > 1000:  # > 1GB
            quality_assessment['recommendations'].append(
                "Large file size - consider chunked processing or data subsampling"
            )

        # Dependency recommendations
        missing_deps = [dep for dep, available in self.dependencies.items() if not available]
        if missing_deps:
            quality_assessment['recommendations'].append(
                f"Install optional dependencies for enhanced validation: {', '.join(missing_deps)}"
            )

        return quality_assessment

    def validate_spatial_data_compatibility(
        self,
        file_paths: List[Union[str, Path]],
        validation_level: ValidationLevel = ValidationLevel.STRUCTURE
    ) -> Dict[str, Any]:
        """Check compatibility between multiple spatial data files.

        Args:
            file_paths: List of file paths to check for compatibility
            validation_level: Level of validation to perform

        Returns:
            Dictionary with compatibility assessment
        """
        results = self.validate_multiple_files(file_paths, validation_level)

        compatibility = {
            'total_files': len(results),
            'valid_files': len([r for r in results if r.is_valid]),
            'format_compatibility': True,
            'dimension_compatibility': True,
            'coordinate_compatibility': True,
            'issues': [],
            'recommendations': []
        }

        # Check format compatibility
        formats = [r.format_type for r in results if r.is_valid]
        unique_formats = set(formats)

        if len(unique_formats) > 1:
            compatibility['format_compatibility'] = False
            compatibility['issues'].append(
                f"Mixed formats detected: {', '.join(unique_formats)}"
            )
            compatibility['recommendations'].append(
                "Convert all files to a common format for joint analysis"
            )

        # Check dimension compatibility
        shapes = []
        for result in results:
            if result.is_valid and 'expression_matrix_shape' in result.metadata:
                shapes.append(result.metadata['expression_matrix_shape'])
            elif result.is_valid and 'n_obs' in result.metadata and 'n_vars' in result.metadata:
                shapes.append([result.metadata['n_obs'], result.metadata['n_vars']])

        if shapes:
            gene_counts = [s[1] if len(s) > 1 else None for s in shapes]
            gene_counts = [g for g in gene_counts if g is not None]

            if len(set(gene_counts)) > 1:
                compatibility['dimension_compatibility'] = False
                compatibility['issues'].append(
                    f"Different gene counts: {gene_counts}"
                )
                compatibility['recommendations'].append(
                    "Harmonize gene sets or filter to common genes"
                )

        # Check spatial coordinate compatibility
        spatial_info = []
        for result in results:
            if result.is_valid:
                has_spatial = (
                    'spatial_keys' in result.metadata or
                    'components' in result.metadata and
                    any(comp in ['images', 'points', 'shapes'] for comp in result.metadata['components'])
                )
                spatial_info.append(has_spatial)

        if spatial_info and not all(spatial_info):
            compatibility['coordinate_compatibility'] = False
            compatibility['issues'].append(
                "Some files lack spatial coordinate information"
            )
            compatibility['recommendations'].append(
                "Ensure all files contain spatial coordinates for spatial analysis"
            )

        # Overall compatibility score
        compatibility_score = sum([
            compatibility['format_compatibility'],
            compatibility['dimension_compatibility'],
            compatibility['coordinate_compatibility']
        ]) / 3.0

        compatibility['compatibility_score'] = compatibility_score
        compatibility['overall_compatible'] = compatibility_score >= 0.67

        return compatibility
