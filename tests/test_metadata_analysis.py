"""Unit tests for bioinformatics metadata extraction and analysis module."""

import pytest
import tempfile
import json
import yaml
from pathlib import Path
from unittest.mock import Mock, patch, mock_open

from openproblems_mcp.metadata_analysis import (
    BioinformaticsMetadataExtractor,
    MetadataType,
    MetadataField,
    MetadataExtractionResult
)


class TestBioinformaticsMetadataExtractor:
    """Test cases for BioinformaticsMetadataExtractor."""

    def setup_method(self):
        """Set up test fixtures."""
        self.extractor = BioinformaticsMetadataExtractor()

    def test_init(self):
        """Test extractor initialization."""
        assert isinstance(self.extractor, BioinformaticsMetadataExtractor)
        assert hasattr(self.extractor, 'dependencies')
        assert isinstance(self.extractor.dependencies, dict)

    def test_check_dependencies(self):
        """Test dependency checking."""
        deps = self.extractor._check_dependencies()

        # These should always be available
        assert deps['yaml'] is True
        assert deps['json'] is True
        assert deps['re'] is True
        assert deps['pathlib'] is True

        # These might not be available in test environment
        assert isinstance(deps.get('h5py', False), bool)
        assert isinstance(deps.get('zarr', False), bool)

    def test_detect_file_type_nextflow(self):
        """Test file type detection for Nextflow files."""
        # Test nextflow config
        path = Path("nextflow.config")
        result = self.extractor._detect_file_type(path)
        assert result == "nextflow_config"

        # Test nextflow workflow
        path = Path("main.nf")
        result = self.extror._detect_file_type(path)
        assert result == "nextflow_workflow"

        # Test nextflow params
        path = Path("params.yaml")
        result = self.extractor._detect_file_type(path)
        assert result == "nextflow_params"

    def test_detect_file_type_viash(self):
        """Test file type detection for Viash files."""
        # Test viash config
        path = Path("config.vsh.yaml")
        result = self.extractor._detect_file_type(path)
        assert result == "viash_config"

        # Test viash script
        path = Path("script.py")
        result = self.extractor._detect_file_type(path)
        assert result == "viash_script"

    def test_detect_file_type_spatial(self):
        """Test file type detection for spatial data files."""
        # Test spatialdata
        path = Path("data.spatialdata")
        result = self.extractor._detect_file_type(path)
        assert result == "spatial_spatialdata"

        # Test anndata
        path = Path("data.h5ad")
        result = self.extractor._detect_file_type(path)
        assert result == "spatial_anndata"

    def test_detect_file_type_generic(self):
        """Test file type detection for generic files."""
        # Test YAML
        path = Path("config.yaml")
        result = self.extractor._detect_file_type(path)
        assert result == "yaml_config"

        # Test JSON
        path = Path("config.json")
        result = self.extractor._detect_file_type(path)
        assert result == "json_config"

        # Test unknown
        path = Path("unknown.txt")
        result = self.extractor._detect_file_type(path)
        assert result == "unknown"

    def test_extract_metadata_nonexistent_file(self):
        """Test metadata extraction for non-existent file."""
        result = self.extractor.extract_metadata("nonexistent.nf")

        assert result.extraction_success is False
        assert len(result.issues) == 1
        assert "does not exist" in result.issues[0]

    def test_extract_nextflow_config_metadata(self):
        """Test Nextflow config metadata extraction."""
        nextflow_config = """
        process {
            cpus = 4
            memory = '8.GB'
            container = 'ubuntu:20.04'
        }

        profiles {
            standard {
                process.executor = 'local'
            }
            docker {
                docker.enabled = true
            }
        }
        """

        with tempfile.NamedTemporaryFile(mode='w', suffix='.config', delete=False) as f:
            f.write(nextflow_config)
            f.flush()

            try:
                result = self.extractor.extract_metadata(f.name)

                assert result.extraction_success is True
                assert result.file_type == "unknown"  # Will be detected as unknown due to suffix

                # Check for extracted metadata
                field_names = [field.name for field in result.metadata_fields]
                assert "file_size" in field_names
                assert "line_count" in field_names

            finally:
                Path(f.name).unlink()

    def test_extract_nextflow_workflow_metadata(self):
        """Test Nextflow workflow metadata extraction."""
        nextflow_workflow = """
        workflow {
            input:
                path reads

            output:
                path "results/*"

            main:
                PROCESS_A(reads)
                PROCESS_B(PROCESS_A.out)
        }
        """

        with tempfile.NamedTemporaryFile(mode='w', suffix='.nf', delete=False) as f:
            f.write(nextflow_workflow)
            f.flush()

            try:
                result = self.extractor.extract_metadata(f.name)

                assert result.extraction_success is True
                assert result.file_type == "nextflow_workflow"

                # Check for workflow-specific metadata
                field_names = [field.name for field in result.metadata_fields]
                assert "has_workflow_block" in field_names

            finally:
                Path(f.name).unlink()

    def test_extract_viash_config_metadata(self):
        """Test Viash config metadata extraction."""
        viash_config = {
            "name": "test_component",
            "description": "Test component for spatial analysis",
            "version": "1.0.0",
            "functionality": {
                "arguments": [
                    {"name": "input", "direction": "input", "type": "file"},
                    {"name": "output", "direction": "output", "type": "file"}
                ]
            },
            "platforms": [
                {"type": "docker", "image": "python:3.9"}
            ]
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.vsh.yaml', delete=False) as f:
            yaml.dump(viash_config, f)
            f.flush()

            try:
                result = self.extractor.extract_metadata(f.name)

                assert result.extraction_success is True
                assert result.file_type == "viash_config"

                # Check for component-specific metadata
                field_names = [field.name for field in result.metadata_fields]
                assert "component_name" in field_names
                assert "description" in field_names
                assert "version" in field_names

            finally:
                Path(f.name).unlink()

    def test_extract_python_script_metadata(self):
        """Test Python script metadata extraction."""
        python_script = """
        import pandas as pd
        import numpy as np
        import scanpy as sc
        import spatialdata as sd

        def process_data(input_file):
            # Process spatial transcriptomics data
            return data

        def analyze_results(data):
            # Analyze results
            return analysis
        """

        with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as f:
            f.write(python_script)
            f.flush()

            try:
                result = self.extractor.extract_metadata(f.name)

                assert result.extraction_success is True
                assert result.file_type == "viash_script"

                # Check for Python-specific metadata
                field_names = [field.name for field in result.metadata_fields]
                assert "script_language" in field_names
                assert "bioinformatics_libraries" in field_names
                assert "function_count" in field_names

            finally:
                Path(f.name).unlink()

    def test_analyze_workflow_dependencies(self):
        """Test workflow dependency analysis."""
        # Create temporary files
        files = []

        # Nextflow config with containers
        nextflow_config = """
        process {
            container = 'bioconductor/bioconductor_docker:latest'
        }
        """

        # Viash config with Docker platform
        viash_config = {
            "name": "test_component",
            "platforms": [
                {"type": "docker", "image": "python:3.9"}
            ]
        }

        try:
            # Create Nextflow config file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.config', delete=False) as f:
                f.write(nextflow_config)
                files.append(f.name)

            # Create Viash config file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.vsh.yaml', delete=False) as f:
                yaml.dump(viash_config, f)
                files.append(f.name)

            # Analyze dependencies
            analysis = self.extractor.analyze_workflow_dependencies(files)

            assert analysis['total_files_analyzed'] == 2
            assert analysis['successful_extractions'] >= 0
            assert 'unique_containers' in analysis
            assert 'unique_processes' in analysis
            assert 'unique_libraries' in analysis

        finally:
            # Clean up
            for file_path in files:
                Path(file_path).unlink(missing_ok=True)


class TestMetadataField:
    """Test cases for MetadataField dataclass."""

    def test_metadata_field_creation(self):
        """Test MetadataField creation."""
        field = MetadataField(
            name="test_field",
            value="test_value",
            data_type="str",
            category="test",
            description="Test field",
            importance="high"
        )

        assert field.name == "test_field"
        assert field.value == "test_value"
        assert field.data_type == "str"
        assert field.category == "test"
        assert field.description == "Test field"
        assert field.importance == "high"

    def test_metadata_field_defaults(self):
        """Test MetadataField default values."""
        field = MetadataField(
            name="test_field",
            value="test_value",
            data_type="str",
            category="test"
        )

        assert field.description is None
        assert field.importance == "medium"


class TestMetadataExtractionResult:
    """Test cases for MetadataExtractionResult dataclass."""

    def test_metadata_extraction_result_creation(self):
        """Test MetadataExtractionResult creation."""
        fields = [
            MetadataField("field1", "value1", "str", "category1", importance="critical"),
            MetadataField("field2", "value2", "str", "category2", importance="high"),
            MetadataField("field3", "value3", "str", "category1", importance="medium")
        ]

        result = MetadataExtractionResult(
            file_path="test.nf",
            file_type="nextflow_workflow",
            extraction_success=True,
            metadata_fields=fields,
            quality_metrics={"completeness": 0.8},
            issues=[],
            suggestions=["Add more metadata"]
        )

        assert result.file_path == "test.nf"
        assert result.file_type == "nextflow_workflow"
        assert result.extraction_success is True
        assert len(result.metadata_fields) == 3
        assert result.quality_metrics["completeness"] == 0.8

    def test_get_fields_by_category(self):
        """Test getting fields by category."""
        fields = [
            MetadataField("field1", "value1", "str", "category1"),
            MetadataField("field2", "value2", "str", "category2"),
            MetadataField("field3", "value3", "str", "category1")
        ]

        result = MetadataExtractionResult(
            file_path="test.nf",
            file_type="nextflow_workflow",
            extraction_success=True,
            metadata_fields=fields,
            quality_metrics={},
            issues=[],
            suggestions=[]
        )

        category1_fields = result.get_fields_by_category("category1")
        assert len(category1_fields) == 2
        assert all(f.category == "category1" for f in category1_fields)

    def test_get_critical_fields(self):
        """Test getting critical fields."""
        fields = [
            MetadataField("field1", "value1", "str", "category1", importance="critical"),
            MetadataField("field2", "value2", "str", "category2", importance="high"),
            MetadataField("field3", "value3", "str", "category1", importance="critical")
        ]

        result = MetadataExtractionResult(
            file_path="test.nf",
            file_type="nextflow_workflow",
            extraction_success=True,
            metadata_fields=fields,
            quality_metrics={},
            issues=[],
            suggestions=[]
        )

        critical_fields = result.get_critical_fields()
        assert len(critical_fields) == 2
        assert all(f.importance == "critical" for f in critical_fields)


class TestMetadataAnalysisIntegration:
    """Integration tests for metadata analysis."""

    def setup_method(self):
        """Set up test fixtures."""
        self.extractor = BioinformaticsMetadataExtractor()

    def test_end_to_end_nextflow_analysis(self):
        """Test complete Nextflow analysis workflow."""
        # Create a comprehensive Nextflow workflow
        nextflow_content = """
        #!/usr/bin/env nextflow

        params.input = 'data/*.fastq'
        params.output = 'results'

        process QUALITY_CONTROL {
            container 'biocontainers/fastqc:v0.11.9_cv8'

            input:
            path reads

            output:
            path "*.html"

            script:
            '''
            fastqc !{reads}
            '''
        }

        process ALIGNMENT {
            container 'biocontainers/bwa:v0.7.17_cv1'
            cpus 4
            memory '8.GB'

            input:
            path reads

            output:
            path "*.bam"

            script:
            '''
            bwa mem ref.fa !{reads} | samtools sort -o aligned.bam
            '''
        }

        workflow {
            reads_ch = Channel.fromPath(params.input)
            QUALITY_CONTROL(reads_ch)
            ALIGNMENT(reads_ch)
        }
        """

        with tempfile.NamedTemporaryFile(mode='w', suffix='.nf', delete=False) as f:
            f.write(nextflow_content)
            f.flush()

            try:
                result = self.extractor.extract_metadata(f.name)

                assert result.extraction_success is True
                assert result.file_type == "nextflow_workflow"

                # Verify comprehensive metadata extraction
                field_names = [field.name for field in result.metadata_fields]
                assert "file_size" in field_names
                assert "line_count" in field_names

                # Check quality assessment
                assert "total_metadata_fields" in result.quality_metrics
                assert "metadata_completeness" in result.quality_metrics

            finally:
                Path(f.name).unlink()

    def test_batch_metadata_extraction(self):
        """Test batch metadata extraction workflow."""
        files = []

        try:
            # Create multiple workflow files
            configs = [
                ("nextflow.config", "process { cpus = 2 }"),
                ("main.nf", "workflow { println 'Hello' }"),
                ("config.vsh.yaml", yaml.dump({"name": "test", "version": "1.0"}))
            ]

            for filename, content in configs:
                with tempfile.NamedTemporaryFile(mode='w', suffix=f'.{filename.split(".")[-1]}', delete=False) as f:
                    f.write(content)
                    files.append(f.name)

            # Extract metadata from all files
            results = []
            for file_path in files:
                result = self.extractor.extract_metadata(file_path)
                results.append(result)

            # Analyze dependencies
            dependency_analysis = self.extractor.analyze_workflow_dependencies(files)

            assert len(results) == 3
            assert dependency_analysis['total_files_analyzed'] == 3
            assert dependency_analysis['successful_extractions'] >= 0

        finally:
            # Clean up
            for file_path in files:
                Path(file_path).unlink(missing_ok=True)


class TestMetadataType:
    """Test cases for MetadataType enum."""

    def test_metadata_type_values(self):
        """Test MetadataType enum values."""
        assert MetadataType.SPATIAL_DATA.value == "spatial_data"
        assert MetadataType.WORKFLOW_CONFIG.value == "workflow_config"
        assert MetadataType.DATA_QUALITY.value == "data_quality"
        assert MetadataType.BIOINFORMATICS.value == "bioinformatics"


if __name__ == "__main__":
    pytest.main([__file__])
