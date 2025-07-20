"""Bioinformatics metadata extraction and analysis module.

This module provides comprehensive metadata extraction and analysis for spatial
transcriptomics data formats and bioinformatics workflow configurations.
"""

import logging
import os
import json
import yaml
import re
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple
from dataclasses import dataclass
from enum import Enum

logger = logging.getLogger(__name__)


class MetadataType(Enum):
    """Types of metadata that can be extracted."""
    SPATIAL_DATA = "spatial_data"
    WORKFLOW_CONFIG = "workflow_config"
    DATA_QUALITY = "data_quality"
    BIOINFORMATICS = "bioinformatics"


@dataclass
class MetadataField:
    """Represents a metadata field with its properties."""
    name: str
    value: Any
    data_type: str
    category: str
    description: Optional[str] = None
    importance: str = "medium"  # low, medium, high, critical


@dataclass
class MetadataExtractionResult:
    """Result of metadata extraction."""
    file_path: str
    file_type: str
    extraction_success: bool
    metadata_fields: List[MetadataField]
    quality_metrics: Dict[str, Any]
    issues: List[str]
    suggestions: List[str]

    def get_fields_by_category(self, category: str) -> List[MetadataField]:
        """Get metadata fields by category."""
        return [field for field in self.metadata_fields if field.category == category]

    def get_critical_fields(self) -> List[MetadataField]:
        """Get critical metadata fields."""
        return [field for field in self.metadata_fields if field.importance == "critical"]


class BioinformaticsMetadataExtractor:
    """Extractor for bioinformatics metadata from various file formats."""

    # File type patterns
    NEXTFLOW_PATTERNS = {
        'config': [r'nextflow\.config$', r'\.nf\.config$'],
        'workflow': [r'main\.nf$', r'.*\.nf$'],
        'params': [r'params\.ya?ml$', r'nextflow_params\.ya?ml$']
    }

    VIASH_PATTERNS = {
        'config': [r'config\.vsh\.ya?ml$', r'.*\.vsh\.ya?ml$'],
        'script': [r'script\.py$', r'script\.R$', r'script\.sh$']
    }

    SPATIAL_DATA_PATTERNS = {
        'spatialdata': [r'.*\.spatialdata$', r'.*\.zarr$'],
        'anndata': [r'.*\.h5ad$', r'.*\.h5$'],
        'images': [r'.*\.(tiff?|png|jpg|jpeg)$'],
        'coordinates': [r'.*coordinates.*\.(csv|tsv|txt)$']
    }

    def __init__(self):
        """Initialize the metadata extractor."""
        self._check_dependencies()

    def _check_dependencies(self) -> Dict[str, bool]:
        """Check availability of optional libraries."""
        dependencies = {}

        for lib in ['yaml', 'json', 're', 'pathlib']:
            try:
                __import__(lib)
                dependencies[lib] = True
            except ImportError:
                dependencies[lib] = False

        # Check for specialized libraries
        for lib in ['h5py', 'zarr', 'anndata', 'spatialdata']:
            try:
                __import__(lib)
                dependencies[lib] = True
            except ImportError:
                dependencies[lib] = False

        self.dependencies = dependencies
        return dependencies

    def extract_metadata(self, file_path: Union[str, Path]) -> MetadataExtractionResult:
        """Extract metadata from a bioinformatics file.

        Args:
            file_path: Path to the file to analyze

        Returns:
            MetadataExtractionResult with extracted information
        """
        file_path = Path(file_path)

        # Initialize result
        result = MetadataExtractionResult(
            file_path=str(file_path),
            file_type="unknown",
            extraction_success=False,
            metadata_fields=[],
            quality_metrics={},
            issues=[],
            suggestions=[]
        )

        # Check file existence
        if not file_path.exists():
            result.issues.append(f"File does not exist: {file_path}")
            return result

        # Detect file type
        file_type = self._detect_file_type(file_path)
        result.file_type = file_type

        try:
            # Extract metadata based on file type
            if file_type.startswith('nextflow'):
                self._extract_nextflow_metadata(file_path, result)
            elif file_type.startswith('viash'):
                self._extract_viash_metadata(file_path, result)
            elif file_type.startswith('spatial'):
                self._extract_spatial_data_metadata(file_path, result)
            else:
                # Try generic extraction
                self._extract_generic_metadata(file_path, result)

            result.extraction_success = True

        except Exception as e:
            result.issues.append(f"Metadata extraction failed: {e}")
            logger.error(f"Metadata extraction failed for {file_path}: {e}")

        # Perform quality assessment
        self._assess_metadata_quality(result)

        return result

    def _detect_file_type(self, file_path: Path) -> str:
        """Detect the type of bioinformatics file."""
        filename = file_path.name.lower()

        # Check Nextflow patterns
        for nf_type, patterns in self.NEXTFLOW_PATTERNS.items():
            for pattern in patterns:
                if re.search(pattern, filename):
                    return f"nextflow_{nf_type}"

        # Check Viash patterns
        for viash_type, patterns in self.VIASH_PATTERNS.items():
            for pattern in patterns:
                if re.search(pattern, filename):
                    return f"viash_{viash_type}"

        # Check spatial data patterns
        for spatial_type, patterns in self.SPATIAL_DATA_PATTERNS.items():
            for pattern in patterns:
                if re.search(pattern, filename):
                    return f"spatial_{spatial_type}"

        # Check by extension
        suffix = file_path.suffix.lower()
        if suffix in ['.yaml', '.yml']:
            return "yaml_config"
        elif suffix == '.json':
            return "json_config"
        elif suffix == '.nf':
            return "nextflow_workflow"

        return "unknown"

    def _extract_nextflow_metadata(self, file_path: Path, result: MetadataExtractionResult) -> None:
        """Extract metadata from Nextflow files."""
        if result.file_type == "nextflow_config":
            self._extract_nextflow_config_metadata(file_path, result)
        elif result.file_type == "nextflow_workflow":
            self._extract_nextflow_workflow_metadata(file_path, result)
        elif result.file_type == "nextflow_params":
            self._extract_nextflow_params_metadata(file_path, result)

    def _extract_nextflow_config_metadata(self, file_path: Path, result: MetadataExtractionResult) -> None:
        """Extract metadata from Nextflow config files."""
        try:
            content = file_path.read_text(encoding='utf-8')

            # Extract basic configuration information
            result.metadata_fields.extend([
                MetadataField("file_size", len(content), "int", "file_info", "Configuration file size"),
                MetadataField("line_count", len(content.splitlines()), "int", "file_info", "Number of lines")
            ])

            # Extract process configurations
            process_matches = re.findall(r'process\s+(\w+)\s*{([^}]+)}', content, re.MULTILINE | re.DOTALL)
            if process_matches:
                process_names = [match[0] for match in process_matches]
                result.metadata_fields.append(
                    MetadataField("processes", process_names, "list", "workflow",
                                "Defined processes", "high")
                )

            # Extract profiles
            profile_matches = re.findall(r'profiles\s*{([^}]+)}', content, re.MULTILINE | re.DOTALL)
            if profile_matches:
                # Extract profile names from the profiles block
                profile_content = profile_matches[0]
                profile_names = re.findall(r'(\w+)\s*{', profile_content)
                result.metadata_fields.append(
                    MetadataField("profiles", profile_names, "list", "configuration",
                                "Available execution profiles", "medium")
                )

            # Extract resource requirements
            memory_matches = re.findall(r'memory\s*=\s*[\'"]([^\'"]+)[\'"]', content)
            cpu_matches = re.findall(r'cpus\s*=\s*(\d+)', content)

            if memory_matches:
                result.metadata_fields.append(
                    MetadataField("memory_requirements", memory_matches, "list", "resources",
                                "Memory requirements found", "medium")
                )

            if cpu_matches:
                result.metadata_fields.append(
                    MetadataField("cpu_requirements", [int(x) for x in cpu_matches], "list", "resources",
                                "CPU requirements found", "medium")
                )

            # Extract container configurations
            container_matches = re.findall(r'container\s*=\s*[\'"]([^\'"]+)[\'"]', content)
            if container_matches:
                result.metadata_fields.append(
                    MetadataField("containers", container_matches, "list", "containers",
                                "Container images specified", "high")
                )

        except Exception as e:
            result.issues.append(f"Failed to parse Nextflow config: {e}")

    def _extract_nextflow_workflow_metadata(self, file_path: Path, result: MetadataExtractionResult) -> None:
        """Extract metadata from Nextflow workflow files."""
        try:
            content = file_path.read_text(encoding='utf-8')

            # Basic file information
            result.metadata_fields.extend([
                MetadataField("file_size", len(content), "int", "file_info", "Workflow file size"),
                MetadataField("line_count", len(content.splitlines()), "int", "file_info", "Number of lines")
            ])

            # Extract workflow definition
            workflow_match = re.search(r'workflow\s*{([^}]+)}', content, re.MULTILINE | re.DOTALL)
            if workflow_match:
                result.metadata_fields.append(
                    MetadataField("has_workflow_block", True, "bool", "structure",
                                "Contains workflow definition", "critical")
                )

            # Extract process calls
            process_calls = re.findall(r'(\w+)\s*\(', content)
            if process_calls:
                result.metadata_fields.append(
                    MetadataField("process_calls", list(set(process_calls)), "list", "workflow",
                                "Process calls in workflow", "high")
                )

            # Extract channel operations
            channel_ops = re.findall(r'\.(\w+)\s*\(', content)
            if channel_ops:
                common_ops = ['map', 'filter', 'collect', 'flatten', 'groupTuple', 'join']
                found_ops = [op for op in set(channel_ops) if op in common_ops]
                if found_ops:
                    result.metadata_fields.append(
                        MetadataField("channel_operations", found_ops, "list", "workflow",
                                    "Channel operations used", "medium")
                    )

            # Extract input/output declarations
            input_matches = re.findall(r'input:\s*([^\n]+)', content)
            output_matches = re.findall(r'output:\s*([^\n]+)', content)

            if input_matches:
                result.metadata_fields.append(
                    MetadataField("input_declarations", len(input_matches), "int", "io",
                                "Number of input declarations", "high")
                )

            if output_matches:
                result.metadata_fields.append(
                    MetadataField("output_declarations", len(output_matches), "int", "io",
                                "Number of output declarations", "high")
                )

        except Exception as e:
            result.issues.append(f"Failed to parse Nextflow workflow: {e}")

    def _extract_nextflow_params_metadata(self, file_path: Path, result: MetadataExtractionResult) -> None:
        """Extract metadata from Nextflow parameter files."""
        try:
            if file_path.suffix.lower() in ['.yaml', '.yml']:
                with open(file_path, 'r', encoding='utf-8') as f:
                    params = yaml.safe_load(f)
            else:
                with open(file_path, 'r', encoding='utf-8') as f:
                    params = json.load(f)

            if isinstance(params, dict):
                result.metadata_fields.append(
                    MetadataField("parameter_count", len(params), "int", "parameters",
                                "Number of parameters", "medium")
                )

                # Categorize parameters
                input_params = [k for k in params.keys() if 'input' in k.lower() or 'file' in k.lower()]
                output_params = [k for k in params.keys() if 'output' in k.lower() or 'result' in k.lower()]

                if input_params:
                    result.metadata_fields.append(
                        MetadataField("input_parameters", input_params, "list", "parameters",
                                    "Input-related parameters", "high")
                    )

                if output_params:
                    result.metadata_fields.append(
                        MetadataField("output_parameters", output_params, "list", "parameters",
                                    "Output-related parameters", "high")
                    )

        except Exception as e:
            result.issues.append(f"Failed to parse parameter file: {e}")

    def _extract_viash_metadata(self, file_path: Path, result: MetadataExtractionResult) -> None:
        """Extract metadata from Viash files."""
        if result.file_type == "viash_config":
            self._extract_viash_config_metadata(file_path, result)
        elif result.file_type == "viash_script":
            self._extract_viash_script_metadata(file_path, result)

    def _extract_viash_config_metadata(self, file_path: Path, result: MetadataExtractionResult) -> None:
        """Extract metadata from Viash config files."""
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                config = yaml.safe_load(f)

            if not isinstance(config, dict):
                result.issues.append("Viash config is not a valid YAML dictionary")
                return

            # Extract basic component information
            if 'name' in config:
                result.metadata_fields.append(
                    MetadataField("component_name", config['name'], "str", "component",
                                "Component name", "critical")
                )

            if 'description' in config:
                result.metadata_fields.append(
                    MetadataField("description", config['description'], "str", "component",
                                "Component description", "high")
                )

            if 'version' in config:
                result.metadata_fields.append(
                    MetadataField("version", config['version'], "str", "component",
                                "Component version", "medium")
                )

            # Extract functionality information
            if 'functionality' in config:
                func = config['functionality']

                if 'arguments' in func:
                    args = func['arguments']
                    result.metadata_fields.append(
                        MetadataField("argument_count", len(args), "int", "interface",
                                    "Number of arguments", "high")
                    )

                    # Categorize arguments
                    input_args = [arg for arg in args if arg.get('direction') == 'input']
                    output_args = [arg for arg in args if arg.get('direction') == 'output']

                    result.metadata_fields.extend([
                        MetadataField("input_arguments", len(input_args), "int", "interface",
                                    "Number of input arguments", "high"),
                        MetadataField("output_arguments", len(output_args), "int", "interface",
                                    "Number of output arguments", "high")
                    ])

                if 'resources' in func:
                    resources = func['resources']
                    if isinstance(resources, list):
                        script_resources = [r for r in resources if r.get('type') == 'bash_script']
                        if script_resources:
                            result.metadata_fields.append(
                                MetadataField("script_resources", len(script_resources), "int", "resources",
                                            "Number of script resources", "medium")
                            )

            # Extract platform information
            if 'platforms' in config:
                platforms = config['platforms']
                if isinstance(platforms, list):
                    platform_types = [p.get('type', 'unknown') for p in platforms]
                    result.metadata_fields.append(
                        MetadataField("platform_types", platform_types, "list", "platforms",
                                    "Supported platform types", "high")
                    )

                    # Check for Docker platforms
                    docker_platforms = [p for p in platforms if p.get('type') == 'docker']
                    if docker_platforms:
                        images = [p.get('image', 'unknown') for p in docker_platforms]
                        result.metadata_fields.append(
                            MetadataField("docker_images", images, "list", "containers",
                                        "Docker images used", "medium")
                        )

        except Exception as e:
            result.issues.append(f"Failed to parse Viash config: {e}")

    def _extract_viash_script_metadata(self, file_path: Path, result: MetadataExtractionResult) -> None:
        """Extract metadata from Viash script files."""
        try:
            content = file_path.read_text(encoding='utf-8')

            # Basic file information
            result.metadata_fields.extend([
                MetadataField("file_size", len(content), "int", "file_info", "Script file size"),
                MetadataField("line_count", len(content.splitlines()), "int", "file_info", "Number of lines"),
                MetadataField("script_language", file_path.suffix[1:], "str", "script", "Script language", "medium")
            ])

            # Language-specific analysis
            if file_path.suffix.lower() == '.py':
                self._analyze_python_script(content, result)
            elif file_path.suffix.lower() == '.r':
                self._analyze_r_script(content, result)
            elif file_path.suffix.lower() == '.sh':
                self._analyze_bash_script(content, result)

        except Exception as e:
            result.issues.append(f"Failed to analyze script: {e}")

    def _analyze_python_script(self, content: str, result: MetadataExtractionResult) -> None:
        """Analyze Python script content."""
        # Extract imports
        import_matches = re.findall(r'^(?:from\s+(\S+)\s+)?import\s+([^\n]+)', content, re.MULTILINE)
        imports = []
        for match in import_matches:
            if match[0]:  # from X import Y
                imports.append(match[0])
            else:  # import X
                imports.extend([imp.strip() for imp in match[1].split(',')])

        if imports:
            # Focus on scientific/bioinformatics libraries
            bio_libs = [imp for imp in imports if any(lib in imp.lower() for lib in
                       ['pandas', 'numpy', 'scipy', 'sklearn', 'anndata', 'scanpy', 'spatialdata'])]
            if bio_libs:
                result.metadata_fields.append(
                    MetadataField("bioinformatics_libraries", bio_libs, "list", "dependencies",
                                "Bioinformatics libraries used", "high")
                )

        # Extract function definitions
        func_matches = re.findall(r'^def\s+(\w+)\s*\(', content, re.MULTILINE)
        if func_matches:
            result.metadata_fields.append(
                MetadataField("function_count", len(func_matches), "int", "code_structure",
                            "Number of functions defined", "medium")
            )

    def _analyze_r_script(self, content: str, result: MetadataExtractionResult) -> None:
        """Analyze R script content."""
        # Extract library calls
        library_matches = re.findall(r'library\s*\(\s*([^)]+)\s*\)', content)
        if library_matches:
            bio_libs = [lib.strip('\'"') for lib in library_matches
                       if any(bio in lib.lower() for bio in ['seurat', 'bioconductor', 'spatstat'])]
            if bio_libs:
                result.metadata_fields.append(
                    MetadataField("r_bioinformatics_libraries", bio_libs, "list", "dependencies",
                                "R bioinformatics libraries used", "high")
                )

    def _analyze_bash_script(self, content: str, result: MetadataExtractionResult) -> None:
        """Analyze bash script content."""
        # Extract command usage
        commands = re.findall(r'^(\w+)\s+', content, re.MULTILINE)
        bio_commands = [cmd for cmd in set(commands)
                       if cmd in ['samtools', 'bcftools', 'bedtools', 'fastqc', 'bwa', 'bowtie2']]
        if bio_commands:
            result.metadata_fields.append(
                MetadataField("bioinformatics_commands", bio_commands, "list", "tools",
                            "Bioinformatics tools used", "high")
            )

    def _extract_spatial_data_metadata(self, file_path: Path, result: MetadataExtractionResult) -> None:
        """Extract metadata from spatial data files."""
        # This leverages the spatial validation module
        try:
            from .spatial_validation import SpatialDataValidator, ValidationLevel

            validator = SpatialDataValidator()
            validation_result = validator.validate_file(file_path, ValidationLevel.DOMAIN)

            # Convert validation metadata to metadata fields
            for key, value in validation_result.metadata.items():
                category = "spatial_data"
                importance = "medium"

                if key in ['n_obs', 'n_vars', 'expression_matrix_shape']:
                    importance = "high"
                elif key in ['components', 'spatial_keys']:
                    importance = "critical"

                result.metadata_fields.append(
                    MetadataField(key, value, type(value).__name__, category,
                                f"Spatial data: {key}", importance)
                )

            # Add validation issues as quality metrics
            if validation_result.issues:
                result.quality_metrics['validation_issues'] = len(validation_result.issues)
                result.quality_metrics['has_errors'] = validation_result.has_errors
                result.quality_metrics['has_warnings'] = validation_result.has_warnings

        except ImportError:
            result.issues.append("Spatial validation module not available")
        except Exception as e:
            result.issues.append(f"Failed to extract spatial metadata: {e}")

    def _extract_generic_metadata(self, file_path: Path, result: MetadataExtractionResult) -> None:
        """Extract generic metadata from unknown file types."""
        try:
            # Basic file information
            stat = file_path.stat()
            result.metadata_fields.extend([
                MetadataField("file_size", stat.st_size, "int", "file_info", "File size in bytes"),
                MetadataField("file_extension", file_path.suffix, "str", "file_info", "File extension"),
                MetadataField("modification_time", stat.st_mtime, "float", "file_info", "Last modification time")
            ])

            # Try to determine if it's a text file
            try:
                content = file_path.read_text(encoding='utf-8')
                result.metadata_fields.append(
                    MetadataField("line_count", len(content.splitlines()), "int", "file_info", "Number of lines")
                )

                # Check for common bioinformatics patterns
                if re.search(r'(gene|cell|expression|spatial|coordinate)', content, re.IGNORECASE):
                    result.metadata_fields.append(
                        MetadataField("contains_bio_keywords", True, "bool", "content",
                                    "Contains bioinformatics keywords", "medium")
                    )

            except UnicodeDecodeError:
                result.metadata_fields.append(
                    MetadataField("is_binary", True, "bool", "file_info", "File appears to be binary")
                )

        except Exception as e:
            result.issues.append(f"Failed to extract generic metadata: {e}")

    def _assess_metadata_quality(self, result: MetadataExtractionResult) -> None:
        """Assess the quality of extracted metadata."""
        # Count metadata by importance
        critical_count = len([f for f in result.metadata_fields if f.importance == "critical"])
        high_count = len([f for f in result.metadata_fields if f.importance == "high"])
        total_count = len(result.metadata_fields)

        result.quality_metrics.update({
            'total_metadata_fields': total_count,
            'critical_fields': critical_count,
            'high_importance_fields': high_count,
            'metadata_completeness': (critical_count + high_count) / max(total_count, 1)
        })

        # Generate suggestions based on file type and metadata
        if result.file_type.startswith('nextflow'):
            if not any(f.name == 'processes' for f in result.metadata_fields):
                result.suggestions.append("Consider adding process definitions to improve workflow structure")

            if not any(f.name == 'containers' for f in result.metadata_fields):
                result.suggestions.append("Consider adding container specifications for reproducibility")

        elif result.file_type.startswith('viash'):
            if not any(f.name == 'component_name' for f in result.metadata_fields):
                result.suggestions.append("Component should have a descriptive name")

            if not any(f.name == 'docker_images' for f in result.metadata_fields):
                result.suggestions.append("Consider adding Docker platform for containerized execution")

        elif result.file_type.startswith('spatial'):
            if not any(f.name == 'spatial_keys' for f in result.metadata_fields):
                result.suggestions.append("Spatial data should include coordinate information")

        # General suggestions
        if total_count < 5:
            result.suggestions.append("Limited metadata extracted - file may need manual review")

        if result.quality_metrics.get('metadata_completeness', 0) < 0.5:
            result.suggestions.append("Consider adding more descriptive metadata to improve analysis")

    def analyze_workflow_dependencies(self, file_paths: List[Union[str, Path]]) -> Dict[str, Any]:
        """Analyze dependencies across multiple workflow files.

        Args:
            file_paths: List of workflow file paths to analyze

        Returns:
            Dictionary with dependency analysis results
        """
        results = []
        for file_path in file_paths:
            result = self.extract_metadata(file_path)
            results.append(result)

        # Aggregate dependency information
        all_containers = []
        all_processes = []
        all_libraries = []

        for result in results:
            # Extract containers
            container_fields = [f for f in result.metadata_fields if 'container' in f.name.lower()]
            for field in container_fields:
                if isinstance(field.value, list):
                    all_containers.extend(field.value)
                else:
                    all_containers.append(field.value)

            # Extract processes
            process_fields = [f for f in result.metadata_fields if 'process' in f.name.lower()]
            for field in process_fields:
                if isinstance(field.value, list):
                    all_processes.extend(field.value)

            # Extract libraries
            lib_fields = [f for f in result.metadata_fields if 'librar' in f.name.lower()]
            for field in lib_fields:
                if isinstance(field.value, list):
                    all_libraries.extend(field.value)

        return {
            'total_files_analyzed': len(results),
            'successful_extractions': len([r for r in results if r.extraction_success]),
            'unique_containers': list(set(all_containers)),
            'unique_processes': list(set(all_processes)),
            'unique_libraries': list(set(all_libraries)),
            'container_count': len(set(all_containers)),
            'process_count': len(set(all_processes)),
            'library_count': len(set(all_libraries)),
            'file_types': list(set(r.file_type for r in results)),
            'total_issues': sum(len(r.issues) for r in results)
        }
