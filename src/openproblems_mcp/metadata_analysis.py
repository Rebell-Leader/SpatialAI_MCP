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

    def extract_spatial_transcriptomics_metadata(self, file_path: Union[str, Path]) -> MetadataExtractionResult:
        """Extract metadata specifically for spatial transcriptomics workflows and data.

        Args:
            file_path: Path to the spatial transcriptomics file

        Returns:
            MetadataExtractionResult with spatial-specific metadata
        """
        result = self.extract_metadata(file_path)

        # Enhance with spatial transcriptomics specific analysis
        if result.extraction_success:
            self._enhance_spatial_metadata(result)
            self._analyze_spatial_workflow_patterns(result)
            self._assess_spatial_data_quality(result)

        return result

    def _enhance_spatial_metadata(self, result: MetadataExtractionResult) -> None:
        """Enhance metadata with spatial transcriptomics specific information."""
        file_path = Path(result.file_path)

        try:
            # Check for spatial transcriptomics specific patterns
            if result.file_type.startswith('spatial'):
                # Add spatial data specific metadata
                self._extract_spatial_data_characteristics(file_path, result)

            elif result.file_type.startswith('nextflow') or result.file_type.startswith('viash'):
                # Check for spatial transcriptomics workflow patterns
                self._identify_spatial_workflow_components(file_path, result)

            # Check for spatial transcriptomics libraries and tools
            self._identify_spatial_tools_and_libraries(result)

        except Exception as e:
            result.issues.append(f"Failed to enhance spatial metadata: {e}")

    def _extract_spatial_data_characteristics(self, file_path: Path, result: MetadataExtractionResult) -> None:
        """Extract characteristics specific to spatial data files."""
        try:
            # Use spatial validation to get detailed characteristics
            from .spatial_validation import SpatialDataValidator, ValidationLevel

            validator = SpatialDataValidator()
            validation_result = validator.validate_file(file_path, ValidationLevel.DOMAIN)

            # Add spatial characteristics as metadata fields
            if validation_result.metadata:
                for key, value in validation_result.metadata.items():
                    importance = "high" if key in ['n_obs', 'n_vars', 'spatial_keys', 'components'] else "medium"

                    result.metadata_fields.append(
                        MetadataField(
                            name=f"spatial_{key}",
                            value=value,
                            data_type=type(value).__name__,
                            category="spatial_characteristics",
                            description=f"Spatial data characteristic: {key}",
                            importance=importance
                        )
                    )

            # Add data quality assessment
            quality_assessment = validator.assess_spatial_data_quality(file_path)
            result.quality_metrics.update({
                'spatial_quality_score': quality_assessment.get('overall_quality', 'unknown'),
                'spatial_issues_count': quality_assessment.get('issues_found', 0),
                'spatial_file_size_mb': quality_assessment.get('file_size_mb', 0)
            })

        except ImportError:
            result.issues.append("Spatial validation module not available for enhanced analysis")
        except Exception as e:
            result.issues.append(f"Failed to extract spatial data characteristics: {e}")

    def _identify_spatial_workflow_components(self, file_path: Path, result: MetadataExtractionResult) -> None:
        """Identify spatial transcriptomics specific workflow components."""
        try:
            content = file_path.read_text(encoding='utf-8')

            # Spatial transcriptomics specific patterns
            spatial_patterns = {
                'spatial_libraries': [
                    'spatialdata', 'squidpy', 'scanpy', 'ata', 'spatstat',
                    'seurat', 'giotto', 'stlearn', 'stereoscope'
                ],
                'spatial_processes': [
                    'spatial_clustering', 'spatial_autocorrelation', 'spatial_deconvolution',
                    'spatial_registration', 'spatial_segmentation', 'spatial_visualization'
                ],
                'spatial_file_formats': [
                    '.h5ad', '.zarr', '.spatialdata', '.loom', '.h5seurat'
                ],
                'spatial_analysis_types': [
                    'visium', '10x_genomics', 'slide_seq', 'merfish', 'seqfish',
                    'osmfish', 'starmap', 'hdst', 'dbit_seq'
                ]
            }

            content_lower = content.lower()

            # Check for spatial libraries
            found_libraries = [lib for lib in spatial_patterns['spatial_libraries']
                             if lib in content_lower]
            if found_libraries:
                result.metadata_fields.append(
                    MetadataField(
                        name="spatial_libraries_detected",
                        value=found_libraries,
                        data_type="list",
                        category="spatial_workflow",
                        description="Spatial transcriptomics libraries found in workflow",
                        importance="high"
                    )
                )

            # Check for spatial processes
            found_processes = [proc for proc in spatial_patterns['spatial_processes']
                             if proc.replace('_', ' ') in content_lower or proc in content_lower]
            if found_processes:
                result.metadata_fields.append(
                    MetadataField(
                        name="spatial_processes_detected",
                        value=found_processes,
                        data_type="list",
                        category="spatial_workflow",
                        description="Spatial analysis processes found in workflow",
                        importance="high"
                    )
                )

            # Check for spatial analysis types
            found_analysis_types = [atype for atype in spatial_patterns['spatial_analysis_types']
                                  if atype in content_lower]
            if found_analysis_types:
                result.metadata_fields.append(
                    MetadataField(
                        name="spatial_analysis_types",
                        value=found_analysis_types,
                        data_type="list",
                        category="spatial_workflow",
                        description="Spatial transcriptomics platforms/methods detected",
                        importance="high"
                    )
                )

            # Check for file format handling
            found_formats = [fmt for fmt in spatial_patterns['spatial_file_formats']
                           if fmt in content_lower]
            if found_formats:
                result.metadata_fields.append(
                    MetadataField(
                        name="spatial_file_formats",
                        value=found_formats,
                        data_type="list",
                        category="spatial_workflow",
                        description="Spatial data formats handled by workflow",
                        importance="medium"
                    )
                )

        except Exception as e:
            result.issues.append(f"Failed to identify spatial workflow components: {e}")

    def _identify_spatial_tools_and_libraries(self, result: MetadataExtractionResult) -> None:
        """Identify spatial transcriptomics tools and libraries from existing metadata."""
        try:
            # Check existing metadata fields for spatial tools
            spatial_tools = set()
            spatial_libraries = set()

            for field in result.metadata_fields:
                if field.category in ['dependencies', 'tools', 'libraries']:
                    if isinstance(field.value, list):
                        values = field.value
                    else:
                        values = [field.value]

                    for value in values:
                        value_str = str(value).lower()

                        # Check for spatial transcriptomics tools
                        spatial_tool_patterns = [
                            'spatial', 'squidpy', 'scanpy', 'anndata', 'spatialdata',
                            'seurat', 'giotto', 'stlearn', 'stereoscope', 'cell2location',
                            'tangram', 'novosparc', 'paste', 'harmony'
                        ]

                        for pattern in spatial_tool_patterns:
                            if pattern in value_str:
                                if 'library' in field.name.lower() or 'import' in field.name.lower():
                                    spatial_libraries.add(value)
                                else:
                                    spatial_tools.add(value)

            # Add spatial tools metadata
            if spatial_tools:
                result.metadata_fields.append(
                    MetadataField(
                        name="spatial_tools_identified",
                        value=list(spatial_tools),
                        data_type="list",
                        category="spatial_tools",
                        description="Spatial transcriptomics tools identified",
                        importance="high"
                    )
                )

            if spatial_libraries:
                result.metadata_fields.append(
                    MetadataField(
                        name="spatial_libraries_identified",
                        value=list(spatial_libraries),
                        data_type="list",
tegory="spatial_tools",
                        description="Spatial transcriptomics libraries identified",
                        importance="high"
                    )
                )

        except Exception as e:
            result.issues.append(f"Failed to identify spatial tools and libraries: {e}")

    def _analyze_spatial_workflow_patterns(self, result: MetadataExtractionResult) -> None:
        """Analyze patterns specific to spatial transcriptomics workflows."""
        try:
            # Count spatial-related metadata fields
            spatial_fields = [f for f in result.metadata_fields if 'spatial' in f.name.lower()]

            if spatial_fields:
                result.quality_metrics['spatial_metadata_richness'] = len(spatial_fields)

                # Categorize spatial workflow complexity
                if len(spatial_fields) >= 5:
                    complexity = "high"
                elif len(spatial_fields) >= 3:
                    complexity = "medium"
                else:
                    complexity = "low"

                result.quality_metrics['spatial_workflow_complexity'] = complexity

                # Add suggestions based on complexity
                if complexity == "low":
                    result.suggestions.append(
                        "Consider adding more spatial-specific metadata for better workflow documentation"
                    )
                elif complexity == "high":
                    result.suggestions.append(
                        "Rich spatial metadata detected - ensure proper documentation of spatial parameters"
                    )

            # Check for common spatial transcriptomics workflow patterns
            workflow_patterns = {
                'preprocessing': ['normalization', 'filtering', 'quality_control'],
                'spatial_analysis': ['clustering', 'neighborhood', 'autocorrelation'],
                'visualization': ['plotting', 'embedding', 'dimensionality_reduction'],
                'integration': ['batch_correction', 'alignment', 'registration']
            }

            detected_patterns = {}
            for pattern_type, keywords in workflow_patterns.items():
                found_keywords = []
                for field in result.metadata_fields:
                    field_text = f"{field.name} {field.value}".lower()
                    for keyword in keywords:
                        if keyword in field_text:
                            found_keywords.append(keyword)

                if found_keywords:
                    detected_patterns[pattern_type] = list(set(found_keywords))

            if detected_patterns:
                result.metadata_fields.append(
                    MetadataField(
                        name="spatial_workflow_patterns",
                        value=detected_patterns,
                        data_type="dict",
                        category="spatial_workflow",
                        description="Spatial transcriptomics workflow patterns detected",
                        importance="medium"
                    )
                )

        except Exception as e:
            result.issues.append(f"Failed to analyze spatial workflow patterns: {e}")

    def _assess_spatial_data_quality(self, result: MetadataExtractionResult) -> None:
        """Assess data quality specific to spatial transcriptomics."""
        try:
            # Check for quality indicators in metadata
            quality_indicators = {
                'has_spatial_coordinates': False,
                'has_expression_data': False,
                'has_quality_metrics': False,
                'has_spatial_images': False,
                'has_proper_annotations': False
            }

            for field in result.metadata_fields:
                field_name_lower = field.name.lower()

                if 'spatial' in field_name_lower and 'coord' in field_name_lower:
                    quality_indicators['has_spatial_coordinates'] = True
                elif 'expression' in field_name_lower or 'matrix' in field_name_lower:
                    quality_indicators['has_expression_data'] = True
                elif 'quality' in field_name_lower or 'qc' in field_name_lower:
                    quality_indicators['has_quality_metrics'] = True
                elif 'image' in field_name_lower:
                    quality_indicators['has_spatial_images'] = True
                elif 'annotation' in field_name_lower or 'metadata' in field_name_lower:
                    quality_indicators['has_proper_annotations'] = True

            # Calculate quality score
            quality_score = sum(quality_indicators.values()) / len(quality_indicators)
            result.quality_metrics['spatial_data_quality_score'] = quality_score

            # Add quality-based suggestions
            if quality_score < 0.5:
                result.suggestions.append(
                    "Spatial data appears to be missing key components - verify completeness"
                )
            elif quality_score >= 0.8:
                result.suggestions.append(
                    "Spatial data appears comprehensive and well-structured"
                )

            # Specific quality checks
            if not quality_indicators['has_spatial_coordinates']:
                result.suggestions.append(
                    "Add spatial coordinate information for spatial analysis"
                )

            if not quality_indicators['has_expression_data']:
                result.suggestions.append(
                    "Ensure expression data is included for transcriptomics analysis"
                )

        except Exception as e:
            result.issues.append(f"Failed to assess spatial data quality: {e}")

    def analyze_workflow_configuration_advanced(
        self,
        file_path: Union[str, Path],
        include_spatial_analysis: bool = True,
        include_dependency_graph: bool = True
    ) -> Dict[str, Any]:
        """Perform advanced workflow configuration analysis.

        Args:
            file_path: Path to the workflow configuration file
            include_spatial_analysis: Whether to include spatial transcriptomics analysis
            include_dependency_graph: Whether to build dependency graph

        Returns:
            Advanced workflow analysis results
        """
        result = self.extract_metadata(file_path)

        analysis = {
            'file_path': str(file_path),
            'file_type': result.file_type,
            'extraction_success': result.extraction_success,
            'workflow_complexity': 'unknown',
            'spatial_compatibility': False,
            'dependency_graph': {},
            'optimization_suggestions': [],
            'compatibility_issues': []
        }

        if not result.extraction_success:
            analysis['compatibility_issues'].extend(result.issues)
            return analysis

        try:
            # Assess workflow complexity
            complexity_indicators = {
                'process_count': 0,
                'parameter_count': 0,
                'dependency_count': 0,
                'container_count': 0
            }

            for field in result.metadata_fields:
                if 'process' in field.name.lower():
                    if isinstance(field.value, list):
                        complexity_indicators['process_count'] += len(field.value)
                    elif isinstance(field.value, int):
                        complexity_indicators['process_count'] += field.value
                elif 'parameter' in field.name.lower() or 'argument' in field.name.lower():
                    if isinstance(field.value, int):
                        complexity_indicators['parameter_count'] += field.value
                elif 'container' in field.name.lower() or 'docker' in field.name.lower():
                    if isinstance(field.value, list):
                        complexity_indicators['container_count'] += len(field.value)
                elif 'librar' in field.name.lower() or 'depend' in field.name.lower():
                    if isinstance(field.value, list):
                        complexity_indicators['dependency_count'] += len(field.value)

            # Calculate complexity score
            complexity_score = (
                complexity_indicators['process_count'] * 2 +
                complexity_indicators['parameter_count'] +
                complexity_indicators['dependency_count'] +
                complexity_indicators['container_count']
            )

            if complexity_score >= 20:
                analysis['workflow_complexity'] = 'high'
            elif complexity_score >= 10:
                analysis['workflow_complexity'] = 'medium'
            else:
                analysis['workflow_complexity'] = 'low'

            # Spatial compatibility analysis
            if include_spatial_analysis:
                spatial_fields = [f for f in result.metadata_fields if 'spatial' in f.name.lower()]
                spatial_libraries = [f for f in result.metadata_fields
                                   if any(lib in str(f.value).lower() for lib in
                                         ['spatialdata', 'squidpy', 'scanpy', 'seurat'])]

                analysis['spatial_compatibility'] = len(spatial_fields) > 0 or len(spatial_libraries) > 0

                if analysis['spatial_compatibility']:
                    analysis['optimization_suggestions'].append(
                        "Workflow appears compatible with spatial transcriptomics data"
                    )
                else:
                    analysis['optimization_suggestions'].append(
                        "Consider adding spatial transcriptomics libraries for spatial analysis"
                    )

            # Dependency graph analysis
            if include_dependency_graph:
                containers = []
                libraries = []
                processes = []

                for field in result.metadata_fields:
                    if 'container' in field.name.lower():
                        if isinstance(field.value, list):
                            containers.extend(field.value)
                        else:
                            containers.append(field.value)
                    elif 'librar' in field.name.lower():
                        if isinstance(field.value, list):
                            libraries.extend(field.value)
                        else:
                            libraries.append(field.value)
                    elif 'process' in field.name.lower():
                        if isinstance(field.value, list):
                            processes.extend(field.value)
                        else:
                            processes.append(field.value)

                analysis['dependency_graph'] = {
                    'containers': list(set(containers)),
                    'libraries': list(set(libraries)),
                    'processes': list(set(processes))
                }

            # Generate optimization suggestions
            if complexity_indicators['container_count'] == 0:
                analysis['optimization_suggestions'].append(
                    "Consider adding container specifications for reproducibility"
                )

            if complexity_indicators['process_count'] > 10:
                analysis['optimization_suggestions'].append(
                    "High process count detected - consider workflow modularization"
                )

            if complexity_indicators['parameter_count'] > 20:
                analysis['optimization_suggestions'].append(
                    "Many parameters detected - consider parameter grouping or defaults"
                )

        except Exception as e:
            analysis['compatibility_issues'].append(f"Advanced analysis failed: {e}")

        return analysis
