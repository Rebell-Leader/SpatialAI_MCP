"""MCP tools for spatial transcriptomics data validation and analysis."""

import logging
import json
from pathlib import Path
from typing import Dict, List, Any, Optional, Union

from .spatial_validation import (
    SpatialDataValidator,
    ValidationLevel,
    ValidationResult,
    ValidationStatus
)
from .metadata_analysis import BioinformaticsMetadataExtractor, MetadataExtractionResult

logger = logging.getLogger(__name__)


class SpatialMCPTools:
    """MCP tools for spatial transcriptomics data operations."""

    def __init__(self):
        """Initialize spatial MCP tools."""
        self.validator = SpatialDataValidator()
        self.metadata_extractor = BioinformaticsMetadataExtractor()

    def validate_spatial_data(
        self,
        file_path: str,
        validation_level: str = "structure",
        return_format: str = "summary"
    ) -> str:
        """Validate spatial transcriptomics data file.

        Args:
            file_path: Path to the spatial data file to validate
            validation_level: Level of validation (basic, structure, integrity, domain)
            return_format: Format of return data (ary, detailed, json)

        Returns:
            Formatted validation results
        """
        try:
            # Parse validation level
            try:
                level = ValidationLevel(validation_level.lower())
            except ValueError:
                return f"❌ Invalid validation level: {validation_level}. Use: basic, structure, integrity, domain"

            # Validate the file
            result = self.validator.validate_file(file_path, level)

            # Format response based on requested format
            if return_format == "json":
                return self._format_result_json(result)
            elif return_format == "detailed":
                return self._format_result_detailed(result)
            else:  # summary
                return self._format_result_summary(result)

        except Exception as e:
            logger.error(f"Spatial data validation failed: {e}")
            return f"❌ Validation failed: {e}"

    def validate_multiple_spatial_files(
        self,
        file_paths: str,  # JSON string of file paths
        validation_level: str = "structure",
        return_format: str = "summary"
    ) -> str:
        """Validate multiple spatial transcriptomics data files.

        Args:
            file_paths: JSON string containing list of file paths
            validation_level: Level of validation (basic, structure, integrity, domain)
            return_format: Format of return data (summary, detailed, json)

        Returns:
            Formatted validation results for all files
        """
        try:
            # Parse file paths
            try:
                paths = json.loads(file_paths)
                if not isinstance(paths, list):
                    return "❌ file_paths must be a JSON list of file paths"
            except json.JSONDecodeError:
                return "❌ Invalid JSON format for file_paths"

            # Parse validation level
            try:
                level = ValidationLevel(validation_level.lower())
            except ValueError:
                return f"❌ Invalid validation level: {validation_level}. Use: basic, structure, integrity, domain"

            # Validate all files
            results = self.validator.validate_multiple_files(paths, level)

            # Generate summary
            summary = self.validator.get_validation_summary(results)

            # Format response
            if return_format == "json":
                return json.dumps({
                    'summary': summary,
                    'results': [self._result_to_dict(r) for r in results]
                }, indent=2)
            elif return_format == "detailed":
                return self._format_multiple_results_detailed(results, summary)
            else:  # summary
                return self._format_multiple_results_summary(results, summary)

        except Exception as e:
            logger.error(f"Multiple file validation failed: {e}")
            return f"❌ Multiple file validation failed: {e}"

    def analyze_spatial_metadata(
        self,
        file_path: str,
        extract_coordinates: bool = True,
        extract_gene_info: bool = True
    ) -> str:
        """Extract and analyze metadata from spatial transcriptomics data.

        Args:
            file_path: Path to the spatial data file
            extract_coordinates: Whether to extract spatial coordinate information
            extract_gene_info: Whether to extract gene/feature information

        Returns:
            Formatted metadata analysis
        """
        try:
            # First validate to get basic metadata
            result = self.validator.validate_file(file_path, ValidationLevel.DOMAIN)

            if not result.is_valid and result.has_errors:
                return f"❌ Cannot analyze metadata - file validation failed:\n{self._format_issues(result.issues)}"

            # Format metadata analysis
            analysis = f"📊 Spatial Data Metadata Analysis\n"
            analysis += f"File: {file_path}\n"
            analysis += f"Format: {result.format_type}\n"
            analysis += f"Size: {self._format_file_size(result.file_size)}\n\n"

            # Add format-specific metadata
            if result.metadata:
                analysis += "**Metadata Details:**\n"
                for key, value in result.metadata.items():
                    if key == 'components' and isinstance(value, list):
                        analysis += f"  • Components: {', '.join(value)}\n"
                    elif key == 'expression_matrix_shape':
                        analysis += f"  • Expression Matrix Shape: {value} (cells × genes)\n"
                    elif key == 'n_obs':
                        analysis += f"  • Number of Observations: {value}\n"
                    elif key == 'n_vars':
                        analysis += f"  • Number of Variables: {value}\n"
                    elif key == 'spatial_keys':
                        analysis += f"  • Spatial Coordinate Keys: {', '.join(value)}\n"
                    elif key == 'obs_columns':
                        analysis += f"  • Observation Columns: {len(value)} columns\n"
                    elif key == 'var_columns':
                        analysis += f"  • Variable Columns: {len(value)} columns\n"
                    else:
                        analysis += f"  • {key}: {value}\n"

            # Add validation issues as metadata insights
            if result.issues:
                analysis += "\n**Data Quality Insights:**\n"
                for issue in result.issues:
                    if issue.level == ValidationStatus.WARNING:
                        analysis += f"  ⚠️ {issue.message}\n"
                    elif issue.level in [ValidationStatus.ERROR, ValidationStatus.CRITICAL]:
                        analysis += f"  ❌ {issue.message}\n"

            return analysis

        except Exception as e:
            logger.error(f"Metadata analysis failed: {e}")
            return f"❌ Metadata analysis failed: {e}"

    def check_spatial_data_compatibility(
        self,
        file_paths: str,  # JSON string of file paths
        check_coordinates: bool = True,
        check_gene_overlap: bool = True
    ) -> str:
        """Check compatibility between multiple spatial data files.

        Args:
            file_paths: JSON string containing list of file paths to compare
            check_coordinates: Whether to check spatial coordinate compatibility
            check_gene_overlap: Whether to check gene/feature overlap

        Returns:
            Compatibility analysis report
        """
        try:
            # Parse file paths
            try:
                paths = json.loads(file_paths)
                if not isinstance(paths, list) or len(paths) < 2:
                    return "❌ Need at least 2 file paths for compatibility checking"
            except json.JSONDecodeError:
                return "❌ Invalid JSON format for file_paths"

            # Validate all files
            results = self.validator.validate_multiple_files(paths, ValidationLevel.DOMAIN)

            # Check if all files are valid enough for comparison
            valid_results = [r for r in results if not r.has_errors]
            if len(valid_results) < 2:
                return f"❌ Need at least 2 valid files for compatibility checking. Valid files: {len(valid_results)}/{len(results)}"

            # Perform compatibility analysis
            compatibility = f"🔗 Spatial Data Compatibility Analysis\n\n"
            compatibility += f"Analyzing {len(valid_results)} files:\n"
            for i, result in enumerate(valid_results, 1):
                compatibility += f"  {i}. {Path(result.file_path).name} ({result.format_type})\n"

            compatibility += "\n**Format Compatibility:**\n"
            formats = [r.format_type for r in valid_results]
            unique_formats = set(formats)
            if len(unique_formats) == 1:
                compatibility += f"  ✅ All files use the same format: {list(unique_formats)[0]}\n"
            else:
                compatibility += f"  ⚠️ Mixed formats detected: {', '.join(unique_formats)}\n"
                compatibility += "     Consider converting to a common format for analysis\n"

            # Check dimensions if available
            compatibility += "\n**Dimension Compatibility:**\n"
            shapes = []
            for result in valid_results:
                if 'expression_matrix_shape' in result.metadata:
                    shapes.append(result.metadata['expression_matrix_shape'])
                elif 'n_obs' in result.metadata and 'n_vars' in result.metadata:
                    shapes.append([result.metadata['n_obs'], result.metadata['n_vars']])

            if shapes:
                gene_counts = [s[1] if len(s) > 1 else None for s in shapes]
                gene_counts = [g for g in gene_counts if g is not None]

                if gene_counts:
                    if len(set(gene_counts)) == 1:
                        compatibility += f"  ✅ Consistent gene count: {gene_counts[0]} genes\n"
                    else:
                        compatibility += f"  ⚠️ Different gene counts: {gene_counts}\n"
                        compatibility += "     May need gene filtering or harmonization\n"

            # Check spatial coordinate availability
            if check_coordinates:
                compatibility += "\n**Spatial Coordinate Compatibility:**\n"
                spatial_info = []
                for result in valid_results:
                    has_spatial = (
                        'spatial_keys' in result.metadata or
                        'components' in result.metadata and
                        any(comp in ['images', 'points', 'shapes'] for comp in result.metadata['components'])
                    )
                    spatial_info.append(has_spatial)

                if all(spatial_info):
                    compatibility += "  ✅ All files contain spatial coordinate information\n"
                elif any(spatial_info):
                    compatibility += "  ⚠️ Some files missing spatial coordinates\n"
                    compatibility += "     Spatial analysis may be limited\n"
                else:
                    compatibility += "  ❌ No spatial coordinate information found\n"
                    compatibility += "     Files may not be suitable for spatial analysis\n"

            return compatibility

        except Exception as e:
            logger.error(f"Compatibility check failed: {e}")
            return f"❌ Compatibility check failed: {e}"

    def extract_bioinformatics_metadata(
        self,
        file_path: str,
        include_quality_assessment: bool = True
    ) -> str:
        """Extract comprehensive metadata from bioinformatics files.

        Args:
            file_path: Path to the bioinformatics file (Nextflow, Viash, spatial data, etc.)
            include_quality_assessment: Whether to include quality assessment and suggestions

        Returns:
            Formatted metadata extraction results
        """
        try:
            result = self.metadata_extractor.extract_metadata(file_path)

            # Format the results
            report = f"🔍 Bioinformatics Metadata Extraction\n\n"
            report += f"**File:** {Path(result.file_path).name}\n"
            report += f"**Type:** {result.file_type}\n"
            report += f"**Extraction Status:** {'✅ Success' if result.extraction_success else '❌ Failed'}\n\n"

            if result.metadata_fields:
                # Group metadata by category
                categories = {}
                for field in result.metadata_fields:
                    if field.category not in categories:
                        categories[field.category] = []
                    categories[field.category].append(field)

                report += "**Extracted Metadata:**\n"
                for category, fields in categories.items():
                    report += f"\n*{category.replace('_', ' ').title()}:*\n"
                    for field in fields:
                        importance_icon = {
                            'critical': '🔴',
                            'high': '🟠',
                            'medium': '🟡',
                            'low': '🟢'
                        }.get(field.importance, '⚪')

                        report += f"  {importance_icon} **{field.name}**: {field.value}\n"
                        if field.description:
                            report += f"     _{field.description}_\n"

            if include_quality_assessment and result.quality_metrics:
                report += f"\n**Quality Assessment:**\n"
                for metric, value in result.quality_metrics.items():
                    if isinstance(value, float):
                        report += f"  • {metric.replace('_', ' ').title()}: {value:.2f}\n"
                    else:
                        report += f"  • {metric.replace('_', ' ').title()}: {value}\n"

            if result.issues:
                report += f"\n**Issues ({len(result.issues)}):**\n"
                for issue in result.issues:
                    report += f"  ❌ {issue}\n"

            if result.suggestions:
                report += f"\n**Suggestion {result.suggestions}:**\n"
                for suggestion in result.suggestions:
                    report += f"  💡 {suggestion}\n"

            return report

        except Exception as e:
            logger.error(f"Metadata extraction failed: {e}")
            return f"❌ Metadata extraction failed: {e}"

    def analyze_workflow_configuration(
        self,
        file_path: str,
        check_dependencies: bool = True,
        validate_structure: bool = True
    ) -> str:
        """Analyze Nextflow or Viash workflow configuration files.

        Args:
            file_path: Path to the workflow configuration file
            check_dependencies: Whether to analyze dependencies and requirements
            validate_structure: Whether to validate configuration structure

        Returns:
            Formatted workflow configuration analysis
        """
        try:
            result = self.metadata_extractor.extract_metadata(file_path)

            if not result.extraction_success:
                return f"❌ Failed to analyze workflow configuration: {', '.join(result.issues)}"

            analysis = f"⚙️ Workflow Configuration Analysis\n\n"
            analysis += f"**File:** {Path(result.file_path).name}\n"
            analysis += f"**Type:** {result.file_type.replace('_', ' ').title()}\n\n"

            # Extract workflow-specific information
            workflow_fields = [f for f in result.metadata_fields if f.category in ['workflow', 'component', 'interface']]
            if workflow_fields:
                analysis += "**Workflow Structure:**\n"
                for field in workflow_fields:
                    if field.name in ['component_name', 'processes', 'process_calls']:
                        analysis += f"  • {field.name.replace('_', ' ').title()}: {field.value}\n"
                    elif field.name in ['argument_count', 'input_arguments', 'output_arguments']:
                        analysis += f"  • {field.name.replace('_', ' ').title()}: {field.value}\n"

            # Extract resource and platform information
            resource_fields = [f for f in result.metadata_fields if f.category in ['resources', 'platforms', 'containers']]
            if resource_fields:
                analysis += "\n**Resources & Platforms:**\n"
                for field in resource_fields:
                    if isinstance(field.value, list) and len(field.value) > 0:
                        analysis += f"  • {field.name.replace('_', ' ').title()}: {', '.join(map(str, field.value))}\n"
                    elif not isinstance(field.value, list):
                        analysis += f"  • {field.name.replace('_', ' ').title()}: {field.value}\n"

            # Extract dependency information
            if check_dependencies:
                dep_fields = [f for f in result.metadata_fields if f.category in ['dependencies', 'tools']]
                if dep_fields:
                    analysis += "\n**Dependencies:**\n"
                    for field in dep_fields:
                        if isinstance(field.value, list):
                            analysis += f"  • {field.name.replace('_', ' ').title()}: {', '.join(field.value)}\n"
                        else:
                            analysis += f"  • {field.name.replace('_', ' ').title()}: {field.value}\n"

            # Add validation results
            if validate_structure:
                critical_fields = result.get_critical_fields()
                if not critical_fields:
                    analysis += "\n⚠️ **Structure Warning:** No critical workflow components detected\n"

                if result.quality_metrics.get('metadata_completeness', 0) < 0.5:
                    analysis += "⚠️ **Completeness Warning:** Configuration may be incomplete\n"

            # Add suggestions
            if result.suggestions:
                analysis += f"\n**Recommendations:**\n"
                for suggestion in result.suggestions:
                    analysis += f"  💡 {suggestion}\n"

            return analysis

        except Exception as e:
            logger.error(f"Workflow configuration analysis failed: {e}")
            return f"❌ Workflow configuration analysis failed: {e}"

    def assess_data_quality(
        self,
        file_paths: str,  # JSON string of file paths
        include_spatial_validation: bool = True,
        include_metadata_analysis: bool = True
    ) -> str:
        """Assess data quality across multiple bioinformatics files.

        Args:
            file_paths: JSON string containing list of file paths to assess
            include_spatial_validation: Whether to include spatial data validation
            include_metadata_analysis: Whether to include metadata quality analysis

        Returns:
            Comprehensive data quality assessment report
        """
        try:
            # Parse file paths
            try:
                paths = json.loads(file_paths)
                if not isinstance(paths, list):
                    return "❌ file_paths must be a JSON list of file paths"
            except json.JSONDecodeError:
                return "❌ Invalid JSON format for file_paths"

            assessment = f"📊 Data Quality Assessment Report\n\n"
            assessment += f"**Files Analyzed:** {len(paths)}\n\n"

            # Perform metadata extraction for all files
            metadata_results = []
            for file_path in paths:
                result = self.metadata_extractor.extract_metadata(file_path)
                metadata_results.append(result)

            # Overall statistics
            successful_extractions = len([r for r in metadata_results if r.extraction_success])
            total_issues = sum(len(r.issues) for r in metadata_results)

            assessment += f"**Overall Statistics:**\n"
            assessment += f"  • Successful Analyses: {successful_extractions}/{len(paths)}\n"
            assessment += f"  • Total Issues Found: {total_issues}\n"

            # File type distribution
            file_types = {}
            for result in metadata_results:
                file_types[result.file_type] = file_types.get(result.file_type, 0) + 1

            if file_types:
                assessment += f"\n**File Type Distribution:**\n"
                for file_type, count in file_types.items():
                    assessment += f"  • {file_type.replace('_', ' ').title()}: {count}\n"

            # Quality metrics aggregation
            if include_metadata_analysis:
                completeness_scores = [r.quality_metrics.get('metadata_completeness', 0)
                                     for r in metadata_results if r.quality_metrics]
                if completeness_scores:
                    avg_completeness = sum(completeness_scores) / len(completeness_scores)
                    assessment += f"\n**Metadata Quality:**\n"
                    assessment += f"  • Average Completeness: {avg_completeness:.1%}\n"

                    low_quality_files = [r for r in metadata_results
                                       if r.quality_metrics.get('metadata_completeness', 0) < 0.3]
                    if low_quality_files:
                        assessment += f"  • Files Needing Attention: {len(low_quality_files)}\n"

            # Spatial data validation if requested
            if include_spatial_validation:
                spatial_files = [r for r in metadata_results if 'spatial' in r.file_type]
                if spatial_files:
                    assessment += f"\n**Spatial Data Quality:**\n"

                    # Validate spatial files
                    spatial_paths = [r.file_path for r in spatial_files]
                    validation_results = self.validator.validate_multiple_files(spatial_paths, ValidationLevel.DOMAIN)

                    valid_spatial = len([r for r in validation_results if r.is_valid])
                    assessment += f"  • Valid Spatial Files: {valid_spatial}/{len(spatial_files)}\n"

                    spatial_issues = sum(len(r.issues) for r in validation_results)
                    assessment += f"  • Spatial Validation Issues: {spatial_issues}\n"

            # Individual file quality summary
            assessment += f"\n**Individual File Quality:**\n"
            for i, result in enumerate(metadata_results, 1):
                status_icon = "✅" if result.extraction_success and len(result.issues) == 0 else "⚠️" if result.extraction_success else "❌"
                quality_score = result.quality_metrics.get('metadata_completeness', 0)

                assessment += f"{i}. {status_icon} {Path(result.file_path).name} "
                assessment += f"(Quality: {quality_score:.1%}, Issues: {len(result.issues)})\n"

            # Aggregate suggestions
            all_suggestions = []
            for result in metadata_results:
                all_suggestions.extend(result.suggestions)

            if all_suggestions:
                # Get unique suggestions
                unique_suggestions = list(set(all_suggestions))
                assessment += f"\n**Overall Recommendations:**\n"
                for suggestion in unique_suggestions[:5]:  # Top 5 suggestions
                    assessment += f"  💡 {suggestion}\n"

            return assessment

        except Exception as e:
            logger.error(f"Data quality assessment failed: {e}")
            return f"❌ Data quality assessment failed: {e}"

    def analyze_workflow_dependencies(
        self,
        file_paths: str,  # JSON string of file paths
        include_containers: bool = True,
        include_libraries: bool = True
    ) -> str:
        """Analyze dependencies across multiple workflow files.

        Args:
            file_paths: JSON string containing list of workflow file paths
            include_containers: Whether to analyze container dependencies
            include_libraries: Whether to analyze library dependencies

        Returns:
            Comprehensive workflow dependency analysis
        """
        try:
            # Parse file paths
            try:
                paths = json.loads(file_paths)
                if not isinstance(paths, list):
                    return "❌ file_paths must be a JSON list of file paths"
            except json.JSONDecodeError:
                return "❌ Invalid JSON format for file_paths"

            # Perform dependency analysis
            dependency_analysis = self.metadata_extractor.analyze_workflow_dependencies(paths)

            report = f"🔗 Workflow Dependency Analysis\n\n"
            report += f"**Analysis Summary:**\n"
            report += f"  • Files Analyzed: {dependency_analysis['total_files_analyzed']}\n"
            report += f"  • Successful Extractions: {dependency_analysis['successful_extractions']}\n"
            report += f"  • File Types: {', '.join(dependency_analysis['file_types'])}\n"
            report += f"  • Total Issues: {dependency_analysis['total_issues']}\n\n"

            if include_containers and dependency_analysis['unique_containers']:
                report += f"**Container Dependencies ({dependency_analysis['container_count']}):**\n"
                for container in dependency_analysis['unique_containers']:
                    report += f"  • {container}\n"
                report += "\n"

            if dependency_analysis['unique_processes']:
                report += f"**Workflow Processes ({dependency_analysis['process_count']}):**\n"
                for process in dependency_analysis['unique_processes']:
                    report += f"  • {process}\n"
                report += "\n"

            if include_libraries and dependency_analysis['unique_libraries']:
                report += f"**Library Dependencies ({dependency_analysis['library_count']}):**\n"
                for library in dependency_analysis['unique_libraries']:
                    report += f"  • {library}\n"
                report += "\n"

            # Add recommendations
            report += "**Dependency Recommendations:**\n"

            if dependency_analysis['container_count'] == 0:
                report += "  💡 Consider adding container specifications for reproducibility\n"

            if dependency_analysis['container_count'] > 10:
                report += "  ⚠️ Large number of containers - consider consolidation\n"

            if dependency_analysis['total_issues'] > 0:
                report += f"  ⚠️ {dependency_analysis['total_issues']} issues found - review individual files\n"

            return report

        except Exception as e:
            logger.error(f"Workflow dependency analysis failed: {e}")
            return f"❌ Workflow dependency analysis failed: {e}"

    def _format_result_summary(self, result: ValidationResult) -> str:
        """Format validation result as a summary."""
        status_icon = "✅" if result.is_valid else "❌"

        summary = f"{status_icon} Spatial Data Validation Summary\n\n"
        summary += f"**File:** {Path(result.file_path).name}\n"
        summary += f"**Format:** {result.format_type}\n"
        summary += f"**Size:** {self._format_file_size(result.file_size)}\n"
        summary += f"**Status:** {'Valid' if result.is_valid else 'Invalid'}\n"
        summary += f"**Validation Level:** {result.validation_level.value}\n\n"

        if result.issues:
            summary += f"**Issues Found:** {len(result.issues)}\n"

            # Group issues by level
            errors = [i for i in result.issues if i.level in [ValidationStatus.ERROR, ValidationStatus.CRITICAL]]
            warnings = [i for i in result.issues if i.level == ValidationStatus.WARNING]

            if errors:
                summary += f"  • Errors: {len(errors)}\n"
            if warnings:
                summary += f"  • Warnings: {len(warnings)}\n"

            summary += "\n**Key Issues:**\n"
            for issue in result.issues[:5]:  # Show first 5 issues
                level_icon = "❌" if issue.level in [ValidationStatus.ERROR, ValidationStatus.CRITICAL] else "⚠️"
                summary += f"{level_icon} {issue.message}\n"
                if issue.suggested_fix:
                    summary += f"   💡 {issue.suggested_fix}\n"

            if len(result.issues) > 5:
                summary += f"   ... and {len(result.issues) - 5} more issues\n"
        else:
            summary += "**Issues Found:** None\n"

        return summary

    def _format_result_detailed(self, result: ValidationResult) -> str:
        """Format validation result with full details."""
        status_icon = "✅" if result.is_valid else "❌"

        detailed = f"{status_icon} Detailed Spatial Data Validation Report\n\n"
        detailed += f"**File Information:**\n"
        detailed += f"  • Path: {result.file_path}\n"
        detailed += f"  • Format: {result.format_type}\n"
        detailed += f"  • Size: {self._format_file_size(result.file_size)}\n"
        detailed += f"  • Status: {'Valid' if result.is_valid else 'Invalid'}\n"
        detailed += f"  • Validation Level: {result.validation_level.value}\n\n"

        if result.metadata:
            detailed += "**Metadata:**\n"
            for key, value in result.metadata.items():
                detailed += f"  • {key}: {value}\n"
            detailed += "\n"

        if result.issues:
            detailed += f"**Issues ({len(result.issues)}):**\n"
            for i, issue in enumerate(result.issues, 1):
                level_icon = {
                    ValidationStatus.CRITICAL: "🔴",
                    ValidationStatus.ERROR: "❌",
                    ValidationStatus.WARNING: "⚠️",
                    ValidationStatus.VALID: "✅"
                }.get(issue.level, "ℹ️")

                detailed += f"{i}. {level_icon} **{issue.category.upper()}**: {issue.message}\n"
                if issue.details:
                    detailed += f"   Details: {issue.details}\n"
                if issue.suggested_fix:
                    detailed += f"   💡 Suggested Fix: {issue.suggested_fix}\n"
                detailed += "\n"
        else:
            detailed += "**Issues:** None found\n"

        return detailed

    def _format_result_json(self, result: ValidationResult) -> str:
        """Format validation result as JSON."""
        return json.dumps(self._result_to_dict(result), indent=2)

    def _result_to_dict(self, result: ValidationResult) -> Dict[str, Any]:
        """Convert ValidationResult to dictionary."""
        return {
            'file_path': result.file_path,
            'format_type': result.format_type,
            'is_valid': result.is_valid,
            'validation_level': result.validation_level.value,
            'file_size': result.file_size,
            'metadata': result.metadata,
            'issues': [
                {
                    'level': issue.level.value,
                    'category': issue.category,
                    'message': issue.message,
                    'details': issue.details,
                    'suggested_fix': issue.suggested_fix,
                    'file_path': issue.file_path,
                    'line_number': issue.line_number
                }
                for issue in result.issues
            ]
        }

    def _format_multiple_results_summary(
        self,
        results: List[ValidationResult],
        summary: Dict[str, Any]
    ) -> str:
        """Format multiple validation results as summary."""
        report = f"📊 Multiple File Validation Summary\n\n"
        report += f"**Overall Statistics:**\n"
        report += f"  • Total Files: {summary['total_files']}\n"
        report += f"  • Valid Files: {summary['valid_files']}\n"
        report += f"  • Files with Warnings: {summary['files_with_warnings']}\n"
        report += f"  • Files with Errors: {summary['files_with_errors']}\n"
        report += f"  • Success Rate: {summary['success_rate']:.1%}\n\n"

        if summary['format_distribution']:
            report += "**Format Distribution:**\n"
            for format_type, count in summary['format_distribution'].items():
                report += f"  • {format_type}: {count} files\n"
            report += "\n"

        report += "**File Results:**\n"
        for i, result in enumerate(results, 1):
            status_icon = "✅" if result.is_valid else "❌"
            issue_count = len(result.issues)
            report += f"{i}. {status_icon} {Path(result.file_path).name} ({result.format_type})"
            if issue_count > 0:
                report += f" - {issue_count} issues"
            report += "\n"

        return report

    def _format_multiple_results_detailed(
        self,
        results: List[ValidationResult],
        summary: Dict[str, Any]
    ) -> str:
        """Format multiple validation results with details."""
        report = self._format_multiple_results_summary(results, summary)

        # Add detailed issues
        report += "\n**Detailed Issues:**\n"
        for result in results:
            if result.issues:
                report += f"\n**{Path(result.file_path).name}:**\n"
                report += self._format_issues(result.issues)

        return report

    def _format_issues(self, issues: List) -> str:
        """Format a list of validation issues."""
        if not issues:
            return "  No issues found\n"

        formatted = ""
        for issue in issues:
            level_icon = {
                ValidationStatus.CRITICAL: "🔴",
                ValidationStatus.ERROR: "❌",
                ValidationStatus.WARNING: "⚠️",
                ValidationStatus.VALID: "✅"
            }.get(issue.level, "ℹ️")

            formatted += f"  {level_icon} {issue.message}\n"
            if issue.suggested_fix:
                formatted += f"     💡 {issue.suggested_fix}\n"

        return formatted

    def _format_file_size(self, size_bytes: int) -> str:
        """Format file size in human readable format."""
        if size_bytes == 0:
            return "0 B"

        for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
            if size_bytes < 1024.0:
                return f"{size_bytes:.1f} {unit}"
            size_bytes /= 1024.0

        return f"{size_bytes:.1f} PB"
