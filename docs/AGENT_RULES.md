# OpenProblems Spatial Transcriptomics Agent Rules

## Build & Development Commands

### Viash Component Development
- **Use `viash run` for executing components**: `viash run src/methods/component_name/config.vsh.yaml -- --input_train data.h5ad --output result.h5ad`
- **Build components with Docker engine**: Always specify `--engine docker` for consistent environments
- **Test individual components**: Use `viash test src/methods/component_name/config.vsh.yaml` before integration
- **Run parallel testing**: Execute `viash ns test --parallel --engine docker` for comprehensive validation
- **Validate configurations**: Every component must have a valid `config.vsh.yaml` file
- **Use test data**: Always test with resources from `resources_test/` directory first

### Nextflow Workflow Commands
- **Run workflows locally**: Use `nextflow run workflow.nf` with proper parameters
- **Validate pipeline syntax**: Execute `nextflow config workflow.nf` to check configuration
- **Use profiles**: Specify appropriate profiles with `-profile docker,test` for development
- **Monitor execution**: Use `nextflow log` to track workflow progress and debug issues
- **Resume failed runs**: Apply `-resume` flag to continue from last successful checkpoint

### Docker Integration Commands
- **Build component images**: Use Docker engine through Viash for consistency
- **Test containerized components**: Verify all dependencies are included in containers
- **Push to registries**: Use standardized tagging conventions for component images
- **Validate environments**: Ensure Python/R environments match OpenProblems specifications

## Testing Guidelines

### Component Testing Strategy
- **Run unit tests first**: Execute `viash test` on individual components before integration
- **Test with multiple datasets**: Validate components work across different spatial datasets
- **Validate input/output formats**: Ensure h5ad files maintain proper structure and metadata
- **Test edge cases**: Include empty datasets, single-cell data, and boundary conditions
- **Verify Docker builds**: Confirm all components build successfully in containerized environments

### Integration Testing Approach
- **Test complete workflows**: Run end-to-end pipelines with realistic data sizes
- **Validate metric calculations**: Ensure accuracy metrics produce expected ranges and distributions
- **Test control methods**: Verify positive and negative controls behave as expected
- **Cross-validate results**: Compare outputs across different methods for consistency
- **Performance benchmarking**: Measure execution time and memory usage for scalability

### Quality Assurance Checklist
- **Check GitHub Actions**: Ensure all CI/CD checks pass before merging
- **Validate test coverage**: Confirm critical code paths are tested
- **Review error handling**: Test failure modes and error message clarity
- **Verify reproducibility**: Ensure identical inputs produce identical outputs
- **Test resource requirements**: Validate memory and compute constraints are met

## Code Style & Guidelines

### Viash Component Structure
- **Follow standard layout**: Organize components with `config.vsh.yaml`, `script.py/R`, and `test.py/R`
- **Use descriptive names**: Component names should clearly indicate their function and scope
- **Define clear inputs/outputs**: Specify all required and optional parameters with types
- **Include comprehensive metadata**: Add author, description, keywords, and version information
- **Implement proper logging**: Use structured logging for debugging and monitoring

### Python Code Standards
- **Follow PEP 8**: Use consistent indentation, naming, and formatting
- **Use type hints**: Annotate function parameters and return types
- **Handle AnnData objects**: Follow scanpy/squidpy conventions for spatial data manipulation
- **Implement error handling**: Use try-catch blocks with informative error messages
- **Document functions**: Include docstrings with parameter descriptions and examples

### R Code Standards
- **Use tidyverse conventions**: Apply consistent data manipulation and visualization patterns
- **Handle Seurat objects**: Follow best practices for spatial transcriptomics analysis
- **Implement proper error handling**: Use tryCatch with meaningful error messages
- **Document functions**: Include roxygen2 documentation for all functions
- **Use consistent naming**: Apply snake_case for functions and variables

### Configuration Management
- **Use YAML for configs**: Structure configuration files with clear hierarchies
- **Define resource requirements**: Specify CPU, memory, and disk requirements accurately
- **Include version constraints**: Pin software versions for reproducibility
- **Document parameters**: Provide clear descriptions and default values
- **Validate inputs**: Implement parameter validation and type checking

## Documentation Guidelines

### Component Documentation
- **Write clear descriptions**: Explain the biological/computational problem being addressed
- **Document algorithm details**: Describe the core methodology and implementation approach
- **Provide usage examples**: Include concrete examples with sample data and parameters
- **List dependencies**: Document all required software, packages, and versions
- **Include references**: Cite relevant papers and methodological sources

### Task Documentation Structure
- **Define task motivation**: Explain the biological significance and research gaps addressed
- **Describe datasets**: Detail input data types, formats, and expected characteristics
- **Outline methods**: List implemented methods with brief algorithmic descriptions
- **Specify metrics**: Define evaluation criteria and interpretation guidelines
- **Document controls**: Explain positive and negative control implementations

### Workflow Documentation
- **Create process diagrams**: Visualize workflow steps and data flow
- **Document parameters**: Explain all configurable options and their effects
- **Provide troubleshooting**: Include common issues and resolution strategies
- **List output formats**: Describe all generated files and their contents
- **Include performance notes**: Document expected runtime and resource usage

### API Documentation Standards
- **Use OpenAPI specifications**: Document REST endpoints with complete schemas
- **Provide request/response examples**: Include realistic data samples
- **Document error codes**: Explain all possible error conditions and responses
- **Include authentication**: Detail security requirements and token usage
- **Maintain versioning**: Document API changes and backwards compatibility

## Collaboration & Review Guidelines

### Pull Request Standards
- **Create focused PRs**: Address single features or bug fixes per request
- **Write descriptive titles**: Clearly summarize changes and their purpose
- **Include comprehensive descriptions**: Explain motivation, changes, and testing performed
- **Add reviewers**: Tag appropriate domain experts and maintainers
- **Respond to feedback**: Address review comments promptly and thoroughly

### Code Review Process
- **Review for correctness**: Verify algorithmic implementation and logic
- **Check for consistency**: Ensure adherence to established patterns and conventions
- **Validate testing**: Confirm adequate test coverage and quality
- **Assess documentation**: Review clarity and completeness of documentation
- **Consider performance**: Evaluate computational efficiency and scalability

### Community Engagement
- **Use GitHub discussions**: Engage in technical discussions and feature planning
- **Participate in Discord**: Join real-time conversations and collaboration
- **Follow issue templates**: Use structured formats for bug reports and feature requests
- **Share knowledge**: Contribute to documentation and community resources
- **Mentor newcomers**: Help onboard new contributors to the ecosystem

## Quality Control & Validation

### Data Quality Standards
- **Validate spatial coordinates**: Ensure x,y coordinates are properly formatted and scaled
- **Check gene expression**: Verify count matrices have appropriate ranges and distributions
- **Assess metadata completeness**: Confirm required annotations and sample information
- **Test data integrity**: Validate file formats and cross-reference identifiers
- **Monitor data provenance**: Track data sources and processing steps

### Results Validation Process
- **Cross-method comparison**: Compare results across different algorithmic approaches
- **Statistical validation**: Apply appropriate statistical tests and multiple comparison corrections
- **Biological interpretation**: Ensure results align with known biological principles
- **Reproducibility testing**: Verify consistent results across multiple runs
- **External validation**: Compare against published benchmarks and literature

### Performance Monitoring
- **Track execution metrics**: Monitor runtime, memory usage, and resource consumption
- **Assess scalability**: Test performance across different data sizes and complexities
- **Monitor quality metrics**: Track accuracy, precision, recall, and domain-specific measures
- **Evaluate user experience**: Gather feedback on usability and documentation quality
- **Continuous improvement**: Regularly review and optimize component performance
