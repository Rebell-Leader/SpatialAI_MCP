#!/usr/bin/env python3
"""Setup script for OpenProblems Spatial Transcriptomics Local MCP Server."""

from setuptools import setup, find_packages
import os

# Read the README file
def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return "OpenProblems Spatial Transcriptomics Local MCP Server"

# Read version from package
def read_version():
    version_file = os.path.join('src', 'openproblems_mcp', '__init__.py')
    if os.path.exists(version_file):
        with open(version_file, 'r') as f:
            for line in f:
                if line.startswith('__version__'):
                    return line.split('=')[1].strip().strip('"').strip("'")
    return "0.1.0"

setup(
    name="openproblems-spatial-mcp",
    version=read_version(),
    description="Local Model Context Protocol server for OpenProblems spatial transcriptomics workflows",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    author="OpenProblems MCP Contributors",
    author_email="info@openproblems.bio",
    url="https://github.com/openproblems-bio/SpatialAI_MCP",
    project_urls={
        "Bug Reports": "https://github.com/openproblems-bio/SpatialAI_MCP/issues",
        "Source": "https://github.com/openproblems-bio/SpatialAI_MCP",
        "Documentation": "https://github.com/openproblems-bio/SpatialAI_MCP/docs",
    },

    # Package configuration
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.8",

    # Dependencies
    install_requires=[
        "mcp>=1.9.2",
        "pyyaml>=6.0",
        "python-dotenv>=1.0.0",
        "psutil>=5.9.0",
        "click>=8.1.0",
        "rich>=13.0.0",
        "aiofiles>=23.0.0",
    ],

    # Optional dependencies
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-asyncio>=0.21.0",
            "black>=23.0.0",
            "flake8>=6.0.0",
            "mypy>=1.0.0",
            "coverage>=7.0.0",
        ],
        "docs": [
            "mkdocs>=1.4.0",
            "mkdocs-material>=9.0.0",
            "mkdocs-mermaid2-plugin>=0.6.0",
        ],
    },

    # Entry points
    entry_points={
        "console_scripts": [
            "openproblems-mcp-server=openproblems_mcp.main:main",
            "openproblems-mcp=openproblems_mcp.cli:main",
        ],
    },

    # Package data
    include_package_data=True,
    package_data={
        "openproblems_mcp": [
            "config/*.yaml",
            "templates/*.yaml",
            "templates/*.py",
            "templates/*.R",
        ],
    },

    # Classifiers
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Operating System :: OS Independent",
    ],

    # Keywords
    keywords="mcp model-context-protocol spatial-transcriptomics bioinformatics nextflow viash openproblems",

    # License
    license="MIT",

    # Zip safe
    zip_safe=False,
)
