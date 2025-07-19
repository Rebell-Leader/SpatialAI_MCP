"""Tool detection and validation for the OpenProblems MCP Server."""

import os
import subprocess
import shutil
import logging
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from pathlib import Path

from .exceptions import DependencyError
from .config import ToolConfig


logger = logging.getLogger(__name__)


@dataclass
class ToolInfo:
    """Information about a detected tool."""
    name: str
    executable: str
    version: Optional[str] = None
    available: bool = False
    error_message: Optional[str] = None


class ToolDetector:
    """Detects and validates local bioinformatics tools."""

    def __init__(self, tool_config: ToolConfig):
        self.tool_config = tool_config
        self._detected_tools: Optional[Dict[str, ToolInfo]] = None

    def detect_all_tools(self) -> Dict[str, ToolInfo]:
        """Detect all configured tools."""
        if self._detected_tools is not None:
            return self._detected_tools

        tools = {}

        # Detect each tool
        tools['nextflow'] = self._detect_nextflow()
        tools['viash'] = self._detect_viash()
        tools['docker'] = self._detect_docker()
        tools['git'] = self._detect_git()
        tools['python'] = self._detect_python()

        self._detected_tools = tools
        return tools

    def _detect_nextflow(self) -> ToolInfo:
        """Detect Nextflow installation."""
        executable = self.tool_config.nextflow_executable

        # Check if executable exists
        if not shutil.which(executable):
            return ToolInfo(
                name="nextflow",
                executable=executable,
                available=False,
                error_message=f"Nextflow executable '{executable}' not found in PATH"
            )

        # Get version
        try:
            result = subprocess.run(
                [executable, '-version'],
                capture_output=True,
                text=True,
                timeout=10
            )

            if result.returncode == 0:
                # Parse version from output
                version = self._parse_nextflow_version(result.stdout)
                return ToolInfo(
                    name="nextflow",
                    executable=executable,
                    version=version,
                    available=True
                )
            else:
                return ToolInfo(
                    name="nextflow",
                    executable=executable,
                    available=False,
                    error_message=f"Nextflow version check failed: {result.stderr}"
                )

        except subprocess.TimeoutExpired:
            return ToolInfo(
                name="nextflow",
                executable=executable,
                available=False,
                error_message="Nextflow version check timed out"
            )
        except Exception as e:
            return ToolInfo(
                name="nextflow",
                executable=executable,
                available=False,
                error_message=f"Error checking Nextflow: {e}"
            )

    def _detect_viash(self) -> ToolInfo:
        """Detect Viash installation."""
        executable = self.tool_config.viash_executable

        if not shutil.which(executable):
            return ToolInfo(
                name="viash",
                executable=executable,
                available=False,
                error_message=f"Viash executable '{executable}' not found in PATH"
            )

        try:
            result = subprocess.run(
                [executable, '--version'],
                capture_output=True,
                text=True,
=10
            )

            if result.returncode == 0:
                version = result.stdout.strip()
                return ToolInfo(
                    name="viash",
                    executable=executable,
                    version=version,
                    available=True
                )
            else:
                return ToolInfo(
                    name="viash",
                    executable=executable,
                    available=False,
                    error_message=f"Viash version check failed: {result.stderr}"
                )

        except subprocess.TimeoutExpired:
            return ToolInfo(
                name="viash",
                executable=executable,
                available=False,
                error_message="Viash version check timed out"
            )
        except Exception as e:
            return ToolInfo(
                name="viash",
                executable=executable,
                available=False,
                error_message=f"Error checking Viash: {e}"
            )

    def _detect_docker(self) -> ToolInfo:
        """Detect Docker installation."""
        executable = self.tool_config.docker_executable

        if not shutil.which(executable):
            return ToolInfo(
                name="docker",
                executable=executable,
                available=False,
                error_message=f"Docker executable '{executable}' not found in PATH"
            )

        try:
            # Check if Docker daemon is running
            result = subprocess.run(
                [executable, 'version', '--format', '{{.Server.Version}}'],
                capture_output=True,
                text=True,
                timeout=10
            )

            if result.returncode == 0:
                version = result.stdout.strip()
                return ToolInfo(
                    name="docker",
                    executable=executable,
                    version=version,
                    available=True
                )
            else:
                return ToolInfo(
                    name="docker",
                    executable=executable,
                    available=False,
                    error_message="Docker daemon not running or permission denied"
                )

        except subprocess.TimeoutExpired:
            return ToolInfo(
                name="docker",
                executable=executable,
                available=False,
                error_message="Docker version check timed out"
            )
        except Exception as e:
            return ToolInfo(
                name="docker",
                executable=executable,
                available=False,
                error_message=f"Error checking Docker: {e}"
            )

    def _detect_git(self) -> ToolInfo:
        """Detect Git installation."""
        executable = self.tool_config.git_executable

        if not shutil.which(executable):
            return ToolInfo(
                name="git",
                executable=executable,
                available=False,
                error_message=f"Git executable '{executable}' not found in PATH"
            )

        try:
            result = subprocess.run(
                [executable, '--version'],
                capture_output=True,
                text=True,
                timeout=10
            )

            if result.returncode == 0:
                # Parse version from "git version X.Y.Z"
                version_line = result.stdout.strip()
                version = version_line.split()[-1] if version_line.startswith('git version') else version_line
                return ToolInfo(
                    name="git",
                    executable=executable,
                    version=version,
                    available=True
                )
            else:
                return ToolInfo(
                    name="git",
                    executable=executable,
                    available=False,
                    error_message=f"Git version check failed: {result.stderr}"
                )

        except subprocess.TimeoutExpired:
            return ToolInfo(
                name="git",
                executable=executable,
                available=False,
                error_message="Git version check timed out"
            )
        except Exception as e:
            return ToolInfo(
                name="git",
                executable=executable,
                available=False,
                error_message=f"Error checking Git: {e}"
            )

    def _detect_python(self) -> ToolInfo:
        """Detect Python installation."""
        executable = self.tool_config.python_executable

        if not shutil.which(executable):
            return ToolInfo(
                name="python",
                executable=executable,
                available=False,
                error_message=f"Python executable '{executable}' not found in PATH"
            )

        try:
            result = subprocess.run(
                [executable, '--version'],
                capture_output=True,
                text=True,
                timeout=10
            )

            if result.returncode == 0:
                # Parse version from "Python X.Y.Z"
                version_line = result.stdout.strip()
                version = version_line.split()[-1] if version_line.startswith('Python') else version_line
                return ToolInfo(
                    name="python",
                    executable=executable,
                    version=version,
                    available=True
                )
            else:
                return ToolInfo(
                    name="python",
                    executable=executable,
                    available=False,
                    error_message=f"Python version check failed: {result.stderr}"
                )

        except subprocess.TimeoutExpired:
            return ToolInfo(
                name="python",
                executable=executable,
                available=False,
                error_message="Python version check timed out"
            )
        except Exception as e:
            return ToolInfo(
                name="python",
                executable=executable,
                available=False,
                error_message=f"Error checking Python: {e}"
            )

    def _parse_nextflow_version(self, version_output: str) -> str:
        """Parse Nextflow version from version output."""
        lines = version_output.strip().split('\n')
        for line in lines:
            if 'version' in line.lower():
                # Extract version number
                parts = line.split()
                for part in parts:
                    if part[0].isdigit():
                        return part
        return "unknown"

    def get_missing_tools(self) -> List[ToolInfo]:
        """Get list of missing or unavailable tools."""
        tools = self.detect_all_tools()
        return [tool for tool in tools.values() if not tool.available]

    def get_available_tools(self) -> List[ToolInfo]:
        """Get list of available tools."""
        tools = self.detect_all_tools()
        return [tool for tool in tools.values() if tool.available]

    def validate_required_tools(self, required_tools: List[str]) -> None:
        """Validate that required tools are available."""
        tools = self.detect_all_tools()
        missing = []

        for tool_name in required_tools:
            if tool_name not in tools or not tools[tool_name].available:
                missing.append(tool_name)

        if missing:
            raise DependencyError(
                f"Required tools not available: {', '.join(missing)}",
                context={"missing_tools": missing},
                suggested_fixes=self._get_installation_suggestions(missing)
            )

    def _get_installation_suggestions(self, missing_tools: List[str]) -> List[str]:
        """Get installation suggestions for missing tools."""
        suggestions = []

        for tool in missing_tools:
            if tool == 'nextflow':
                suggestions.append("Install Nextflow: curl -s https://get.nextflow.io | bash")
            elif tool == 'viash':
                suggestions.append("Install Viash: https://viash.io/installation/")
            elif tool == 'docker':
                suggestions.append("Install Docker: https://docs.docker.com/get-docker/")
            elif tool == 'git':
                suggestions.append("Install Git: https://git-scm.com/downloads")
            elif tool == 'python':
                suggestions.append("Install Python: https://www.python.org/downloads/")

        return suggestions

    def get_health_status(self) -> Dict[str, any]:
        """Get overall health status of tool dependencies."""
        tools = self.detect_all_tools()
        available_count = sum(1 for tool in tools.values() if tool.available)
        total_count = len(tools)

        return {
            "tools_available": available_count,
            "tools_total": total_count,
            "all_tools_available": available_count == total_count,
            "tools": {name: {
                "available": tool.available,
                "version": tool.version,
                "error": tool.error_message
            } for name, tool in tools.items()}
        }
