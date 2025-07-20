"""Custom exceptions for the OpenProblems MCP Server."""

from typing import Dict, Any, List, Optional
from enum import Enum


class ErrorType(Enum):
    """Classification of error types for structured error handling."""
    VALIDATION_ERROR = "validation_error"
    EXECUTION_ERROR = "execution_error"
    RESOURCE_ERROR = "resource_error"
    PERMISSION_ERROR = "permission_error"
    TIMEOUT_ERROR = "timeout_error"
    DEPENDENCY_ERROR = "dependency_error"
    CONFIGURATION_ERROR = "configuration_error"


class MCPServerError(Exception):
    """Base exception for all MCP server errors."""

    def __init__(
        self,
        message: str,
        error_type: ErrorType = ErrorType.EXECUTION_ERROR,
        error_code: Optional[str] = None,
        context: Optional[Dict[str, Any]] = None,
        suggested_fixes: Optional[List[str]] = None,
        retry_possible: bool = False
    ):
        super().__init__(message)
        self.message = message
        self.error_type = error_type
        self.error_code = error_code or error_type.value
        self.context = context or {}
        self.suggested_fixes = suggested_fixes or []
        self.retry_possible = retry_possible

    def to_dict(self) -> Dict[str, Any]:
        """Convert exception to dictionary for structured error responses."""
        return {
            "error_type": self.error_type.value,
            "error_code": self.error_code,
            "message": self.message,
            "context": self.context,
            "suggested_fixes": self.suggested_fixes,
            "retry_possible": self.retry_possible
        }


class ValidationError(MCPServerError):
    """Raised when input validation fails."""

    def __init__(self, message: str, **kwargs):
        super().__init__(
            message,
            error_type=ErrorType.VALIDATION_ERROR,
            **kwargs
        )


class ExecutionError(MCPServerError):
    """Raised when tool execution fails."""

    def __init__(self, message: str, **kwargs):
        super().__init__(
            message,
            error_type=ErrorType.EXECUTION_ERROR,
            **kwargs
        )


class ResourceError(MCPServerError):
    """Raised when resource constraints are exceeded."""

    def __init__(self, message: str, **kwargs):
        super().__init__(
            message,
            error_type=ErrorType.RESOURCE_ERROR,
            **kwargs
        )


class PermissionError(MCPServerError):
    """Raised when file system permission errors occur."""

    def __init__(self, message: str, **kwargs):
        super().__init__(
            message,
            error_type=ErrorType.PERMISSION_ERROR,
            **kwargs
        )


class TimeoutError(MCPServerError):
    """Raised when operations exceed timeout limits."""

    def __init__(self, message: str, **kwargs):
        super().__init__(
            message,
            error_type=ErrorType.TIMEOUT_ERROR,
            retry_possible=True,
            **kwargs
        )


class DependencyError(MCPServerError):
    """Raised when external dependencies are missing or invalid."""

    def __init__(self, message: str, **kwargs):
        super().__init__(
            message,
            error_type=ErrorType.DEPENDENCY_ERROR,
            **kwargs
        )


class ConfigurationError(MCPServerError):
    """Raised when configuration is invalid or missing."""

    def __init__(self, message: str, **kwargs):
        super().__init__(
            message,
            error_type=ErrorType.CONFIGURATION_ERROR,
            **kwargs
        )
