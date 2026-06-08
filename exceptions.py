"""
Pipeline Exceptions Module
Defines custom exception classes for the metabolic pathway pipeline.
"""

from typing import Any, Dict, Optional


class PipelineError(Exception):
    """Base exception class for all pipeline errors."""
    
    def __init__(
        self,
        message: str,
        stage: Optional[int] = None,
        function_name: Optional[str] = None,
        input_json: Optional[Dict[str, Any]] = None,
        original_exception: Optional[Exception] = None
    ):
        super().__init__(message)
        self.message = message
        self.stage = stage
        self.function_name = function_name
        self.input_json = input_json or {}
        self.original_exception = original_exception
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert exception to dictionary for logging."""
        return {
            "exception_type": self.__class__.__name__,
            "message": self.message,
            "stage": self.stage,
            "function_name": self.function_name,
            "input_json_summary": str(self.input_json)[:500],
            "original_exception": str(self.original_exception) if self.original_exception else None
        }


class DataIngestionError(PipelineError):
    """Raised when data ingestion fails."""
    
    def __init__(
        self,
        message: str,
        data_source: Optional[str] = None,
        **kwargs
    ):
        super().__init__(message, **kwargs)
        self.data_source = data_source


class SchemaValidationError(PipelineError):
    """Raised when JSON schema validation fails."""
    
    def __init__(
        self,
        message: str,
        schema_name: Optional[str] = None,
        validation_errors: Optional[list] = None,
        **kwargs
    ):
        super().__init__(message, **kwargs)
        self.schema_name = schema_name
        self.validation_errors = validation_errors or []
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert exception to dictionary for logging."""
        base_dict = super().to_dict()
        base_dict["schema_name"] = self.schema_name
        base_dict["validation_errors"] = [str(e) for e in self.validation_errors]
        return base_dict


class ModelInferenceError(PipelineError):
    """Raised when ML model inference fails."""
    
    def __init__(
        self,
        message: str,
        model_name: Optional[str] = None,
        fallback_used: bool = False,
        **kwargs
    ):
        super().__init__(message, **kwargs)
        self.model_name = model_name
        self.fallback_used = fallback_used


class FBAConvergenceError(PipelineError):
    """Raised when Flux Balance Analysis fails to converge."""
    
    def __init__(
        self,
        message: str,
        objective_value: Optional[float] = None,
        status: Optional[str] = None,
        **kwargs
    ):
        super().__init__(message, **kwargs)
        self.objective_value = objective_value
        self.status = status


class FermentationSimulationError(PipelineError):
    """Raised when fermentation simulation fails."""
    
    def __init__(
        self,
        message: str,
        simulation_time: Optional[float] = None,
        ode_error: Optional[str] = None,
        **kwargs
    ):
        super().__init__(message, **kwargs)
        self.simulation_time = simulation_time
        self.ode_error = ode_error


class RetrosynthesisError(PipelineError):
    """Raised when retrosynthesis pathway generation fails."""
    
    def __init__(
        self,
        message: str,
        target_molecule: Optional[str] = None,
        num_pathways_found: int = 0,
        **kwargs
    ):
        super().__init__(message, **kwargs)
        self.target_molecule = target_molecule
        self.num_pathways_found = num_pathways_found


class EnzymeSelectionError(PipelineError):
    """Raised when enzyme selection fails."""
    
    def __init__(
        self,
        message: str,
        reaction_id: Optional[str] = None,
        organism: Optional[str] = None,
        **kwargs
    ):
        super().__init__(message, **kwargs)
        self.reaction_id = reaction_id
        self.organism = organism


class CodonOptimizationError(PipelineError):
    """Raised when codon optimization fails."""
    
    def __init__(
        self,
        message: str,
        gene_name: Optional[str] = None,
        sequence_length: Optional[int] = None,
        **kwargs
    ):
        super().__init__(message, **kwargs)
        self.gene_name = gene_name
        self.sequence_length = sequence_length


class ToxicityPredictionError(PipelineError):
    """Raised when toxicity prediction fails."""
    
    def __init__(
        self,
        message: str,
        metabolite: Optional[str] = None,
        **kwargs
    ):
        super().__init__(message, **kwargs)
        self.metabolite = metabolite


class ScaleUpError(PipelineError):
    """Raised when scale-up prediction fails."""
    
    def __init__(
        self,
        message: str,
        scale_level: Optional[float] = None,
        parameter: Optional[str] = None,
        **kwargs
    ):
        super().__init__(message, **kwargs)
        self.scale_level = scale_level
        self.parameter = parameter


class DownstreamProcessingError(PipelineError):
    """Raised when downstream processing simulation fails."""
    
    def __init__(
        self,
        message: str,
        step_name: Optional[str] = None,
        recovery_percent: Optional[float] = None,
        **kwargs
    ):
        super().__init__(message, **kwargs)
        self.step_name = step_name
        self.recovery_percent = recovery_percent


class RegulatoryComplianceError(PipelineError):
    """Raised when regulatory compliance check fails."""
    
    def __init__(
        self,
        message: str,
        regulation_type: Optional[str] = None,
        severity: str = "WARNING",
        **kwargs
    ):
        super().__init__(message, **kwargs)
        self.regulation_type = regulation_type
        self.severity = severity


def create_error_context(
    exception: Exception,
    stage: int,
    function_name: str,
    input_data: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    Create a standardized error context dictionary for logging.
    
    Args:
        exception: The exception that was raised
        stage: The pipeline stage number (1-5)
        function_name: Name of the function where error occurred
        input_data: Input data that caused the error
    
    Returns:
        Dictionary with error context for logging
    """
    import traceback
    from datetime import datetime
    
    return {
        "timestamp": datetime.utcnow().isoformat(),
        "error_type": exception.__class__.__name__,
        "error_message": str(exception),
        "stage": stage,
        "function_name": function_name,
        "input_data_summary": str(input_data)[:1000] if input_data else None,
        "traceback": traceback.format_exc(),
        "pipeline_status": "FAILED"
    }
