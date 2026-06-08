"""
Logger Setup Module
Provides logging infrastructure for the metabolic pathway pipeline.
"""

import logging
import os
import sys
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional
from dataclasses import dataclass, asdict


@dataclass
class ErrorSummary:
    """Structured error reporting dataclass."""
    timestamp: str
    stage: int
    function_name: str
    exception_type: str
    exception_message: str
    input_json_summary: Optional[str] = None
    traceback: Optional[str] = None
    recovery_action: Optional[str] = None
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return asdict(self)
    
    def to_json(self) -> str:
        """Convert to JSON string."""
        return json.dumps(self.to_dict(), indent=2)


class StageFilter(logging.Filter):
    """Logging filter that adds stage context."""
    
    def __init__(self, stage: int):
        super().__init__()
        self.stage = stage
    
    def filter(self, record: logging.LogRecord) -> bool:
        """Add stage information to log record."""
        record.stage = self.stage
        return True


class PipelineLogger:
    """
    Stage-aware logger for the metabolic pathway pipeline.
    
    Provides:
    - Dual output to console and file
    - Stage context in all log messages
    - JSON contract logging
    - Error summary generation
    """
    
    _instances: Dict[str, 'PipelineLogger'] = {}
    
    def __new__(cls, pipeline_id: str, stage: int, log_dir: str = "logs") -> 'PipelineLogger':
        """Singleton pattern per pipeline_id."""
        if pipeline_id not in cls._instances:
            instance = super().__new__(cls)
            instance._initialized = False
            cls._instances[pipeline_id] = instance
        return cls._instances[pipeline_id]
    
    def __init__(self, pipeline_id: str, stage: int, log_dir: str = "logs"):
        """Initialize the logger."""
        if self._initialized:
            return
        
        self.pipeline_id = pipeline_id
        self.stage = stage
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        
        # Create logger
        self.logger = logging.getLogger(f"pipeline_{pipeline_id}_stage{stage}")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        
        # Create formatter with stage context
        formatter = logging.Formatter(
            fmt="%(asctime)s | %(levelname)-8s | STAGE:%(stage)s | %(funcName)s | %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
        
        # Console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(formatter)
        console_handler.addFilter(StageFilter(stage))
        self.logger.addHandler(console_handler)
        
        # File handler
        timestamp = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
        log_file = self.log_dir / f"pipeline_{timestamp}_stage{stage}.log"
        file_handler = logging.FileHandler(log_file, mode='w')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        file_handler.addFilter(StageFilter(stage))
        self.logger.addHandler(file_handler)
        
        self.log_file = log_file
        self._initialized = True
        
        # Log initialization
        self.logger.info(f"=== STAGE {stage} START === pipeline_id={pipeline_id}")
    
    def debug(self, msg: str) -> None:
        """Log debug message."""
        self.logger.debug(msg)
    
    def info(self, msg: str) -> None:
        """Log info message."""
        self.logger.info(msg)
    
    def warning(self, msg: str) -> None:
        """Log warning message."""
        self.logger.warning(msg)
    
    def error(self, msg: str) -> None:
        """Log error message."""
        self.logger.error(msg)
    
    def critical(self, msg: str) -> None:
        """Log critical message."""
        self.logger.critical(msg)
    
    def log_json_contract(
        self,
        data: Dict[str, Any],
        schema_name: str,
        direction: str = "output"
    ) -> None:
        """
        Log JSON contract details.
        
        Args:
            data: The JSON data being passed
            schema_name: Name of the schema
            direction: "input" or "output"
        """
        self.info(f"JSON Contract [{direction.upper()}]: {schema_name}")
        
        # Log key fields
        if "pipeline_id" in data:
            self.debug(f"  Pipeline ID: {data['pipeline_id']}")
        if "stage_1_status" in data:
            self.debug(f"  Status: {data['stage_1_status']}")
        if "pathway_candidates" in data:
            num_pathways = len(data.get("pathway_candidates", []))
            self.debug(f"  Pathway candidates: {num_pathways}")
        if "fba_results" in data:
            fba = data.get("fba_results", {})
            self.debug(f"  FBA objective: {fba.get('objective_value', 'N/A')}")
        
        # Log full JSON at DEBUG level (truncated if too long)
        json_str = json.dumps(data, indent=2, default=str)
        if len(json_str) > 2000:
            json_str = json_str[:2000] + "\n... [truncated]"
        self.debug(f"Full JSON payload:\n{json_str}")
    
    def log_error_summary(
        self,
        exception: Exception,
        function_name: str,
        input_data: Optional[Dict[str, Any]] = None,
        recovery_action: Optional[str] = None
    ) -> ErrorSummary:
        """
        Log a structured error summary.
        
        Args:
            exception: The exception that occurred
            function_name: Name of the function where error occurred
            input_data: Input data that caused the error
            recovery_action: Action taken to recover
        
        Returns:
            ErrorSummary dataclass instance
        """
        import traceback
        
        error_summary = ErrorSummary(
            timestamp=datetime.utcnow().isoformat(),
            stage=self.stage,
            function_name=function_name,
            exception_type=exception.__class__.__name__,
            exception_message=str(exception),
            input_json_summary=json.dumps(input_data, default=str)[:500] if input_data else None,
            traceback=traceback.format_exc(),
            recovery_action=recovery_action
        )
        
        self.error(f"Exception in {function_name}: {exception.__class__.__name__}: {exception}")
        if input_data:
            self.error(f"Input JSON that caused error: {json.dumps(input_data, default=str)[:500]}")
        self.error(f"Traceback:\n{error_summary.traceback}")
        if recovery_action:
            self.info(f"Recovery action: {recovery_action}")
        
        # Save error summary to file
        error_file = self.log_dir / f"errors_stage{self.stage}_{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}.json"
        try:
            with open(error_file, 'w') as f:
                json.dump(error_summary.to_dict(), f, indent=2)
            self.debug(f"Error summary saved to: {error_file}")
        except Exception as e:
            self.error(f"Failed to save error summary: {e}")
        
        return error_summary
    
    def log_stage_start(
        self,
        organism: str,
        molecule: str,
        config: Optional[Dict[str, Any]] = None
    ) -> None:
        """
        Log stage start with configuration.
        
        Args:
            organism: Organism name
            molecule: Target molecule name
            config: Optional configuration dict
        """
        self.info(f"=== STAGE {self.stage} START ===")
        self.info(f"Organism: {organism}")
        self.info(f"Target molecule: {molecule}")
        if config:
            self.debug(f"Configuration: {json.dumps(config, indent=2)}")
    
    def log_stage_complete(
        self,
        status: str,
        duration_seconds: float,
        output_summary: Optional[Dict[str, Any]] = None
    ) -> None:
        """
        Log stage completion.
        
        Args:
            status: PASS, FAIL, or WARN
            duration_seconds: Time taken for the stage
            output_summary: Optional summary of outputs
        """
        self.info(f"Duration: {duration_seconds:.2f}s")
        if output_summary:
            self.debug(f"Output summary: {json.dumps(output_summary, indent=2)}")
        self.info(f"=== STAGE {self.stage} COMPLETE === status={status}")
    
    def save_stage_summary(self, summary_data: Dict[str, Any]) -> Path:
        """
        Save stage summary to JSON file.
        
        Args:
            summary_data: Dictionary containing stage summary
        
        Returns:
            Path to the saved summary file
        """
        summary_file = self.log_dir / f"stage_{self.stage}_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary_data, f, indent=2, default=str)
        self.info(f"Stage summary saved to: {summary_file}")
        return summary_file
    
    def get_logger(self) -> logging.Logger:
        """Return the underlying logger instance."""
        return self.logger


def setup_logger(
    pipeline_id: str,
    stage: int,
    log_dir: str = "logs",
    log_level: str = "DEBUG"
) -> PipelineLogger:
    """
    Set up and return a PipelineLogger instance.
    
    Args:
        pipeline_id: Unique identifier for the pipeline run
        stage: Pipeline stage number (1-5)
        log_dir: Directory for log files
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    
    Returns:
        Configured PipelineLogger instance
    """
    logger = PipelineLogger(pipeline_id, stage, log_dir)
    logger.get_logger().setLevel(getattr(logging, log_level.upper(), logging.DEBUG))
    return logger


def get_existing_logger(pipeline_id: str) -> Optional[PipelineLogger]:
    """
    Get an existing logger for a pipeline.
    
    Args:
        pipeline_id: Unique identifier for the pipeline run
    
    Returns:
        Existing PipelineLogger or None if not found
    """
    return PipelineLogger._instances.get(pipeline_id)


if __name__ == "__main__":
    # Test logger functionality
    import uuid
    from exceptions import PipelineError
    
    print("Testing PipelineLogger...")
    
    pipeline_id = str(uuid.uuid4())
    logger = setup_logger(pipeline_id, stage=1)
    
    # Test basic logging
    logger.log_stage_start(
        organism="Escherichia coli",
        molecule="lycopene",
        config={"test": "config"}
    )
    
    logger.debug("This is a DEBUG message")
    logger.info("This is an INFO message")
    logger.warning("This is a WARNING message")
    
    # Test error logging
    try:
        raise PipelineError(
            message="Test error",
            stage=1,
            function_name="test_function",
            input_json={"test": "data"}
        )
    except PipelineError as e:
        error_summary = logger.log_error_summary(
            exception=e,
            function_name="test_function",
            input_data={"test": "data"},
            recovery_action="Using fallback method"
        )
        print(f"\nError summary created: {error_summary.to_dict()}")
    
    # Test JSON contract logging
    test_json = {
        "pipeline_id": pipeline_id,
        "stage_1_status": "PASS",
        "organism": {"name": "E. coli"},
        "target_molecule": {"name": "lycopene"}
    }
    logger.log_json_contract(test_json, "stage_1_output", "output")
    
    # Test stage completion
    logger.log_stage_complete(
        status="PASS",
        duration_seconds=1.234,
        output_summary={"result": "success"}
    )
    
    # Save stage summary
    summary = {
        "pipeline_id": pipeline_id,
        "stage": 1,
        "status": "PASS",
        "duration_seconds": 1.234,
        "tests_passed": True
    }
    summary_file = logger.save_stage_summary(summary)
    
    print(f"\n✓ Logger tests complete!")
    print(f"Log file: {logger.log_file}")
    print(f"Summary file: {summary_file}")
