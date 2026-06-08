#!/usr/bin/env python3
"""
Master Pipeline Runner for Industrial Metabolic Pathway Design.
Orchestrates all 5 stages of the DBTL pipeline.

Usage:
    python main_pipeline_runner.py --organism ecoli --molecule lycopene --cycles 3
    python main_pipeline_runner.py --organism scerevisiae --molecule vanillin --cycles 5
"""

import argparse
import json
import logging
import sys
import traceback
import uuid
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional

# Import pipeline modules
from pipeline_config import PipelineConfig, OrganismType, MoleculeType
from logger_setup import setup_logger
import logging

def get_logger(name: str):
    """Get a standard Python logger for the given module name."""
    return logging.getLogger(name)
from exceptions import PipelineError, DataIngestionError, SchemaValidationError
from data_layer import run_stage_1
from schema_validator import validate_stage_output

# Stage 2 imports (may be in stage_2_pathway_ai.py)
try:
    from stage_2_pathway_ai import run_stage_2
    STAGE_2_AVAILABLE = True
except ImportError:
    STAGE_2_AVAILABLE = False
    print("Warning: Stage 2 module not found. Will skip pathway prediction.")

# Stage 3 imports
try:
    from flux_analysis_orchestrator import run_stage_3
    STAGE_3_AVAILABLE = True
except ImportError:
    STAGE_3_AVAILABLE = False
    print("Warning: Stage 3 module not found. Will skip flux analysis.")

# Stage 4 imports
try:
    from dbtl_orchestrator import run_stage_4
    STAGE_4_AVAILABLE = True
except ImportError:
    STAGE_4_AVAILABLE = False
    print("Warning: Stage 4 module not found. Will skip DBTL simulation.")

# Stage 5 imports
try:
    from scaleup_engine import ScaleUpEngine
    from downstream_processor import DownstreamProcessor
    from regulatory_checker import RegulatoryChecker
    from pipeline_report_generator import run_stage_5_reporting
    STAGE_5_AVAILABLE = True
except ImportError:
    STAGE_5_AVAILABLE = False
    print("Warning: Stage 5 modules not found. Will skip scale-up/downstream/regulatory.")


class PipelineRunner:
    """
    Master orchestrator for the complete metabolic pathway design pipeline.
    Manages JSON handoffs between stages, error handling, and reporting.
    """
    
    def __init__(self, config: Optional[PipelineConfig] = None):
        self.config = config or PipelineConfig()
        self.logger = get_logger(__name__)
        
        # Tracking
        self.errors: List[str] = []
        self.warnings: List[str] = []
        self.stage_outputs: Dict[str, Any] = {}
        
        # Timing
        self.start_time: Optional[datetime] = None
        self.stage_times: Dict[str, float] = {}

    def run_full_pipeline(self, 
                         organism: str, 
                         target_molecule: str,
                         dbtl_cycles: int = 3,
                         output_dir: str = "./pipeline_output") -> Dict[str, Any]:
        """
        Execute the complete 5-stage pipeline.
        
        Args:
            organism: Organism name (e.g., 'ecoli', 'scerevisiae')
            target_molecule: Target product (e.g., 'lycopene', 'vanillin')
            dbtl_cycles: Number of DBTL cycles to simulate
            output_dir: Directory for output files
            
        Returns:
            Final pipeline report dictionary
        """
        self.start_time = datetime.now()
        self.logger.info("=" * 80)
        self.logger.info("STARTING FULL METABOLIC PATHWAY PIPELINE")
        self.logger.info(f"Organism: {organism}")
        self.logger.info(f"Target Molecule: {target_molecule}")
        self.logger.info(f"DBTL Cycles: {dbtl_cycles}")
        self.logger.info("=" * 80)
        
        try:
            # ========== STAGE 1: Data Layer ==========
            self.logger.info("\n" + "=" * 60)
            self.logger.info("EXECUTING STAGE 1: Core Infrastructure + Data Layer")
            self.logger.info("=" * 60)
            
            stage_1_start = datetime.now()
            try:
                stage_1_json = run_stage_1(organism, target_molecule)
                self.stage_outputs['stage_1'] = stage_1_json
                self.stage_times['stage_1'] = (datetime.now() - stage_1_start).total_seconds()
                self.logger.info(f"Stage 1 COMPLETE in {self.stage_times['stage_1']:.2f}s")
                
                # Validate output
                if not validate_stage_output(stage_1_json, 'stage_1_output'):
                    self.warnings.append("Stage 1 output schema validation failed")
                    
            except Exception as e:
                error_msg = f"Stage 1 failed: {str(e)}"
                self.logger.error(error_msg, exc_info=True)
                self.errors.append(error_msg)
                if self.config.stop_on_error:
                    raise
                # Continue with limited functionality
            
            # ========== STAGE 2: Pathway Prediction ==========
            if STAGE_2_AVAILABLE and 'stage_1' in self.stage_outputs:
                self.logger.info("\n" + "=" * 60)
                self.logger.info("EXECUTING STAGE 2: Pathway Prediction + AI Engine")
                self.logger.info("=" * 60)
                
                stage_2_start = datetime.now()
                try:
                    stage_2_json = run_stage_2(self.stage_outputs['stage_1'])
                    self.stage_outputs['stage_2'] = stage_2_json
                    self.stage_times['stage_2'] = (datetime.now() - stage_2_start).total_seconds()
                    self.logger.info(f"Stage 2 COMPLETE in {self.stage_times['stage_2']:.2f}s")
                    
                    # Check for warnings
                    if stage_2_json.get('stage_2_status') == 'WARN':
                        self.warnings.append("Stage 2 completed with warnings - limited pathway options")
                        
                except Exception as e:
                    error_msg = f"Stage 2 failed: {str(e)}"
                    self.logger.error(error_msg, exc_info=True)
                    self.errors.append(error_msg)
                    self.warnings.append("Pathway prediction failed - using default pathways")
                    # Create minimal stage 2 output for continuity
                    self.stage_outputs['stage_2'] = self._create_minimal_stage_2_output()
            else:
                self.warnings.append("Stage 2 not available - skipping pathway prediction")
                self.stage_outputs['stage_2'] = self._create_minimal_stage_2_output()
            
            # ========== STAGE 3: Flux Analysis ==========
            if STAGE_3_AVAILABLE and 'stage_2' in self.stage_outputs:
                self.logger.info("\n" + "=" * 60)
                self.logger.info("EXECUTING STAGE 3: Flux Analysis + Strain Optimization")
                self.logger.info("=" * 60)
                
                stage_3_start = datetime.now()
                try:
                    stage_3_json = run_stage_3(self.stage_outputs['stage_2'])
                    self.stage_outputs['stage_3'] = stage_3_json
                    self.stage_times['stage_3'] = (datetime.now() - stage_3_start).total_seconds()
                    self.logger.info(f"Stage 3 COMPLETE in {self.stage_times['stage_3']:.2f}s")
                    
                except Exception as e:
                    error_msg = f"Stage 3 failed: {str(e)}"
                    self.logger.error(error_msg, exc_info=True)
                    self.errors.append(error_msg)
                    self.warnings.append("Flux analysis failed - using estimated values")
                    self.stage_outputs['stage_3'] = self._create_minimal_stage_3_output()
            else:
                self.warnings.append("Stage 3 not available - skipping flux analysis")
                self.stage_outputs['stage_3'] = self._create_minimal_stage_3_output()
            
            # ========== STAGE 4: DBTL + Fermentation ==========
            if STAGE_4_AVAILABLE and 'stage_3' in self.stage_outputs:
                self.logger.info("\n" + "=" * 60)
                self.logger.info("EXECUTING STAGE 4: DBTL Loop + Fermentation Simulation")
                self.logger.info("=" * 60)
                
                stage_4_start = datetime.now()
                try:
                    stage_4_json = run_stage_4(self.stage_outputs['stage_3'], cycles=dbtl_cycles)
                    self.stage_outputs['stage_4'] = stage_4_json
                    self.stage_times['stage_4'] = (datetime.now() - stage_4_start).total_seconds()
                    self.logger.info(f"Stage 4 COMPLETE in {self.stage_times['stage_4']:.2f}s")
                    
                except Exception as e:
                    error_msg = f"Stage 4 failed: {str(e)}"
                    self.logger.error(error_msg, exc_info=True)
                    self.errors.append(error_msg)
                    self.warnings.append("DBTL simulation failed - using estimated fermentation data")
                    self.stage_outputs['stage_4'] = self._create_minimal_stage_4_output()
            else:
                self.warnings.append("Stage 4 not available - skipping DBTL simulation")
                self.stage_outputs['stage_4'] = self._create_minimal_stage_4_output()
            
            # ========== STAGE 5: Scale-up + Downstream + Regulatory ==========
            if STAGE_5_AVAILABLE and 'stage_4' in self.stage_outputs:
                self.logger.info("\n" + "=" * 60)
                self.logger.info("EXECUTING STAGE 5: Scale-up + Downstream + Regulatory")
                self.logger.info("=" * 60)
                
                stage_5_start = datetime.now()
                try:
                    # Run Scale-up
                    scaleup_engine = ScaleUpEngine(self.config)
                    stage_5_scaleup = scaleup_engine.run_stage_5_scaleup(self.stage_outputs['stage_4'])
                    
                    # Run Downstream
                    downstream_proc = DownstreamProcessor(self.config)
                    stage_5_downstream = downstream_proc.run_stage_5_downstream(self.stage_outputs['stage_4'])
                    
                    # Run Regulatory
                    reg_checker = RegulatoryChecker(self.config)
                    stage_5_regulatory = reg_checker.run_stage_5_regulatory(self.stage_outputs['stage_4'])
                    
                    self.stage_outputs['stage_5'] = {
                        'scaleup': stage_5_scaleup,
                        'downstream': stage_5_downstream,
                        'regulatory': stage_5_regulatory
                    }
                    self.stage_times['stage_5'] = (datetime.now() - stage_5_start).total_seconds()
                    self.logger.info(f"Stage 5 COMPLETE in {self.stage_times['stage_5']:.2f}s")
                    
                except Exception as e:
                    error_msg = f"Stage 5 failed: {str(e)}"
                    self.logger.error(error_msg, exc_info=True)
                    self.errors.append(error_msg)
                    self.warnings.append("Stage 5 partial failure - some assessments incomplete")
            else:
                self.warnings.append("Stage 5 not available - skipping scale-up/downstream/regulatory")
            
            # ========== GENERATE FINAL REPORT ==========
            self.logger.info("\n" + "=" * 60)
            self.logger.info("GENERATING FINAL PIPELINE REPORT")
            self.logger.info("=" * 60)
            
            total_time = (datetime.now() - self.start_time).total_seconds()
            
            if STAGE_5_AVAILABLE:
                report_files = run_stage_5_reporting(
                    stage_1_json=self.stage_outputs.get('stage_1', {}),
                    stage_2_json=self.stage_outputs.get('stage_2'),
                    stage_3_json=self.stage_outputs.get('stage_3'),
                    stage_4_json=self.stage_outputs.get('stage_4'),
                    stage_5_scaleup=self.stage_outputs.get('stage_5', {}).get('scaleup') if 'stage_5' in self.stage_outputs else None,
                    stage_5_downstream=self.stage_outputs.get('stage_5', {}).get('downstream') if 'stage_5' in self.stage_outputs else None,
                    stage_5_regulatory=self.stage_outputs.get('stage_5', {}).get('regulatory') if 'stage_5' in self.stage_outputs else None,
                    errors=self.errors,
                    warnings=self.warnings,
                    output_dir=output_dir
                )
                
                self.logger.info(f"Reports saved to: {report_files}")
            else:
                # Save raw JSON outputs
                output_path = Path(output_dir)
                output_path.mkdir(parents=True, exist_ok=True)
                
                with open(output_path / "pipeline_outputs.json", 'w') as f:
                    json.dump(self.stage_outputs, f, indent=2, default=str)
                self.logger.info(f"Raw outputs saved to: {output_path / 'pipeline_outputs.json'}")
            
            # Final summary
            self.logger.info("\n" + "=" * 80)
            self.logger.info("PIPELINE EXECUTION COMPLETE")
            self.logger.info(f"Total Duration: {total_time:.2f}s")
            self.logger.info(f"Errors: {len(self.errors)}")
            self.logger.info(f"Warnings: {len(self.warnings)}")
            self.logger.info("=" * 80)
            
            return {
                'status': 'SUCCESS' if not self.errors else 'PARTIAL',
                'pipeline_id': self.stage_outputs.get('stage_1', {}).get('pipeline_id', 'unknown'),
                'total_time_seconds': total_time,
                'stage_times': self.stage_times,
                'errors': self.errors,
                'warnings': self.warnings,
                'output_directory': output_dir
            }
            
        except Exception as e:
            self.logger.critical(f"Pipeline execution failed catastrophically: {e}", exc_info=True)
            self.errors.append(f"Catastrophic failure: {str(e)}")
            return {
                'status': 'FAILED',
                'errors': self.errors,
                'traceback': traceback.format_exc()
            }

    def _create_minimal_stage_2_output(self) -> Dict[str, Any]:
        """Create minimal Stage 2 output for pipeline continuity."""
        return {
            'pipeline_id': self.stage_outputs.get('stage_1', {}).get('pipeline_id', 'unknown'),
            'stage_1_output': self.stage_outputs.get('stage_1', {}),
            'pathway_candidates': [],
            'gene_modifications': {'knockouts': [], 'overexpressions': [], 'heterologous_insertions': []},
            'stage_2_status': 'FAIL'
        }

    def _create_minimal_stage_3_output(self) -> Dict[str, Any]:
        """Create minimal Stage 3 output for pipeline continuity."""
        return {
            'pipeline_id': self.stage_outputs.get('stage_1', {}).get('pipeline_id', 'unknown'),
            'stage_2_output': self.stage_outputs.get('stage_2', {}),
            'fba_results': {},
            'strain_design': {'final_knockouts': [], 'final_overexpressions': []},
            'stage_3_status': 'FAIL'
        }

    def _create_minimal_stage_4_output(self) -> Dict[str, Any]:
        """Create minimal Stage 4 output for pipeline continuity."""
        return {
            'pipeline_id': self.stage_outputs.get('stage_1', {}).get('pipeline_id', 'unknown'),
            'stage_3_output': self.stage_outputs.get('stage_3', {}),
            'fermentation_simulation': {'final_titer_g_per_l': 0.0},
            'stage_4_status': 'FAIL'
        }


def main():
    """Main entry point with CLI argument parsing."""
    parser = argparse.ArgumentParser(
        description="Metabolic Pathway Design Pipeline for Industrial Microorganisms",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main_pipeline_runner.py --organism ecoli --molecule lycopene
  python main_pipeline_runner.py --organism scerevisiae --molecule vanillin --cycles 5
  python main_pipeline_runner.py --organism cglutamicum --molecule l_lysine --output-dir ./results
        """
    )
    
    parser.add_argument(
        '--organism', '-o',
        type=str,
        required=True,
        choices=['ecoli', 'scerevisiae', 'bsubtilis', 'cglutamicum', 'pputida'],
        help='Target organism for production'
    )
    
    parser.add_argument(
        '--molecule', '-m',
        type=str,
        required=True,
        choices=['lycopene', 'vanillin', 'l_lysine', 'l_glutamate', 'riboflavin', 'pha', 'artemisinic_acid'],
        help='Target molecule to produce'
    )
    
    parser.add_argument(
        '--cycles', '-c',
        type=int,
        default=3,
        help='Number of DBTL cycles to simulate (default: 3)'
    )
    
    parser.add_argument(
        '--output-dir',
        type=str,
        default='./pipeline_output',
        help='Directory for output files (default: ./pipeline_output)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging (DEBUG level)'
    )
    
    parser.add_argument(
        '--stop-on-error',
        action='store_true',
        help='Stop pipeline execution on first error'
    )
    
    args = parser.parse_args()
    
    # Generate pipeline ID
    pipeline_id = str(uuid.uuid4())
    
    # Setup logging
    log_level = 'DEBUG' if args.verbose else 'INFO'
    logger = setup_logger(pipeline_id=pipeline_id, stage=1, log_dir='logs', log_level=log_level)
    
    # Create config
    config = PipelineConfig(stop_on_error=args.stop_on_error)
    
    # Run pipeline
    runner = PipelineRunner(config)
    result = runner.run_full_pipeline(
        organism=args.organism,
        target_molecule=args.molecule,
        dbtl_cycles=args.cycles,
        output_dir=args.output_dir
    )
    
    # Print summary
    print("\n" + "=" * 80)
    print("PIPELINE EXECUTION SUMMARY")
    print("=" * 80)
    print(f"Status: {result['status']}")
    print(f"Pipeline ID: {result.get('pipeline_id', 'N/A')}")
    print(f"Total Time: {result.get('total_time_seconds', 0):.2f}s")
    print(f"Errors: {len(result.get('errors', []))}")
    print(f"Warnings: {len(result.get('warnings', []))}")
    print(f"Output Directory: {result.get('output_directory')}")
    print("=" * 80)
    
    # Exit code
    sys.exit(0 if result['status'] == 'SUCCESS' else 1)


if __name__ == "__main__":
    main()
