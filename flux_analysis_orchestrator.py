"""
Flux Analysis Orchestrator - Main entry point for Stage 3.

This module orchestrates FBA, strain optimization, and toxicity prediction
to produce the complete Stage 3 output JSON.
"""

import logging
import json
from datetime import datetime
from typing import Dict, List, Any, Optional
import uuid

from logger_setup import get_existing_logger, setup_logger
from schema_validator import validate_stage_output
from exceptions import PipelineError, SchemaValidationError
from pipeline_config import PipelineConfig
from fba_engine import FBAEngine, FBAModel
from strain_optimizer import StrainOptimizer, StrainDesign
from toxicity_predictor import ToxicityPredictor

# Global logger placeholder - will be initialized in run_stage_3
logger = None


def run_stage_3(stage_2_json: Dict[str, Any], config: Optional[PipelineConfig] = None) -> Dict[str, Any]:
    """
    Run Stage 3: Flux Analysis + Strain Optimization.
    
    This function:
    1. Validates Stage 2 input JSON against schema
    2. Extracts pathway and organism information
    3. Runs FBA to calculate fluxes
    4. Optimizes strain design using NSGA-III
    5. Assesses toxicity of intermediates
    6. Returns validated Stage 3 output JSON
    
    Args:
        stage_2_json: Validated Stage 2 output dictionary
        config: Optional pipeline configuration
        
    Returns:
        Stage 3 output dictionary matching STAGE_3_OUTPUT_SCHEMA
    """
    # Initialize
    config = config or PipelineConfig()
    start_time = datetime.now()
    pipeline_id = stage_2_json.get("pipeline_id", str(uuid.uuid4()))
    
    # Initialize logger for this pipeline run
    global logger
    logger = get_existing_logger(pipeline_id)
    if not logger:
        logger = setup_logger(pipeline_id, stage=3, log_dir="logs")
    
    logger.info(f"=== STAGE 3 START === pipeline_id={pipeline_id}")
    logger.debug(f"Input JSON received: {json.dumps(stage_2_json, default=str)[:500]}...")
    
    try:
        # Validate Stage 2 input
        logger.info("Validating Stage 2 input schema")
        validation_result = validate_stage_output(stage_2_json, "stage_2_output")
        # validation_result is a tuple: (is_valid: bool, errors: list)
        is_valid, errors = validation_result
        if not is_valid:
            logger.warning(f"Schema validation warnings: {errors}")
        
        # Extract key information from Stage 2
        stage_1_output = stage_2_json.get("stage_1_output", {})
        organism_info = stage_1_output.get("organism", {})
        target_info = stage_1_output.get("target_molecule", {})
        pathway_candidates = stage_2_json.get("pathway_candidates", [])
        gene_modifications = stage_2_json.get("gene_modifications", {})
        
        organism_name = organism_info.get("name", "ecoli")
        target_molecule = target_info.get("name", "unknown")
        
        logger.info(f"Running Stage 3 for {organism_name} producing {target_molecule}")
        
        if not pathway_candidates:
            logger.error("No pathway candidates found in Stage 2 output")
            return _create_error_stage3_output(pipeline_id, "No pathway candidates from Stage 2")
        
        # Select best pathway (rank 1)
        best_pathway = sorted(pathway_candidates, key=lambda p: p.get("rank", 999))[0]
        pathway_steps = best_pathway.get("steps", [])
        pathway_id = best_pathway.get("pathway_id", "PATHWAY_001")
        
        logger.info(f"Using pathway: {pathway_id} with {len(pathway_steps)} steps")
        
        # Initialize engines
        fba_engine = FBAEngine(config)
        strain_optimizer = StrainOptimizer(config)
        toxicity_predictor = ToxicityPredictor(config)
        
        # Step 1: Build FBA model and run flux analysis
        logger.info("Building stoichiometric matrix and running FBA")
        try:
            fba_model = fba_engine.build_stoichiometric_matrix(pathway_steps, organism_name)
            
            # Apply organism-specific constraints
            fba_model = fba_engine.organism_specific_constraints(fba_model, organism_name)
            
            # Run FBA
            product_reaction = pathway_steps[-1].get("reaction_id") if pathway_steps else None
            fba_result = fba_engine.run_fba(fba_model, product_reaction_id=product_reaction)
            
            # Run pFBA for refined fluxes
            pfba_result = fba_engine.run_pfba(fba_model)
            
            # Run FVA on key reactions
            reaction_ids = [step.get("reaction_id") for step in pathway_steps[:5]]
            fva_results = fba_engine.run_fva(fba_model, reaction_ids)
            
            logger.info(f"FBA complete: growth={fba_result.growth_rate_per_hour:.4f}/h, "
                       f"product_flux={fba_result.product_flux_mmol_per_gdw_per_hour:.4f}")
            
            fba_success = True
            fba_results_dict = fba_result.to_dict()
            
        except Exception as e:
            logger.error(f"FBA failed: {type(e).__name__}: {e}")
            logger.error("Attempting fallback: using estimated flux values")
            fba_success = False
            fba_results_dict = _create_fallback_fba_results(organism_name, target_molecule)
        
        # Step 2: Strain optimization
        logger.info("Running strain optimization")
        try:
            # Run OptKnock simulation for knockout targets
            optknock_knockouts = strain_optimizer.run_optknock_simulated(
                pathway_steps, organism_name, target_molecule, max_knockouts=5
            )
            
            # Run NSGA-III multi-objective optimization
            pareto_designs = strain_optimizer.run_nsga3_optimization(
                pathway_steps, organism_name, target_molecule,
                population_size=30, generations=20
            )
            
            # Rank designs and select best
            ranked_designs = strain_optimizer.rank_strain_designs(pareto_designs)
            best_design = ranked_designs[0] if ranked_designs else None
            
            if best_design:
                # Combine OptKnock and NSGA-III results
                final_knockouts = list(set(optknock_knockouts + best_design.knockouts))[:6]
                final_overexpressions = best_design.overexpressions[:5]
                
                strain_design_dict = {
                    "algorithm_used": "NSGA-III + OptKnock",
                    "final_knockouts": final_knockouts,
                    "final_overexpressions": final_overexpressions,
                    "predicted_titer_g_per_l": best_design.predicted_titer_g_per_l,
                    "predicted_productivity_g_per_l_per_h": best_design.predicted_productivity_g_per_l_per_h,
                    "metabolic_burden_score": best_design.metabolic_burden_score,
                    "genetic_stability_score": best_design.genetic_stability_score
                }
                
                logger.info(f"Best strain design: titer={best_design.predicted_titer_g_per_l:.2f} g/L")
            else:
                strain_design_dict = _create_fallback_strain_design(organism_name, target_molecule)
                
        except Exception as e:
            logger.error(f"Strain optimization failed: {type(e).__name__}: {e}")
            logger.error("Using fallback strain design")
            strain_design_dict = _create_fallback_strain_design(organism_name, target_molecule)
        
        # Step 3: Toxicity assessment
        logger.info("Assessing pathway toxicity")
        try:
            toxicity_result = toxicity_predictor.assess_pathway_toxicity(pathway_steps, organism_name)
            
            toxicity_assessment = {
                "intermediate_toxicity_scores": toxicity_result["intermediate_toxicity_scores"],
                "overall_toxicity_risk": toxicity_result["overall_toxicity_risk"],
                "flagged_intermediates": toxicity_result["flagged_intermediates"]
            }
            
            logger.info(f"Toxicity assessment: risk={toxicity_result['overall_toxicity_risk']}, "
                       f"flagged={len(toxicity_result['flagged_intermediates'])}")
            
        except Exception as e:
            logger.error(f"Toxicity assessment failed: {type(e).__name__}: {e}")
            toxicity_assessment = {
                "intermediate_toxicity_scores": {},
                "overall_toxicity_risk": "MEDIUM",
                "flagged_intermediates": []
            }
        
        # Construct Stage 3 output
        stage_3_output = {
            "pipeline_id": pipeline_id,
            "stage_2_output": stage_2_json,
            "fba_results": fba_results_dict,
            "strain_design": strain_design_dict,
            "toxicity_assessment": toxicity_assessment,
            "stage_3_status": "PASS" if fba_success else "WARN"
        }
        
        # Validate output schema
        logger.info("Validating Stage 3 output schema")
        output_validation = validate_stage_output(stage_3_output, "stage_3_output")
        if not output_validation["valid"]:
            logger.warning(f"Output schema validation issues: {output_validation.get('errors', [])}")
            if stage_3_output["stage_3_status"] == "PASS":
                stage_3_output["stage_3_status"] = "WARN"
        
        # Log output summary
        duration = (datetime.now() - start_time).total_seconds()
        logger.info(f"Stage 3 output JSON: pipeline_id={pipeline_id}, status={stage_3_output['stage_3_status']}")
        logger.info(f"=== STAGE 3 COMPLETE === duration={duration:.2f}s status={stage_3_output['stage_3_status']}")
        
        # Write stage summary
        _write_stage_summary(stage_3_output, pipeline_id, duration)
        
        return stage_3_output
        
    except Exception as e:
        logger.error(f"Stage 3 critical failure: {type(e).__name__}: {e}")
        logger.error(f"Input JSON that caused error: {json.dumps(stage_2_json, default=str)[:1000]}")
        
        # Return error output
        error_output = _create_error_stage3_output(pipeline_id, str(e))
        _write_stage_summary(error_output, pipeline_id, (datetime.now() - start_time).total_seconds())
        return error_output


def _create_fallback_fba_results(organism: str, target: str) -> Dict[str, Any]:
    """Create fallback FBA results when FBA fails."""
    
    # Biologically plausible fallback values
    growth_rates = {
        "ecoli": 0.7, "scerevisiae": 0.35, "bsubtilis": 0.5,
        "cglutamicum": 0.4, "pputida": 0.45
    }
    
    base_fluxes = {
        "lycopene": 0.5, "vanillin": 0.3, "l_lysine": 1.0,
        "riboflavin": 0.6, "pha": 0.8
    }
    
    growth = growth_rates.get(organism, 0.5)
    product_flux = base_fluxes.get(target.lower(), 0.3)
    substrate_uptake = 5.0
    
    return {
        "objective_value": growth,
        "growth_rate_per_hour": growth,
        "product_flux_mmol_per_gdw_per_hour": product_flux,
        "substrate_uptake_mmol_per_gdw_per_hour": substrate_uptake,
        "theoretical_max_yield": product_flux / substrate_uptake if substrate_uptake > 0 else 0.0,
        "flux_map": {
            "EX_glucose": -substrate_uptake,
            "BIOMASS": growth,
            "PRODUCT": product_flux
        },
        "success": False,
        "status": "FALLBACK",
        "message": "FBA failed, using estimated values"
    }


def _create_fallback_strain_design(organism: str, target: str) -> Dict[str, Any]:
    """Create fallback strain design when optimization fails."""
    
    default_knockouts = {
        "lycopene": ["ldhA", "adhE"],
        "vanillin": ["fadR"],
        "l_lysine": ["ldhA"],
        "riboflavin": ["ribR"],
        "pha": ["fadR"]
    }
    
    default_overexpressions = {
        "lycopene": ["dxs", "crtE"],
        "vanillin": ["aroG"],
        "l_lysine": ["dapA"],
        "riboflavin": ["ribA"],
        "pha": ["phaA"]
    }
    
    return {
        "algorithm_used": "Rule-based fallback",
        "final_knockouts": default_knockouts.get(target.lower(), ["ldhA"]),
        "final_overexpressions": default_overexpressions.get(target.lower(), ["dxs"]),
        "predicted_titer_g_per_l": 1.0,
        "predicted_productivity_g_per_l_per_h": 0.04,
        "metabolic_burden_score": 0.3,
        "genetic_stability_score": 0.85
    }


def _create_error_stage3_output(pipeline_id: str, error_message: str) -> Dict[str, Any]:
    """Create error Stage 3 output."""
    return {
        "pipeline_id": pipeline_id,
        "stage_2_output": {},
        "fba_results": {
            "objective_value": 0.0,
            "growth_rate_per_hour": 0.0,
            "product_flux_mmol_per_gdw_per_hour": 0.0,
            "substrate_uptake_mmol_per_gdw_per_hour": 0.0,
            "theoretical_max_yield": 0.0,
            "flux_map": {},
            "success": False,
            "status": "ERROR",
            "message": error_message
        },
        "strain_design": {
            "algorithm_used": "None",
            "final_knockouts": [],
            "final_overexpressions": [],
            "predicted_titer_g_per_l": 0.0,
            "predicted_productivity_g_per_l_per_h": 0.0,
            "metabolic_burden_score": 0.0,
            "genetic_stability_score": 0.0
        },
        "toxicity_assessment": {
            "intermediate_toxicity_scores": {},
            "overall_toxicity_risk": "UNKNOWN",
            "flagged_intermediates": []
        },
        "stage_3_status": "FAIL"
    }


def _write_stage_summary(output: Dict[str, Any], pipeline_id: str, duration: float) -> None:
    """Write stage summary to logs directory."""
    import os
    from pathlib import Path
    
    logs_dir = Path("logs")
    logs_dir.mkdir(exist_ok=True)
    
    summary = {
        "pipeline_id": pipeline_id,
        "stage": 3,
        "duration_seconds": duration,
        "status": output.get("stage_3_status", "UNKNOWN"),
        "fba_success": output.get("fba_results", {}).get("success", False),
        "growth_rate": output.get("fba_results", {}).get("growth_rate_per_hour", 0.0),
        "predicted_titer": output.get("strain_design", {}).get("predicted_titer_g_per_l", 0.0),
        "toxicity_risk": output.get("toxicity_assessment", {}).get("overall_toxicity_risk", "UNKNOWN"),
        "timestamp": datetime.now().isoformat()
    }
    
    summary_file = logs_dir / "stage_3_summary.json"
    with open(summary_file, "w") as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"Stage 3 summary written to {summary_file}")


if __name__ == "__main__":
    print("=== Testing Stage 3 Orchestrator ===")
    
    # Create mock Stage 2 output
    mock_stage_2 = {
        "pipeline_id": "test-pipeline-001",
        "stage_1_output": {
            "organism": {
                "name": "ecoli",
                "strain": "K-12 MG1655",
                "doubling_time_min": 20.0
            },
            "target_molecule": {
                "name": "lycopene",
                "smiles": "CC(C=CCC=C(C)C=CC=C(C)C=CC=C(C)C)C",
                "target_titer_g_per_l": 5.0
            }
        },
        "pathway_candidates": [
            {
                "pathway_id": "LYC_PATHWAY_001",
                "rank": 1,
                "pathway_name": "MEP pathway to lycopene",
                "steps": [
                    {"step_number": 1, "reaction_id": "R1", "substrate": "glucose", "product": "g6p"},
                    {"step_number": 2, "reaction_id": "R2", "substrate": "g6p", "product": "pyruvate"},
                    {"step_number": 3, "reaction_id": "R3", "substrate": "pyruvate", "product": "dxp"},
                    {"step_number": 4, "reaction_id": "R4", "substrate": "dxp", "product": "ipp"},
                    {"step_number": 5, "reaction_id": "R5", "substrate": "ipp", "product": "lycopene"}
                ],
                "total_steps": 5,
                "predicted_yield_mol_per_mol": 0.3
            }
        ],
        "gene_modifications": {
            "knockouts": ["ldhA", "adhE"],
            "overexpressions": ["dxs", "idi"],
            "heterologous_insertions": ["crtE", "crtB", "crtI"]
        },
        "codon_optimized_sequences": {
            "dxs": "ATGAGTCGT...",
            "crtE": "ATGCCGAAA...",
            "crtB": "ATGGCTTAC...",
            "crtI": "ATGTTCAAG..."
        },
        "stage_2_status": "PASS"
    }
    
    # Run Stage 3
    result = run_stage_3(mock_stage_2)
    
    print(f"\nStage 3 Result:")
    print(f"  Status: {result['stage_3_status']}")
    print(f"  Growth rate: {result['fba_results']['growth_rate_per_hour']:.4f}/h")
    print(f"  Predicted titer: {result['strain_design']['predicted_titer_g_per_l']:.2f} g/L")
    print(f"  Toxicity risk: {result['toxicity_assessment']['overall_toxicity_risk']}")
    print(f"  Knockouts: {result['strain_design']['final_knockouts']}")
    print(f"  Overexpressions: {result['strain_design']['final_overexpressions']}")
    
    print("\n=== Stage 3 Orchestrator Test Complete ===")
