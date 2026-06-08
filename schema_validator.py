"""
Schema Validator Module
Defines JSON schemas for inter-stage data contracts and validation functions.
"""

import json
from typing import Any, Dict, List, Optional
from datetime import datetime

try:
    import jsonschema
    from jsonschema import validate, ValidationError as JsonSchemaValidationError
    JSONSCHEMA_AVAILABLE = True
except ImportError:
    JSONSCHEMA_AVAILABLE = False
    print("WARNING: jsonschema not installed. Validation will use basic type checking.")


# ============================================================================
# STAGE 1 OUTPUT SCHEMA (Input to Stage 2)
# ============================================================================

STAGE_1_OUTPUT_SCHEMA = {
    "$schema": "http://json-schema.org/draft-07/schema#",
    "$id": "pipeline_stage_1_output_v1",
    "title": "Stage 1 Output - Core Infrastructure Data",
    "description": "Output from Stage 1 containing organism and molecule data",
    "type": "object",
    "required": [
        "pipeline_id", "timestamp", "organism", "target_molecule",
        "genomic_data", "data_quality_report", "stage_1_status"
    ],
    "properties": {
        "pipeline_id": {"type": "string", "format": "uuid"},
        "timestamp": {"type": "string", "format": "date-time"},
        "organism": {
            "type": "object",
            "required": ["name", "strain", "gem_model_id", "doubling_time_min",
                        "optimal_ph", "optimal_temp_c", "gram_stain"],
            "properties": {
                "name": {"type": "string"},
                "strain": {"type": "string"},
                "gem_model_id": {"type": "string"},
                "doubling_time_min": {"type": "number", "minimum": 0},
                "optimal_ph": {"type": "number", "minimum": 0, "maximum": 14},
                "optimal_temp_c": {"type": "number", "minimum": -20, "maximum": 100},
                "gram_stain": {"type": "string", "enum": ["positive", "negative"]},
                "codon_table_id": {"type": "integer"},
                "gc_content_percent": {"type": "number", "minimum": 0, "maximum": 100}
            }
        },
        "target_molecule": {
            "type": "object",
            "required": ["name", "smiles", "chebi_id", "target_titer_g_per_l", "target_yield_mol_per_mol"],
            "properties": {
                "name": {"type": "string"},
                "smiles": {"type": "string"},
                "chebi_id": {"type": "string"},
                "target_titer_g_per_l": {"type": "number", "minimum": 0},
                "target_yield_mol_per_mol": {"type": "number", "minimum": 0, "maximum": 1},
                "molecular_weight": {"type": "number"},
                "category": {"type": "string"},
                "toxicity_level": {"type": "string", "enum": ["LOW", "MEDIUM", "HIGH"]},
                "secretion_type": {"type": "string", "enum": ["intracellular", "secreted", "periplasmic"]}
            }
        },
        "genomic_data": {
            "type": "object",
            "required": ["total_genes", "essential_genes", "available_promoters", "codon_table_id", "gc_content_percent"],
            "properties": {
                "total_genes": {"type": "integer", "minimum": 0},
                "essential_genes": {"type": "array", "items": {"type": "string"}},
                "available_promoters": {"type": "array", "items": {"type": "string"}},
                "codon_table_id": {"type": "integer"},
                "gc_content_percent": {"type": "number", "minimum": 0, "maximum": 100}
            }
        },
        "data_quality_report": {
            "type": "object",
            "required": ["completeness_score", "warnings", "errors"],
            "properties": {
                "completeness_score": {"type": "number", "minimum": 0, "maximum": 1},
                "warnings": {"type": "array", "items": {"type": "string"}},
                "errors": {"type": "array", "items": {"type": "string"}}
            }
        },
        "stage_1_status": {"type": "string", "enum": ["PASS", "FAIL", "WARN"]}
    }
}


# ============================================================================
# STAGE 2 OUTPUT SCHEMA (Input to Stage 3)
# ============================================================================

STAGE_2_OUTPUT_SCHEMA = {
    "$schema": "http://json-schema.org/draft-07/schema#",
    "$id": "pipeline_stage_2_output_v1",
    "title": "Stage 2 Output - Pathway Prediction Results",
    "description": "Output from Stage 2 containing pathway candidates and gene modifications",
    "type": "object",
    "required": [
        "pipeline_id", "stage_1_output", "pathway_candidates",
        "gene_modifications", "codon_optimized_sequences", "stage_2_status"
    ],
    "properties": {
        "pipeline_id": {"type": "string", "format": "uuid"},
        "stage_1_output": {"$ref": "#/definitions/stage_1_output"},
        "pathway_candidates": {
            "type": "array",
            "minItems": 1,
            "items": {"$ref": "#/definitions/pathway_candidate"}
        },
        "gene_modifications": {
            "type": "object",
            "required": ["knockouts", "overexpressions", "heterologous_insertions"],
            "properties": {
                "knockouts": {"type": "array", "items": {"type": "string"}},
                "overexpressions": {"type": "array", "items": {"type": "string"}},
                "heterologous_insertions": {"type": "array", "items": {"type": "string"}}
            }
        },
        "codon_optimized_sequences": {
            "type": "object",
            "additionalProperties": {"type": "string"}
        },
        "stage_2_status": {"type": "string", "enum": ["PASS", "FAIL", "WARN"]}
    },
    "definitions": {
        "stage_1_output": STAGE_1_OUTPUT_SCHEMA,
        "pathway_candidate": {
            "type": "object",
            "required": [
                "pathway_id", "rank", "pathway_name", "steps",
                "total_steps", "predicted_yield_mol_per_mol",
                "thermodynamic_feasibility_score", "gnn_viability_score",
                "host_compatibility_score"
            ],
            "properties": {
                "pathway_id": {"type": "string"},
                "rank": {"type": "integer", "minimum": 1},
                "pathway_name": {"type": "string"},
                "steps": {
                    "type": "array",
                    "items": {"$ref": "#/definitions/pathway_step"}
                },
                "total_steps": {"type": "integer", "minimum": 1},
                "predicted_yield_mol_per_mol": {"type": "number", "minimum": 0, "maximum": 1},
                "thermodynamic_feasibility_score": {"type": "number", "minimum": 0, "maximum": 1},
                "gnn_viability_score": {"type": "number", "minimum": 0, "maximum": 1},
                "host_compatibility_score": {"type": "number", "minimum": 0, "maximum": 1}
            }
        },
        "pathway_step": {
            "type": "object",
            "required": [
                "step_number", "reaction_id", "enzyme_name", "gene_name",
                "ec_number", "substrate", "product", "delta_g_kj_per_mol",
                "is_heterologous"
            ],
            "properties": {
                "step_number": {"type": "integer", "minimum": 1},
                "reaction_id": {"type": "string"},
                "enzyme_name": {"type": "string"},
                "gene_name": {"type": "string"},
                "ec_number": {"type": "string"},
                "substrate": {"type": "string"},
                "product": {"type": "string"},
                "delta_g_kj_per_mol": {"type": "number"},
                "kcat_per_sec": {"type": "number", "minimum": 0},
                "km_mm": {"type": "number", "minimum": 0},
                "is_heterologous": {"type": "boolean"},
                "source_organism": {"type": "string"}
            }
        }
    }
}


# ============================================================================
# STAGE 3 OUTPUT SCHEMA (Input to Stage 4)
# ============================================================================

STAGE_3_OUTPUT_SCHEMA = {
    "$schema": "http://json-schema.org/draft-07/schema#",
    "$id": "pipeline_stage_3_output_v1",
    "title": "Stage 3 Output - Flux Analysis and Strain Design",
    "description": "Output from Stage 3 containing FBA results and strain optimization",
    "type": "object",
    "required": [
        "pipeline_id", "stage_2_output", "fba_results", "strain_design",
        "toxicity_assessment", "stage_3_status"
    ],
    "properties": {
        "pipeline_id": {"type": "string", "format": "uuid"},
        "stage_2_output": {"$ref": "#/definitions/stage_2_output"},
        "fba_results": {"$ref": "#/definitions/fba_results"},
        "strain_design": {"$ref": "#/definitions/strain_design"},
        "toxicity_assessment": {"$ref": "#/definitions/toxicity_assessment"},
        "stage_3_status": {"type": "string", "enum": ["PASS", "FAIL", "WARN"]}
    },
    "definitions": {
        "stage_2_output": STAGE_2_OUTPUT_SCHEMA,
        "fba_results": {
            "type": "object",
            "required": [
                "objective_value", "growth_rate_per_hour",
                "product_flux_mmol_per_gdw_per_hour", "substrate_uptake_mmol_per_gdw_per_hour",
                "theoretical_max_yield", "flux_map"
            ],
            "properties": {
                "objective_value": {"type": "number"},
                "growth_rate_per_hour": {"type": "number", "minimum": 0},
                "product_flux_mmol_per_gdw_per_hour": {"type": "number"},
                "substrate_uptake_mmol_per_gdw_per_hour": {"type": "number"},
                "theoretical_max_yield": {"type": "number", "minimum": 0, "maximum": 1},
                "flux_map": {"type": "object", "additionalProperties": {"type": "number"}}
            }
        },
        "strain_design": {
            "type": "object",
            "required": [
                "algorithm_used", "final_knockouts", "final_overexpressions",
                "predicted_titer_g_per_l", "predicted_productivity_g_per_l_per_h",
                "metabolic_burden_score", "genetic_stability_score"
            ],
            "properties": {
                "algorithm_used": {"type": "string"},
                "final_knockouts": {"type": "array", "items": {"type": "string"}},
                "final_overexpressions": {"type": "array", "items": {"type": "string"}},
                "predicted_titer_g_per_l": {"type": "number", "minimum": 0},
                "predicted_productivity_g_per_l_per_h": {"type": "number", "minimum": 0},
                "metabolic_burden_score": {"type": "number", "minimum": 0, "maximum": 1},
                "genetic_stability_score": {"type": "number", "minimum": 0, "maximum": 1}
            }
        },
        "toxicity_assessment": {
            "type": "object",
            "required": ["intermediate_toxicity_scores", "overall_toxicity_risk", "flagged_intermediates"],
            "properties": {
                "intermediate_toxicity_scores": {
                    "type": "object",
                    "additionalProperties": {"type": "number", "minimum": 0, "maximum": 1}
                },
                "overall_toxicity_risk": {"type": "string", "enum": ["LOW", "MEDIUM", "HIGH"]},
                "flagged_intermediates": {"type": "array", "items": {"type": "string"}}
            }
        }
    }
}


# ============================================================================
# STAGE 4 OUTPUT SCHEMA (Input to Stage 5)
# ============================================================================

STAGE_4_OUTPUT_SCHEMA = {
    "$schema": "http://json-schema.org/draft-07/schema#",
    "$id": "pipeline_stage_4_output_v1",
    "title": "Stage 4 Output - DBTL Loop and Fermentation Results",
    "description": "Output from Stage 4 containing DBTL cycles and fermentation simulation",
    "type": "object",
    "required": [
        "pipeline_id", "stage_3_output", "dbtl_cycles", "fermentation_simulation",
        "optimal_fermentation_conditions", "stage_4_status"
    ],
    "properties": {
        "pipeline_id": {"type": "string", "format": "uuid"},
        "stage_3_output": {"$ref": "#/definitions/stage_3_output"},
        "dbtl_cycles": {
            "type": "array",
            "minItems": 1,
            "items": {"$ref": "#/definitions/dbtl_cycle"}
        },
        "fermentation_simulation": {"$ref": "#/definitions/fermentation_simulation"},
        "optimal_fermentation_conditions": {"$ref": "#/definitions/fermentation_conditions"},
        "stage_4_status": {"type": "string", "enum": ["PASS", "FAIL", "WARN"]}
    },
    "definitions": {
        "stage_3_output": STAGE_3_OUTPUT_SCHEMA,
        "dbtl_cycle": {
            "type": "object",
            "required": [
                "cycle_number", "constructs_tested", "best_titer_g_per_l",
                "best_construct_id", "improvement_fold", "bo_next_candidates"
            ],
            "properties": {
                "cycle_number": {"type": "integer", "minimum": 1},
                "constructs_tested": {"type": "integer", "minimum": 1},
                "best_titer_g_per_l": {"type": "number", "minimum": 0},
                "best_construct_id": {"type": "string"},
                "improvement_fold": {"type": "number", "minimum": 0},
                "bo_next_candidates": {"type": "array", "items": {"type": "string"}}
            }
        },
        "fermentation_simulation": {
            "type": "object",
            "required": [
                "mode", "duration_hours", "final_titer_g_per_l",
                "final_yield_g_per_g", "final_productivity_g_per_l_per_h",
                "ode_convergence", "organism_specific_events"
            ],
            "properties": {
                "mode": {"type": "string", "enum": ["batch", "fed-batch", "continuous"]},
                "duration_hours": {"type": "number", "minimum": 0},
                "final_titer_g_per_l": {"type": "number", "minimum": 0},
                "final_yield_g_per_g": {"type": "number", "minimum": 0},
                "final_productivity_g_per_l_per_h": {"type": "number", "minimum": 0},
                "ode_convergence": {"type": "boolean"},
                "organism_specific_events": {"type": "array", "items": {"type": "string"}}
            }
        },
        "fermentation_conditions": {
            "type": "object",
            "required": [
                "temperature_c", "ph", "do_percent_saturation",
                "glucose_feed_g_per_l_per_h", "agitation_rpm", "aeration_vvm"
            ],
            "properties": {
                "temperature_c": {"type": "number", "minimum": 0, "maximum": 100},
                "ph": {"type": "number", "minimum": 0, "maximum": 14},
                "do_percent_saturation": {"type": "number", "minimum": 0, "maximum": 100},
                "glucose_feed_g_per_l_per_h": {"type": "number", "minimum": 0},
                "agitation_rpm": {"type": "number", "minimum": 0},
                "aeration_vvm": {"type": "number", "minimum": 0}
            }
        }
    }
}


# ============================================================================
# VALIDATION FUNCTIONS
# ============================================================================

SCHEMA_VERSION = "1.0.0"
SCHEMA_REGISTRY = {
    "stage_1_output": STAGE_1_OUTPUT_SCHEMA,
    "stage_2_output": STAGE_2_OUTPUT_SCHEMA,
    "stage_3_output": STAGE_3_OUTPUT_SCHEMA,
    "stage_4_output": STAGE_4_OUTPUT_SCHEMA
}


def validate_stage_output(
    data: Dict[str, Any],
    schema_name: str,
    raise_on_error: bool = True
) -> tuple[bool, List[str]]:
    """
    Validate output data against a predefined JSON schema.
    
    Args:
        data: The data dictionary to validate
        schema_name: Name of the schema to validate against
                    (e.g., "stage_1_output", "stage_2_output")
        raise_on_error: If True, raise SchemaValidationError on failure
    
    Returns:
        Tuple of (is_valid: bool, error_messages: List[str])
    
    Raises:
        SchemaValidationError: If validation fails and raise_on_error=True
    """
    from exceptions import SchemaValidationError
    
    errors = []
    
    # Check if schema exists
    if schema_name not in SCHEMA_REGISTRY:
        error_msg = f"Unknown schema: {schema_name}. Available: {list(SCHEMA_REGISTRY.keys())}"
        errors.append(error_msg)
        if raise_on_error:
            raise SchemaValidationError(
                message=error_msg,
                schema_name=schema_name
            )
        return False, errors
    
    schema = SCHEMA_REGISTRY[schema_name]
    
    if JSONSCHEMA_AVAILABLE:
        try:
            validate(instance=data, schema=schema)
            return True, []
        except JsonSchemaValidationError as e:
            error_msg = f"Schema validation failed for {schema_name}: {e.message}"
            errors.append(error_msg)
            
            if raise_on_error:
                raise SchemaValidationError(
                    message=error_msg,
                    schema_name=schema_name,
                    validation_errors=[e],
                    input_json=data
                )
    else:
        # Fallback: basic type checking
        warnings = ["jsonschema not available, using basic validation"]
        if not isinstance(data, dict):
            errors.append(f"Expected dict, got {type(data).__name__}")
        
        required_fields = schema.get("required", [])
        for field in required_fields:
            if field not in data:
                errors.append(f"Missing required field: {field}")
        
        if errors and raise_on_error:
            raise SchemaValidationError(
                message="Basic validation failed",
                schema_name=schema_name,
                validation_errors=errors,
                input_json=data
            )
    
    return len(errors) == 0, errors


def get_schema_version() -> str:
    """Return the current schema version."""
    return SCHEMA_VERSION


def log_json_contract(
    data: Dict[str, Any],
    schema_name: str,
    logger: Any,
    stage: int
) -> None:
    """
    Log JSON contract details for debugging.
    
    Args:
        data: The data being validated
        schema_name: Name of the schema
        logger: Logger instance
        stage: Pipeline stage number
    """
    import json
    
    logger.info(f"Validating against schema: {schema_name} (v{SCHEMA_VERSION})")
    
    # Log summary (first 500 chars of JSON)
    json_str = json.dumps(data, indent=2, default=str)
    summary = json_str[:500] + "..." if len(json_str) > 500 else json_str
    logger.debug(f"JSON payload preview: {summary}")
    
    # Log key fields based on schema type
    if "pipeline_id" in data:
        logger.debug(f"Pipeline ID: {data['pipeline_id']}")
    if "stage_1_status" in data:
        logger.debug(f"Stage status: {data['stage_1_status']}")
    if "pathway_candidates" in data:
        num_pathways = len(data.get("pathway_candidates", []))
        logger.debug(f"Number of pathway candidates: {num_pathways}")
    if "fba_results" in data:
        fba = data.get("fba_results", {})
        logger.debug(f"FBA objective: {fba.get('objective_value', 'N/A')}")


def create_sample_pipeline_id() -> str:
    """Create a sample pipeline ID (UUID format)."""
    import uuid
    return str(uuid.uuid4())


def create_timestamp() -> str:
    """Create an ISO8601 timestamp."""
    return datetime.utcnow().isoformat() + "Z"


if __name__ == "__main__":
    # Test schema validation
    print("Testing Schema Validator...")
    print(f"Schema version: {get_schema_version()}")
    print(f"Available schemas: {list(SCHEMA_REGISTRY.keys())}")
    print(f"jsonschema available: {JSONSCHEMA_AVAILABLE}")
    
    # Test with sample Stage 1 output
    sample_stage_1 = {
        "pipeline_id": create_sample_pipeline_id(),
        "timestamp": create_timestamp(),
        "organism": {
            "name": "Escherichia coli",
            "strain": "K-12 MG1655",
            "gem_model_id": "iML1515",
            "doubling_time_min": 20.0,
            "optimal_ph": 7.0,
            "optimal_temp_c": 37.0,
            "gram_stain": "negative",
            "codon_table_id": 11,
            "gc_content_percent": 50.8
        },
        "target_molecule": {
            "name": "lycopene",
            "smiles": "CC1=CCC2(C1)CC=C3C(=CC=C4C(=C3)C(=CC=C5C(=C4)CC=C6C(=C5)CC(=CC=C6)C)C)C2",
            "chebi_id": "CHEBI:15956",
            "target_titer_g_per_l": 5.0,
            "target_yield_mol_per_mol": 0.3,
            "molecular_weight": 536.87,
            "category": "carotenoid",
            "toxicity_level": "LOW",
            "secretion_type": "intracellular"
        },
        "genomic_data": {
            "total_genes": 4377,
            "essential_genes": ["dnaA", "dnaN", "gyrA"],
            "available_promoters": ["J23100", "Ptac", "PBAD"],
            "codon_table_id": 11,
            "gc_content_percent": 50.8
        },
        "data_quality_report": {
            "completeness_score": 0.95,
            "warnings": [],
            "errors": []
        },
        "stage_1_status": "PASS"
    }
    
    print("\n--- Testing Stage 1 Schema Validation ---")
    is_valid, errors = validate_stage_output(sample_stage_1, "stage_1_output", raise_on_error=False)
    print(f"Validation result: {'PASSED' if is_valid else 'FAILED'}")
    if errors:
        print(f"Errors: {errors}")
    else:
        print("✓ No validation errors!")
    
    # Test invalid data
    print("\n--- Testing Invalid Data ---")
    invalid_data = {"pipeline_id": "test"}  # Missing required fields
    is_valid, errors = validate_stage_output(invalid_data, "stage_1_output", raise_on_error=False)
    print(f"Validation result: {'PASSED' if is_valid else 'FAILED'}")
    print(f"Expected errors detected: {len(errors)}")
    
    print("\n✓ Schema validator tests complete!")
