"""
FBA Engine - Flux Balance Analysis for metabolic pathway optimization.

This module provides flux balance analysis capabilities using scipy.optimize
without requiring COBRApy dependency.
"""

import logging
import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any
from scipy.optimize import linprog
import json

from logger_setup import setup_logger, PipelineLogger
from exceptions import PipelineError, FBAConvergenceError
from pipeline_config import PipelineConfig

logger = setup_logger("STAGE_3", "fba_engine")


@dataclass
class FBAModel:
    """Dataclass representing a Flux Balance Analysis model."""
    
    model_id: str
    organism: str
    reactions: List[str] = field(default_factory=list)
    metabolites: List[str] = field(default_factory=list)
    stoichiometric_matrix: np.ndarray = field(default_factory=lambda: np.array([]))
    lower_bounds: np.ndarray = field(default_factory=lambda: np.array([]))
    upper_bounds: np.ndarray = field(default_factory=lambda: np.array([]))
    objective_coefficients: np.ndarray = field(default_factory=lambda: np.array([]))
    biomass_reaction_id: str = ""
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert FBA model to dictionary for JSON serialization."""
        return {
            "model_id": self.model_id,
            "organism": self.organism,
            "num_reactions": len(self.reactions),
            "num_metabolites": len(self.metabolites),
            "biomass_reaction_id": self.biomass_reaction_id,
            "stoichiometric_matrix_shape": list(self.stoichiometric_matrix.shape) if self.stoichiometric_matrix.size > 0 else [0, 0]
        }


@dataclass
class FBAResult:
    """Dataclass representing FBA computation results."""
    
    success: bool
    objective_value: float
    growth_rate_per_hour: float
    product_flux_mmol_per_gdw_per_hour: float
    substrate_uptake_mmol_per_gdw_per_hour: float
    theoretical_max_yield: float
    flux_map: Dict[str, float]
    status: str
    message: str
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert FBA result to dictionary."""
        return {
            "success": self.success,
            "objective_value": self.objective_value,
            "growth_rate_per_hour": self.growth_rate_per_hour,
            "product_flux_mmol_per_gdw_per_hour": self.product_flux_mmol_per_gdw_per_hour,
            "substrate_uptake_mmol_per_gdw_per_hour": self.substrate_uptake_mmol_per_gdw_per_hour,
            "theoretical_max_yield": self.theoretical_max_yield,
            "flux_map": self.flux_map,
            "status": self.status,
            "message": self.message
        }


class FBAEngine:
    """
    Flux Balance Analysis Engine for metabolic network analysis.
    
    Provides methods for FBA, pFBA, and FVA using scipy.optimize.linprog
    without external COBRApy dependency.
    """
    
    # Organism-specific growth parameters
    ORGANISM_PARAMS = {
        "ecoli": {
            "mu_max": 0.9,  # 1/h
            "glucose_uptake_max": 10.0,  # mmol/gDW/h
            "atp_maintenance": 8.39,  # mmol/gDW/h
            "biomass_coefficient": 1.0
        },
        "scerevisiae": {
            "mu_max": 0.45,
            "glucose_uptake_max": 8.0,
            "atp_maintenance": 6.0,
            "biomass_coefficient": 1.0
        },
        "bsubtilis": {
            "mu_max": 0.7,
            "glucose_uptake_max": 7.5,
            "atp_maintenance": 7.0,
            "biomass_coefficient": 1.0
        },
        "cglutamicum": {
            "mu_max": 0.55,
            "glucose_uptake_max": 6.5,
            "atp_maintenance": 5.5,
            "biomass_coefficient": 1.0
        },
        "pputida": {
            "mu_max": 0.6,
            "glucose_uptake_max": 7.0,
            "atp_maintenance": 6.5,
            "biomass_coefficient": 1.0
        }
    }
    
    def __init__(self, config: Optional[PipelineConfig] = None):
        """
        Initialize FBA Engine.
        
        Args:
            config: Pipeline configuration object
        """
        self.config = config or PipelineConfig()
        self.logger = setup_logger("STAGE_3", "fba_engine")
        self.current_model: Optional[FBAModel] = None
        
    def build_stoichiometric_matrix(
        self, 
        pathway_steps: List[Dict[str, Any]],
        organism: str
    ) -> FBAModel:
        """
        Build stoichiometric matrix from pathway steps.
        
        Args:
            pathway_steps: List of pathway step dictionaries with substrate/product info
            organism: Organism identifier
            
        Returns:
            FBAModel with constructed stoichiometric matrix
        """
        self.logger.info(f"Building stoichiometric matrix for {len(pathway_steps)} pathway steps")
        self.logger.debug(f"Pathway steps: {pathway_steps[:3]}...")
        
        # Extract unique metabolites and reactions
        metabolites_set = set()
        reactions_list = []
        
        # Add glucose uptake and biomass reactions
        reactions_list = ["EX_glucose", "BIOMASS"]
        metabolites_set.update(["glucose", "biomass"])
        
        for step in pathway_steps:
            reaction_id = step.get("reaction_id", f"R{step.get('step_number', 0)}")
            substrate = step.get("substrate", "")
            product = step.get("product", "")
            
            reactions_list.append(reaction_id)
            if substrate:
                metabolites_set.add(substrate)
            if product:
                metabolites_set.add(product)
        
        # Add ATP/ADP/NADH/NADPH balancing
        metabolites_set.update(["atp", "adp", "nadph", "nadp", "pi"])
        
        metabolites_list = sorted(list(metabolites_set))
        num_metabolites = len(metabolites_list)
        num_reactions = len(reactions_list)
        
        # Initialize stoichiometric matrix (metabolites x reactions)
        S = np.zeros((num_metabolites, num_reactions))
        
        # Create index mappings
        met_index = {m: i for i, m in enumerate(metabolites_list)}
        rxn_index = {r: i for i, r in enumerate(reactions_list)}
        
        # Add glucose uptake (negative = consumption)
        if "glucose" in met_index:
            S[met_index["glucose"], rxn_index["EX_glucose"]] = -1.0
        
        # Add biomass formation
        if "biomass" in met_index:
            S[met_index["biomass"], rxn_index["BIOMASS"]] = 1.0
        # Biomass consumes precursors (simplified)
        if "glucose" in met_index:
            S[met_index["glucose"], rxn_index["BIOMASS"]] = -0.5
        if "atp" in met_index:
            S[met_index["atp"], rxn_index["BIOMASS"]] = -30.0
            S[met_index["adp"], rxn_index["BIOMASS"]] = 30.0
        
        # Add pathway reactions
        for step in pathway_steps:
            reaction_id = step.get("reaction_id", f"R{step.get('step_number', 0)}")
            substrate = step.get("substrate", "")
            product = step.get("product", "")
            
            if reaction_id in rxn_index:
                rxn_idx = rxn_index[reaction_id]
                if substrate and substrate in met_index:
                    S[met_index[substrate], rxn_idx] = -1.0
                if product and product in met_index:
                    S[met_index[product], rxn_idx] = 1.0
        
        # Set bounds
        lb = np.full(num_reactions, -1000.0)
        ub = np.full(num_reactions, 1000.0)
        
        # Glucose uptake constraint
        org_params = self.ORGANISM_PARAMS.get(organism, self.ORGANISM_PARAMS["ecoli"])
        lb[rxn_index["EX_glucose"]] = -org_params["glucose_uptake_max"]
        ub[rxn_index["EX_glucose"]] = 0.0
        
        # Biomass must be positive
        lb[rxn_index["BIOMASS"]] = 0.0
        
        # Objective: maximize biomass
        obj = np.zeros(num_reactions)
        obj[rxn_index["BIOMASS"]] = -1.0  # Negative for maximization
        
        model = FBAModel(
            model_id=f"FBA_{organism}_{len(pathway_steps)}steps",
            organism=organism,
            reactions=reactions_list,
            metabolites=metabolites_list,
            stoichiometric_matrix=S,
            lower_bounds=lb,
            upper_bounds=ub,
            objective_coefficients=obj,
            biomass_reaction_id="BIOMASS"
        )
        
        self.current_model = model
        self.logger.info(f"Built stoichiometric matrix: {num_metabolites} metabolites x {num_reactions} reactions")
        self.logger.debug(f"Model summary: {model.to_dict()}")
        
        return model
    
    def organism_specific_constraints(
        self, 
        model: FBAModel, 
        organism: str,
        additional_constraints: Optional[Dict[str, Tuple[float, float]]] = None
    ) -> FBAModel:
        """
        Apply organism-specific constraints to the model.
        
        Args:
            model: FBA model to constrain
            organism: Organism identifier
            additional_constraints: Optional dict of reaction_id -> (lb, ub) tuples
            
        Returns:
            Constrained FBAModel
        """
        self.logger.info(f"Applying organism-specific constraints for {organism}")
        
        org_params = self.ORGANISM_PARAMS.get(organism, self.ORGANISM_PARAMS["ecoli"])
        
        # Adjust ATP maintenance based on organism
        if "atp" in model.metabolites and "BIOMASS" in model.reactions:
            met_idx = model.metabolites.index("atp")
            rxn_idx = model.reactions.index("BIOMASS")
            # Already set in build_stoichiometric_matrix
            
        # Apply additional constraints
        if additional_constraints:
            for rxn_id, (lb, ub) in additional_constraints.items():
                if rxn_id in model.reactions:
                    rxn_idx = model.reactions.index(rxn_id)
                    model.lower_bounds[rxn_idx] = max(model.lower_bounds[rxn_idx], lb)
                    model.upper_bounds[rxn_idx] = min(model.upper_bounds[rxn_idx], ub)
                    self.logger.debug(f"Constrained {rxn_id}: [{lb}, {ub}]")
        
        return model
    
    def run_fba(
        self, 
        model: Optional[FBAModel] = None,
        product_reaction_id: Optional[str] = None
    ) -> FBAResult:
        """
        Run Flux Balance Analysis.
        
        Args:
            model: FBA model (uses current_model if None)
            product_reaction_id: Optional product reaction for yield calculation
            
        Returns:
            FBAResult with flux distribution
        """
        model = model or self.current_model
        if model is None:
            raise FBAConvergenceError("No FBA model available. Call build_stoichiometric_matrix first.")
        
        self.logger.info("Running Flux Balance Analysis")
        self.logger.debug(f"Model: {model.model_id}, {len(model.reactions)} reactions")
        
        try:
            # Solve linear programming problem
            result = linprog(
                c=model.objective_coefficients,
                A_eq=model.stoichiometric_matrix,
                b_eq=np.zeros(model.stoichiometric_matrix.shape[0]),
                bounds=list(zip(model.lower_bounds, model.upper_bounds)),
                method='highs',
                options={'presolve': True}
            )
            
            if not result.success:
                self.logger.warning(f"FBA did not converge: {result.message}")
                # Try fallback with different method
                result = linprog(
                    c=model.objective_coefficients,
                    A_eq=model.stoichiometric_matrix,
                    b_eq=np.zeros(model.stoichiometric_matrix.shape[0]),
                    bounds=list(zip(model.lower_bounds, model.upper_bounds)),
                    method='revised simplex'
                )
            
            # Extract flux values
            flux_values = result.x if hasattr(result, 'x') else np.zeros(len(model.reactions))
            flux_map = {rxn: float(flux_values[i]) for i, rxn in enumerate(model.reactions)}
            
            # Calculate metrics
            biomass_idx = model.reactions.index(model.biomass_reaction_id) if model.biomass_reaction_id in model.reactions else -1
            growth_rate = abs(flux_values[biomass_idx]) if biomass_idx >= 0 else 0.0
            
            # Product flux
            product_flux = 0.0
            if product_reaction_id and product_reaction_id in flux_map:
                product_flux = abs(flux_map[product_reaction_id])
            
            # Substrate uptake
            substrate_uptake = 0.0
            if "EX_glucose" in flux_map:
                substrate_uptake = abs(flux_map["EX_glucose"])
            
            # Theoretical max yield
            theoretical_yield = product_flux / substrate_uptake if substrate_uptake > 0 else 0.0
            
            fba_result = FBAResult(
                success=result.success,
                objective_value=-result.fun if hasattr(result, 'fun') else growth_rate,
                growth_rate_per_hour=growth_rate,
                product_flux_mmol_per_gdw_per_hour=product_flux,
                substrate_uptake_mmol_per_gdw_per_hour=substrate_uptake,
                theoretical_max_yield=theoretical_yield,
                flux_map=flux_map,
                status="OPTIMAL" if result.success else "SUBOPTIMAL",
                message=result.message if hasattr(result, 'message') else "FBA completed"
            )
            
            self.logger.info(f"FBA complete: growth={growth_rate:.4f}/h, product_flux={product_flux:.4f}")
            self.logger.debug(f"Flux map (first 5): {dict(list(flux_map.items())[:5])}")
            
            return fba_result
            
        except Exception as e:
            self.logger.error(f"FBA failed: {type(e).__name__}: {e}")
            self.logger.error(f"Model causing error: {model.to_dict()}")
            raise FBAConvergenceError(f"FBA failed: {e}")
    
    def run_pfba(
        self, 
        model: Optional[FBAModel] = None,
        target_growth_fraction: float = 0.9
    ) -> FBAResult:
        """
        Run parsimonious FBA (minimize total flux while maintaining growth).
        
        Args:
            model: FBA model
            target_growth_fraction: Fraction of optimal growth to maintain
            
        Returns:
            FBAResult with minimized flux distribution
        """
        model = model or self.current_model
        if model is None:
            raise FBAConvergenceError("No FBA model available.")
        
        self.logger.info("Running parsimonious FBA")
        
        # First get optimal growth
        optimal_result = self.run_fba(model)
        target_growth = optimal_result.growth_rate_per_hour * target_growth_fraction
        
        # Add growth constraint
        if model.biomass_reaction_id in model.reactions:
            biomass_idx = model.reactions.index(model.biomass_reaction_id)
            model.lower_bounds[biomass_idx] = target_growth
        
        # Minimize sum of absolute fluxes (approximated)
        new_obj = np.ones(len(model.reactions))
        
        try:
            result = linprog(
                c=new_obj,
                A_eq=model.stoichiometric_matrix,
                b_eq=np.zeros(model.stoichiometric_matrix.shape[0]),
                bounds=list(zip(model.lower_bounds, model.upper_bounds)),
                method='highs'
            )
            
            flux_values = result.x if hasattr(result, 'x') else np.zeros(len(model.reactions))
            flux_map = {rxn: float(flux_values[i]) for i, rxn in enumerate(model.reactions)}
            
            pfba_result = FBAResult(
                success=result.success,
                objective_value=sum(abs(flux_values)),
                growth_rate_per_hour=target_growth,
                product_flux_mmol_per_gdw_per_hour=optimal_result.product_flux_mmol_per_gdw_per_hour,
                substrate_uptake_mmol_per_gdw_per_hour=optimal_result.substrate_uptake_mmol_per_gdw_per_hour,
                theoretical_max_yield=optimal_result.theoretical_max_yield,
                flux_map=flux_map,
                status="PFBA_OPTIMAL" if result.success else "PFBA_SUBOPTIMAL",
                message="pFBA completed"
            )
            
            self.logger.info(f"pFBA complete: total_flux={pfba_result.objective_value:.4f}")
            return pfba_result
            
        except Exception as e:
            self.logger.error(f"pFBA failed: {e}")
            return optimal_result  # Fallback to standard FBA
    
    def run_fva(
        self, 
        model: Optional[FBAModel] = None,
        reaction_ids: Optional[List[str]] = None
    ) -> Dict[str, Dict[str, float]]:
        """
        Run Flux Variability Analysis.
        
        Args:
            model: FBA model
            reaction_ids: Reactions to analyze (all if None)
            
        Returns:
            Dict mapping reaction_id -> {min, max} flux values
        """
        model = model or self.current_model
        if model is None:
            raise FBAConvergenceError("No FBA model available.")
        
        self.logger.info("Running Flux Variability Analysis")
        
        if reaction_ids is None:
            reaction_ids = model.reactions
        
        fva_results = {}
        
        for rxn_id in reaction_ids:
            if rxn_id not in model.reactions:
                continue
                
            rxn_idx = model.reactions.index(rxn_id)
            
            # Maximize this reaction
            obj_max = np.zeros(len(model.reactions))
            obj_max[rxn_idx] = -1.0
            
            try:
                result_max = linprog(
                    c=obj_max,
                    A_eq=model.stoichiometric_matrix,
                    b_eq=np.zeros(model.stoichiometric_matrix.shape[0]),
                    bounds=list(zip(model.lower_bounds, model.upper_bounds)),
                    method='highs'
                )
                max_val = -result_max.fun if result_max.success else 0.0
                
                # Minimize this reaction
                obj_min = np.zeros(len(model.reactions))
                obj_min[rxn_idx] = 1.0
                
                result_min = linprog(
                    c=obj_min,
                    A_eq=model.stoichiometric_matrix,
                    b_eq=np.zeros(model.stoichiometric_matrix.shape[0]),
                    bounds=list(zip(model.lower_bounds, model.upper_bounds)),
                    method='highs'
                )
                min_val = result_min.fun if result_min.success else 0.0
                
                fva_results[rxn_id] = {"min": float(min_val), "max": float(max_val)}
                
            except Exception as e:
                self.logger.warning(f"FVA failed for {rxn_id}: {e}")
                fva_results[rxn_id] = {"min": 0.0, "max": 0.0}
        
        self.logger.debug(f"FVA results for {len(fva_results)} reactions")
        return fva_results


if __name__ == "__main__":
    # Test FBA Engine
    import json
    
    print("=== Testing FBA Engine ===")
    
    engine = FBAEngine()
    
    # Create test pathway
    test_pathway = [
        {"step_number": 1, "reaction_id": "R1", "substrate": "glucose", "product": "g6p"},
        {"step_number": 2, "reaction_id": "R2", "substrate": "g6p", "product": "f6p"},
        {"step_number": 3, "reaction_id": "R3", "substrate": "f6p", "product": "pyruvate"},
        {"step_number": 4, "reaction_id": "R4", "substrate": "pyruvate", "product": "acetate"}
    ]
    
    # Build model
    model = engine.build_stoichiometric_matrix(test_pathway, "ecoli")
    print(f"Model built: {model.to_dict()}")
    
    # Run FBA
    result = engine.run_fba(model, product_reaction_id="R4")
    print(f"\nFBA Result:")
    print(json.dumps(result.to_dict(), indent=2))
    
    # Run pFBA
    pfba_result = engine.run_pfba(model)
    print(f"\npFBA total flux: {pfba_result.objective_value:.4f}")
    
    # Run FVA
    fva_results = engine.run_fva(model, ["R1", "R2", "R3", "R4"])
    print(f"\nFVA Results:")
    for rxn, bounds in fva_results.items():
        print(f"  {rxn}: [{bounds['min']:.4f}, {bounds['max']:.4f}]")
    
    print("\n=== FBA Engine Test Complete ===")
