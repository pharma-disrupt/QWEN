"""
STAGE 2: Pathway Prediction + AI Engine

This module provides pathway prediction, enzyme selection, codon optimization,
and expression cassette design for metabolic engineering.
"""

import logging
import random
import math
import json
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Optional, Tuple, Any
from datetime import datetime
import uuid

# Stage 1 imports
from pipeline_config import PipelineConfig, OrganismType, MoleculeType
from exceptions import PipelineError, ModelInferenceError, SchemaValidationError
from schema_validator import validate_stage_output, STAGE_1_OUTPUT_SCHEMA, STAGE_2_OUTPUT_SCHEMA
from logger_setup import PipelineLogger, setup_logger

# Global logger - will be initialized when run_stage_2 is called
logger: Optional[PipelineLogger] = None


def get_logger(pipeline_id: str = "test", stage: int = 2) -> PipelineLogger:
    """Get or create logger for Stage 2."""
    global logger
    if logger is None:
        logger = setup_logger(pipeline_id, stage)
    return logger


@dataclass
class ReactionRule:
    """Represents a biochemical reaction rule for retrosynthesis."""
    rule_id: str
    name: str
    substrate_pattern: str
    product_pattern: str
    enzyme_class: str
    ec_number: str
    delta_g_kj_per_mol: float
    reversibility: str  # "reversible", "irreversible_forward", "irreversible_backward"
    cofactors: List[str] = field(default_factory=list)
    organisms: List[str] = field(default_factory=list)


@dataclass
class PathwayStep:
    """Represents a single step in a metabolic pathway."""
    step_number: int
    reaction_id: str
    enzyme_name: str
    gene_name: str
    ec_number: str
    substrate: str
    product: str
    delta_g_kj_per_mol: float
    kcat_per_sec: float
    km_mm: float
    is_heterologous: bool
    source_organism: str


@dataclass
class PathwayCandidate:
    """Represents a complete pathway candidate."""
    pathway_id: str
    rank: int
    pathway_name: str
    steps: List[PathwayStep]
    total_steps: int
    predicted_yield_mol_per_mol: float
    thermodynamic_feasibility_score: float
    gnn_viability_score: float
    host_compatibility_score: float


@dataclass
class EnzymeCandidate:
    """Represents a candidate enzyme for a pathway step."""
    enzyme_id: str
    name: str
    gene_name: str
    source_organism: str
    ec_number: str
    sequence_length: int
    kcat_per_sec: float
    km_mm: float
    optimal_ph: float
    optimal_temp_c: float
    stability_score: float
    expression_score: float


@dataclass
class ExpressionCassette:
    """Represents a designed expression cassette."""
    cassette_id: str
    promoter: str
    promoter_strength: float
    rbs: str
    rbs_strength: float
    gene_name: str
    optimized_sequence: str
    terminator: str
    predicted_tpm: float
    cai: float


class RetrosynthesisEngine:
    """
    Engine for predicting metabolic pathways using retrosynthesis analysis.
    Uses Monte Carlo Tree Search (MCTS) to explore pathway space.
    """
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.reaction_rules: List[ReactionRule] = []
        self._load_reaction_rules()
        get_logger().info(f"RetrosynthesisEngine initialized with {len(self.reaction_rules)} reaction rules")
    
    def _load_reaction_rules(self) -> None:
        """Load simulated reaction rules from biochemistry knowledge base."""
        try:
            # Simulated reaction rules covering major metabolic transformations
            rules_data = [
                # Glycolysis/TCA related
                {"rule_id": "RXN001", "name": "Hexokinase", "substrate_pattern": "glucose", "product_pattern": "glucose-6-phosphate", "enzyme_class": "transferase", "ec_number": "2.7.1.1", "delta_g": -16.7, "reversibility": "irreversible_forward", "cofactors": ["ATP", "Mg2+"]},
                {"rule_id": "RXN002", "name": "Phosphoglucose isomerase", "substrate_pattern": "glucose-6-phosphate", "product_pattern": "fructose-6-phosphate", "enzyme_class": "isomerase", "ec_number": "5.3.1.9", "delta_g": 1.7, "reversibility": "reversible", "cofactors": []},
                {"rule_id": "RXN003", "name": "Phosphofructokinase", "substrate_pattern": "fructose-6-phosphate", "product_pattern": "fructose-1,6-bisphosphate", "enzyme_class": "transferase", "ec_number": "2.7.1.11", "delta_g": -14.2, "reversibility": "irreversible_forward", "cofactors": ["ATP", "Mg2+"]},
                {"rule_id": "RXN004", "name": "Aldolase", "substrate_pattern": "fructose-1,6-bisphosphate", "product_pattern": "glyceraldehyde-3-phosphate", "enzyme_class": "lyase", "ec_number": "4.1.2.13", "delta_g": 23.8, "reversibility": "reversible", "cofactors": []},
                {"rule_id": "RXN005", "name": "GAPDH", "substrate_pattern": "glyceraldehyde-3-phosphate", "product_pattern": "1,3-bisphosphoglycerate", "enzyme_class": "oxidoreductase", "ec_number": "1.2.1.12", "delta_g": 6.3, "reversibility": "reversible", "cofactors": ["NAD+", "Pi"]},
                
                # MEP pathway (for terpenes like lycopene)
                {"rule_id": "RXN010", "name": "DXS", "substrate_pattern": "pyruvate", "product_pattern": "DXP", "enzyme_class": "synthase", "ec_number": "2.2.1.7", "delta_g": -8.5, "reversibility": "irreversible_forward", "cofactors": ["TPP"]},
                {"rule_id": "RXN011", "name": "DXR", "substrate_pattern": "DXP", "product_pattern": "MEP", "enzyme_class": "reductoisomerase", "ec_number": "1.1.1.267", "delta_g": -15.2, "reversibility": "irreversible_forward", "cofactors": ["NADPH"]},
                {"rule_id": "RXN012", "name": "IspD", "substrate_pattern": "MEP", "product_pattern": "CDP-ME", "enzyme_class": "cytidylyltransferase", "ec_number": "2.7.7.60", "delta_g": -12.1, "reversibility": "irreversible_forward", "cofactors": ["CTP"]},
                {"rule_id": "RXN013", "name": "IspE", "substrate_pattern": "CDP-ME", "product_pattern": "CDP-MEP", "enzyme_class": "kinase", "ec_number": "2.7.1.148", "delta_g": -18.3, "reversibility": "irreversible_forward", "cofactors": ["ATP"]},
                {"rule_id": "RXN014", "name": "IspF", "substrate_pattern": "CDP-MEP", "product_pattern": "ME-cPP", "enzyme_class": "cyclase", "ec_number": "4.6.1.12", "delta_g": -5.4, "reversibility": "reversible", "cofactors": []},
                {"rule_id": "RXN015", "name": "IspG", "substrate_pattern": "ME-cPP", "product_pattern": "HMBPP", "enzyme_class": "reductase", "ec_number": "1.17.7.1", "delta_g": -22.1, "reversibility": "irreversible_forward", "cofactors": ["Fe-S", "NADPH"]},
                {"rule_id": "RXN016", "name": "IspH", "substrate_pattern": "HMBPP", "product_pattern": "IPP", "enzyme_class": "reductase", "ec_number": "1.17.7.4", "delta_g": -19.8, "reversibility": "irreversible_forward", "cofactors": ["Fe-S", "NADPH"]},
                {"rule_id": "RXN017", "name": "IDI", "substrate_pattern": "IPP", "product_pattern": "DMAPP", "enzyme_class": "isomerase", "ec_number": "5.3.3.2", "delta_g": 3.2, "reversibility": "reversible", "cofactors": []},
                
                # Carotenoid pathway (lycopene)
                {"rule_id": "RXN020", "name": "IspB", "substrate_pattern": "IPP", "product_pattern": "FPP", "enzyme_class": "prenyltransferase", "ec_number": "2.5.1.1", "delta_g": -7.8, "reversibility": "irreversible_forward", "cofactors": ["Mg2+"]},
                {"rule_id": "RXN021", "name": "CrtE", "substrate_pattern": "FPP", "product_pattern": "GGPP", "enzyme_class": "prenyltransferase", "ec_number": "2.5.1.29", "delta_g": -6.2, "reversibility": "irreversible_forward", "cofactors": ["Mg2+"]},
                {"rule_id": "RXN022", "name": "CrtB", "substrate_pattern": "GGPP", "product_pattern": "phytoene", "enzyme_class": "desaturase", "ec_number": "2.5.1.32", "delta_g": -4.5, "reversibility": "irreversible_forward", "cofactors": []},
                {"rule_id": "RXN023", "name": "CrtI", "substrate_pattern": "phytoene", "product_pattern": "lycopene", "enzyme_class": "desaturase", "ec_number": "1.3.5.5", "delta_g": -28.4, "reversibility": "irreversible_forward", "cofactors": ["FAD"]},
                
                # Shikimate pathway (for vanillin, aromatics)
                {"rule_id": "RXN030", "name": "DAHP synthase", "substrate_pattern": "PEP", "product_pattern": "DAHP", "enzyme_class": "synthase", "ec_number": "2.5.1.54", "delta_g": -12.3, "reversibility": "irreversible_forward", "cofactors": ["Mn2+"]},
                {"rule_id": "RXN031", "name": "DHQ synthase", "substrate_pattern": "DAHP", "product_pattern": "DHQ", "enzyme_class": "cyclase", "ec_number": "4.6.1.3", "delta_g": -8.7, "reversibility": "irreversible_forward", "cofactors": ["NAD+", "Co2+"]},
                {"rule_id": "RXN032", "name": "DHS dehydratase", "substrate_pattern": "DHQ", "product_pattern": "shikimate", "enzyme_class": "dehydratase", "ec_number": "4.2.1.10", "delta_g": 2.1, "reversibility": "reversible", "cofactors": []},
                {"rule_id": "RXN033", "name": "Shikimate kinase", "substrate_pattern": "shikimate", "product_pattern": "shikimate-3-phosphate", "enzyme_class": "kinase", "ec_number": "2.7.1.71", "delta_g": -10.5, "reversibility": "irreversible_forward", "cofactors": ["ATP"]},
                {"rule_id": "RXN034", "name": "EPSP synthase", "substrate_pattern": "shikimate-3-phosphate", "product_pattern": "chorismate", "enzyme_class": "transferase", "ec_number": "2.5.1.19", "delta_g": -5.2, "reversibility": "irreversible_forward", "cofactors": []},
                
                # Vanillin specific
                {"rule_id": "RXN040", "name": "Chorismate lyase", "substrate_pattern": "chorismate", "product_pattern": "prephenate", "enzyme_class": "lyase", "ec_number": "4.1.3.27", "delta_g": -3.8, "reversibility": "reversible", "cofactors": []},
                {"rule_id": "RXN041", "name": "Prephenate dehydratase", "substrate_pattern": "prephenate", "product_pattern": "phenylpyruvate", "enzyme_class": "dehydratase", "ec_number": "4.2.1.51", "delta_g": -15.6, "reversibility": "irreversible_forward", "cofactors": []},
                {"rule_id": "RXN042", "name": "Aminotransferase", "substrate_pattern": "phenylpyruvate", "product_pattern": "phenylalanine", "enzyme_class": "transferase", "ec_number": "2.6.1.57", "delta_g": -2.1, "reversibility": "reversible", "cofactors": ["PLP"]},
                {"rule_id": "RXN043", "name": "PAL", "substrate_pattern": "phenylalanine", "product_pattern": "cinnamate", "enzyme_class": "ammonia-lyase", "ec_number": "4.3.1.24", "delta_g": -8.9, "reversibility": "irreversible_forward", "cofactors": []},
                {"rule_id": "RXN044", "name": "C4H", "substrate_pattern": "cinnamate", "product_pattern": "p-coumarate", "enzyme_class": "monooxygenase", "ec_number": "1.14.13.11", "delta_g": -45.2, "reversibility": "irreversible_forward", "cofactors": ["NADPH", "O2"]},
                {"rule_id": "RXN045", "name": "4CL", "substrate_pattern": "p-coumarate", "product_pattern": "coumaroyl-CoA", "enzyme_class": "ligase", "ec_number": "6.2.1.12", "delta_g": -18.7, "reversibility": "irreversible_forward", "cofactors": ["ATP", "CoA"]},
                
                # Lysine biosynthesis (DAP pathway)
                {"rule_id": "RXN050", "name": "Citrate synthase", "substrate_pattern": "oxaloacetate", "product_pattern": "citrate", "enzyme_class": "synthase", "ec_number": "2.3.3.1", "delta_g": -31.4, "reversibility": "irreversible_forward", "cofactors": []},
                {"rule_id": "RXN051", "name": "Aconitase", "substrate_pattern": "citrate", "product_pattern": "isocitrate", "enzyme_class": "isomerase", "ec_number": "4.2.1.3", "delta_g": 5.2, "reversibility": "reversible", "cofactors": ["Fe-S"]},
                {"rule_id": "RXN052", "name": "Isocitrate dehydrogenase", "substrate_pattern": "isocitrate", "product_pattern": "alpha-ketoglutarate", "enzyme_class": "dehydrogenase", "ec_number": "1.1.1.42", "delta_g": -21.3, "reversibility": "irreversible_forward", "cofactors": ["NADP+"]},
                {"rule_id": "RXN053", "name": "Alpha-KG dehydrogenase", "substrate_pattern": "alpha-ketoglutarate", "product_pattern": "succinyl-CoA", "enzyme_class": "dehydrogenase", "ec_number": "1.2.4.2", "delta_g": -33.5, "reversibility": "irreversible_forward", "cofactors": ["NAD+", "CoA", "TPP"]},
                {"rule_id": "RXN054", "name": "Succinyl-CoA synthetase", "substrate_pattern": "succinyl-CoA", "product_pattern": "succinate", "enzyme_class": "ligase", "ec_number": "6.2.1.5", "delta_g": -2.9, "reversibility": "reversible", "cofactors": ["GDP", "Pi"]},
                {"rule_id": "RXN055", "name": "Fumarase", "substrate_pattern": "succinate", "product_pattern": "fumarate", "enzyme_class": "hydratase", "ec_number": "4.2.1.2", "delta_g": -3.4, "reversibility": "reversible", "cofactors": []},
                {"rule_id": "RXN056", "name": "Malate dehydrogenase", "substrate_pattern": "fumarate", "product_pattern": "malate", "enzyme_class": "dehydrogenase", "ec_number": "1.3.5.4", "delta_g": -12.1, "reversibility": "reversible", "cofactors": ["FAD"]},
                {"rule_id": "RXN057", "name": "Malate dehydrogenase", "substrate_pattern": "malate", "product_pattern": "oxaloacetate", "enzyme_class": "dehydrogenase", "ec_number": "1.1.1.37", "delta_g": 29.7, "reversibility": "reversible", "cofactors": ["NAD+"]},
                
                # PHA biosynthesis
                {"rule_id": "RXN060", "name": "Beta-ketothiolase", "substrate_pattern": "acetyl-CoA", "product_pattern": "acetoacetyl-CoA", "enzyme_class": "thiolase", "ec_number": "2.3.1.9", "delta_g": -28.5, "reversibility": "reversible", "cofactors": []},
                {"rule_id": "RXN061", "name": "Acetoacetyl-CoA reductase", "substrate_pattern": "acetoacetyl-CoA", "product_pattern": "3HB-CoA", "enzyme_class": "reductase", "ec_number": "1.1.1.36", "delta_g": -25.3, "reversibility": "reversible", "cofactors": ["NADPH"]},
                {"rule_id": "RXN062", "name": "PHA synthase", "substrate_pattern": "3HB-CoA", "product_pattern": "PHA", "enzyme_class": "polymerase", "ec_number": "2.3.1.308", "delta_g": -15.8, "reversibility": "irreversible_forward", "cofactors": []},
            ]
            
            for rule in rules_data:
                reaction_rule = ReactionRule(
                    rule_id=rule["rule_id"],
                    name=rule["name"],
                    substrate_pattern=rule["substrate_pattern"],
                    product_pattern=rule["product_pattern"],
                    enzyme_class=rule["enzyme_class"],
                    ec_number=rule["ec_number"],
                    delta_g_kj_per_mol=rule["delta_g"],
                    reversibility=rule["reversibility"],
                    cofactors=rule.get("cofactors", []),
                    organisms=["E. coli", "S. cerevisiae", "B. subtilis", "C. glutamicum", "P. putida"]
                )
                self.reaction_rules.append(reaction_rule)
            
            get_logger().debug(f"Loaded {len(self.reaction_rules)} reaction rules")
            
        except Exception as e:
            get_logger().error(f"Failed to load reaction rules: {type(e).__name__}: {e}")
            raise ModelInferenceError(f"Reaction rule loading failed: {e}")
    
    def run_mcts(self, target_molecule: str, starting_metabolites: List[str], 
                 max_iterations: int = 100) -> List[List[ReactionRule]]:
        """
        Run Monte Carlo Tree Search to find pathways from starting metabolites to target.
        
        Args:
            target_molecule: Target product name
            starting_metabolites: List of available starting metabolites
            max_iterations: Maximum MCTS iterations
            
        Returns:
            List of pathway candidates (each pathway is a list of ReactionRules)
        """
        try:
            get_logger().info(f"Running MCTS for target={target_molecule} with {max_iterations} iterations")
            pathways = []
            
            # Simplified MCTS simulation - in reality this would be much more complex
            # Find rules that produce the target
            target_rules = [r for r in self.reaction_rules 
                          if target_molecule.lower() in r.product_pattern.lower()]
            
            for rule in target_rules[:3]:  # Take top 3 producing rules
                pathway = [rule]
                current_substrates = [rule.substrate_pattern]
                
                # Backtrack to find precursors
                depth = 0
                while current_substrates and depth < 10:
                    new_substrates = []
                    for sub in current_substrates:
                        # Check if substrate is in starting metabolites
                        is_starting = any(start.lower() in sub.lower() for start in starting_metabolites)
                        if not is_starting:
                            # Find rules that produce this substrate
                            precursor_rules = [r for r in self.reaction_rules 
                                             if sub.lower() in r.product_pattern.lower() 
                                             and r not in pathway]
                            if precursor_rules:
                                selected_rule = random.choice(precursor_rules[:2])
                                pathway.insert(0, selected_rule)
                                new_substrates.append(selected_rule.substrate_pattern)
                    
                    current_substrates = new_substrates if new_substrates else []
                    depth += 1
                
                if pathway:
                    pathways.append(pathway)
            
            get_logger().debug(f"MCTS found {len(pathways)} pathway candidates")
            return pathways
            
        except Exception as e:
            get_logger().error(f"MCTS failed: {type(e).__name__}: {e}")
            # Fallback: return simple direct pathways
            return self._fallback_pathway_generation(target_molecule, starting_metabolites)
    
    def _fallback_pathway_generation(self, target_molecule: str, 
                                     starting_metabolites: List[str]) -> List[List[ReactionRule]]:
        """Fallback method when MCTS fails - use rule-based pathway generation."""
        get_logger().warning("MCTS failed, using fallback rule-based pathway generation")
        
        # Predefined pathways for common targets
        predefined = {
            "lycopene": ["RXN010", "RXN011", "RXN012", "RXN013", "RXN014", "RXN015", "RXN016", "RXN017", "RXN020", "RXN021", "RXN022", "RXN023"],
            "vanillin": ["RXN030", "RXN031", "RXN032", "RXN033", "RXN034", "RXN040", "RXN041", "RXN043", "RXN044"],
            "l_lysine": ["RXN050", "RXN051", "RXN052", "RXN053", "RXN054", "RXN055", "RXN056", "RXN057"],
            "lysine": ["RXN050", "RXN051", "RXN052", "RXN053", "RXN054", "RXN055", "RXN056", "RXN057"],
            "pha": ["RXN060", "RXN061", "RXN062"],
            "riboflavin": ["RXN001", "RXN002", "RXN003", "RXN004", "RXN005"],
            "artemisinic_acid": ["RXN030", "RXN031", "RXN032", "RXN010", "RXN011", "RXN070", "RXN071"],
            "l_glutamate": ["RXN050", "RXN051", "RXN080", "RXN081"],
            "l_threonine": ["RXN090", "RXN091", "RXN092", "RXN093"],
            "hyaluronic_acid": ["RXN001", "RXN100", "RXN101", "RXN102"],
        }
        
        pathways = []
        target_key = target_molecule.lower()
        
        if target_key in predefined:
            pathway_rules = []
            for rule_id in predefined[target_key]:
                matching_rules = [r for r in self.reaction_rules if r.rule_id == rule_id]
                if matching_rules:
                    pathway_rules.append(matching_rules[0])
            if pathway_rules:
                pathways.append(pathway_rules)
        
        return pathways
    
    def score_thermodynamics(self, pathway: List[ReactionRule]) -> float:
        """
        Calculate thermodynamic feasibility score for a pathway.
        
        Args:
            pathway: List of ReactionRules forming a pathway
            
        Returns:
            Thermodynamic feasibility score (0-1, higher is better)
        """
        try:
            total_delta_g = sum(rule.delta_g_kj_per_mol for rule in pathway)
            
            # Score based on overall exergonicity and individual step feasibility
            irreversible_steps = sum(1 for r in pathway if "irreversible_forward" in r.reversibility)
            
            # Favor pathways with negative overall delta_G
            energy_score = max(0, min(1, 1 - (total_delta_g / 200)))
            
            # Favor pathways with some irreversible steps (driving force)
            irreversibility_score = min(1, irreversible_steps / len(pathway)) if pathway else 0
            
            final_score = 0.7 * energy_score + 0.3 * irreversibility_score
            final_score = max(0, min(1, final_score))  # Ensure score is in [0, 1]
            
            get_logger().debug(f"Thermodynamic score: {final_score:.3f} (ΔG={total_delta_g:.1f} kJ/mol)")
            return final_score
            
        except Exception as e:
            get_logger().error(f"Thermodynamic scoring failed: {type(e).__name__}: {e}")
            return 0.5  # Neutral fallback score
    
    def rank_pathways(self, pathways: List[List[ReactionRule]], 
                      organism: str) -> List[PathwayCandidate]:
        """
        Rank pathways based on multiple criteria.
        
        Args:
            pathways: List of pathway candidates (each is a list of ReactionRules)
            organism: Target organism name
            
        Returns:
            Ranked list of PathwayCandidates
        """
        try:
            ranked = []
            
            for idx, pathway in enumerate(pathways):
                thermo_score = self.score_thermodynamics(pathway)
                
                # Simulate GNN viability score (would use graph neural network in reality)
                gnn_score = 0.6 + 0.3 * random.random()
                
                # Host compatibility based on organism presence in rules
                host_matches = sum(1 for r in pathway if organism in r.organisms)
                host_score = host_matches / len(pathway) if pathway else 0
                
                # Predicted yield based on pathway length and thermodynamics
                length_penalty = 0.95 ** len(pathway)  # Longer pathways have lower yield
                predicted_yield = 0.8 * thermo_score * length_penalty
                
                pathway_candidate = PathwayCandidate(
                    pathway_id=f"PATH_{uuid.uuid4().hex[:8]}",
                    rank=idx + 1,
                    pathway_name=f"{pathway[0].name} → ... → {pathway[-1].product_pattern}",
                    steps=[],  # Will be filled by enzyme selector
                    total_steps=len(pathway),
                    predicted_yield_mol_per_mol=predicted_yield,
                    thermodynamic_feasibility_score=thermo_score,
                    gnn_viability_score=gnn_score,
                    host_compatibility_score=host_score
                )
                ranked.append((pathway_candidate, pathway))
            
            # Sort by combined score
            ranked.sort(key=lambda x: (x[0].thermodynamic_feasibility_score * 0.4 + 
                                       x[0].gnn_viability_score * 0.3 + 
                                       x[0].host_compatibility_score * 0.3), 
                         reverse=True)
            
            # Update ranks
            for idx, (candidate, _) in enumerate(ranked):
                candidate.rank = idx + 1
            
            get_logger().info(f"Ranked {len(ranked)} pathways for {organism}")
            return [c for c, _ in ranked]
            
        except Exception as e:
            get_logger().error(f"Pathway ranking failed: {type(e).__name__}: {e}")
            raise ModelInferenceError(f"Pathway ranking failed: {e}")


class EnzymeSelector:
    """
    Selects optimal enzymes for each pathway step based on kinetics,
    host compatibility, and sequence features.
    """
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.enzyme_database: Dict[str, List[EnzymeCandidate]] = {}
        self._load_enzyme_database()
        get_logger().info("EnzymeSelector initialized")
    
    def _load_enzyme_database(self) -> None:
        """Load simulated enzyme database with kinetic parameters."""
        try:
            # Simulated enzyme candidates for key reactions
            enzymes_data = {
                "DXS": [
                    {"name": "dxs_Ec", "gene": "dxs", "source": "E. coli", "ec": "2.2.1.7", "length": 712, "kcat": 2.5, "km": 0.15, "ph": 7.5, "temp": 37, "stability": 0.85, "expression": 0.9},
                    {"name": "dxs_Bs", "gene": "dxs", "source": "B. subtilis", "ec": "2.2.1.7", "length": 698, "kcat": 1.8, "km": 0.22, "ph": 7.0, "temp": 37, "stability": 0.78, "expression": 0.75},
                    {"name": "dxs_Pp", "gene": "dxs", "source": "P. putida", "ec": "2.2.1.7", "length": 705, "kcat": 2.1, "km": 0.18, "ph": 7.2, "temp": 30, "stability": 0.82, "expression": 0.85},
                ],
                "CrtE": [
                    {"name": "crtE_Pa", "gene": "crtE", "source": "Pantoea ananatis", "ec": "2.5.1.29", "length": 295, "kcat": 5.2, "km": 0.08, "ph": 7.5, "temp": 30, "stability": 0.88, "expression": 0.82},
                    {"name": "crtE_Ec", "gene": "idi", "source": "E. coli", "ec": "2.5.1.29", "length": 288, "kcat": 3.8, "km": 0.12, "ph": 7.5, "temp": 37, "stability": 0.85, "expression": 0.88},
                ],
                "CrtB": [
                    {"name": "crtB_Er", "gene": "crtB", "source": "Erwinia uredovora", "ec": "2.5.1.32", "length": 442, "kcat": 1.5, "km": 0.05, "ph": 7.5, "temp": 30, "stability": 0.82, "expression": 0.78},
                    {"name": "crtB_Pa", "gene": "crtB", "source": "Pantoea ananatis", "ec": "2.5.1.32", "length": 438, "kcat": 1.8, "km": 0.04, "ph": 7.5, "temp": 30, "stability": 0.85, "expression": 0.82},
                ],
                "CrtI": [
                    {"name": "crtI_Pa", "gene": "crtI", "source": "Pantoea ananatis", "ec": "1.3.5.5", "length": 535, "kcat": 0.8, "km": 0.02, "ph": 7.5, "temp": 30, "stability": 0.80, "expression": 0.75},
                    {"name": "crtI_Er", "gene": "crtI", "source": "Erwinia uredovora", "ec": "1.3.5.5", "length": 532, "kcat": 0.6, "km": 0.03, "ph": 7.5, "temp": 30, "stability": 0.78, "expression": 0.72},
                ],
                "PAL": [
                    {"name": "pal_At", "gene": "PAL1", "source": "Arabidopsis thaliana", "ec": "4.3.1.24", "length": 712, "kcat": 0.5, "km": 0.8, "ph": 8.5, "temp": 30, "stability": 0.72, "expression": 0.65},
                    {"name": "pal_Rg", "gene": "pal", "source": "Rhodotorula glutinis", "ec": "4.3.1.24", "length": 658, "kcat": 0.8, "km": 0.5, "ph": 8.0, "temp": 30, "stability": 0.78, "expression": 0.72},
                ],
                "C4H": [
                    {"name": "c4h_At", "gene": "CYP73A1", "source": "Arabidopsis thaliana", "ec": "1.14.13.11", "length": 502, "kcat": 0.3, "km": 0.05, "ph": 7.5, "temp": 30, "stability": 0.68, "expression": 0.60},
                ],
                "4CL": [
                    {"name": "4cl_At", "gene": "4CL1", "source": "Arabidopsis thaliana", "ec": "6.2.1.12", "length": 532, "kcat": 2.5, "km": 0.12, "ph": 7.5, "temp": 30, "stability": 0.82, "expression": 0.78},
                ],
            }
            
            for enzyme_name, candidates in enzymes_data.items():
                self.enzyme_database[enzyme_name] = [
                    EnzymeCandidate(
                        enzyme_id=f"ENZ_{uuid.uuid4().hex[:8]}",
                        name=c["name"],
                        gene_name=c["gene"],
                        source_organism=c["source"],
                        ec_number=c["ec"],
                        sequence_length=c["length"],
                        kcat_per_sec=c["kcat"],
                        km_mm=c["km"],
                        optimal_ph=c["ph"],
                        optimal_temp_c=c["temp"],
                        stability_score=c["stability"],
                        expression_score=c["expression"]
                    )
                    for c in candidates
                ]
            
            get_logger().debug(f"Loaded enzyme database with {sum(len(v) for v in self.enzyme_database.values())} candidates")
            
        except Exception as e:
            get_logger().error(f"Failed to load enzyme database: {type(e).__name__}: {e}")
            raise ModelInferenceError(f"Enzyme database loading failed: {e}")
    
    def query_brenda_simulated(self, enzyme_name: str, ec_number: str) -> List[EnzymeCandidate]:
        """
        Simulate querying BRENDA database for enzyme kinetic data.
        
        Args:
            enzyme_name: Name of the enzyme
            ec_number: EC number of the enzyme
            
        Returns:
            List of EnzymeCandidates with kinetic parameters
        """
        try:
            # First check our local database
            if enzyme_name in self.enzyme_database:
                candidates = self.enzyme_database[enzyme_name]
                get_logger().debug(f"Found {len(candidates)} candidates for {enzyme_name}")
                return candidates
            
            # Generate simulated candidates if not in database
            num_candidates = random.randint(2, 5)
            candidates = []
            
            for i in range(num_candidates):
                source_org = random.choice(["E. coli", "S. cerevisiae", "B. subtilis", 
                                           "C. glutamicum", "P. putida", 
                                           "Pantoea ananatis", "Arabidopsis thaliana"])
                
                candidate = EnzymeCandidate(
                    enzyme_id=f"ENZ_{uuid.uuid4().hex[:8]}",
                    name=f"{enzyme_name.lower()}_{source_org[:2]}",
                    gene_name=enzyme_name.lower(),
                    source_organism=source_org,
                    ec_number=ec_number,
                    sequence_length=random.randint(300, 800),
                    kcat_per_sec=round(random.uniform(0.1, 10.0), 2),
                    km_mm=round(random.uniform(0.01, 1.0), 3),
                    optimal_ph=round(random.uniform(6.5, 8.5), 1),
                    optimal_temp_c=round(random.uniform(25, 40), 1),
                    stability_score=round(random.uniform(0.6, 0.95), 2),
                    expression_score=round(random.uniform(0.5, 0.95), 2)
                )
                candidates.append(candidate)
            
            get_logger().debug(f"Generated {len(candidates)} simulated candidates for {enzyme_name}")
            return candidates
            
        except Exception as e:
            get_logger().error(f"BRENDA query failed for {enzyme_name}: {type(e).__name__}: {e}")
            # Return minimal fallback candidate
            return [EnzymeCandidate(
                enzyme_id=f"ENZ_FALLBACK_{uuid.uuid4().hex[:8]}",
                name=f"{enzyme_name}_fallback",
                gene_name=enzyme_name.lower(),
                source_organism="E. coli",
                ec_number=ec_number,
                sequence_length=500,
                kcat_per_sec=1.0,
                km_mm=0.1,
                optimal_ph=7.5,
                optimal_temp_c=37.0,
                stability_score=0.7,
                expression_score=0.7
            )]
    
    def predict_host_compatibility(self, enzyme: EnzymeCandidate, 
                                   host_organism: str) -> float:
        """
        Predict enzyme compatibility with host organism.
        
        Args:
            enzyme: EnzymeCandidate to evaluate
            host_organism: Target host organism
            
        Returns:
            Compatibility score (0-1)
        """
        try:
            score = 0.5  # Base score
            
            # Same organism = highest compatibility
            if enzyme.source_organism == host_organism:
                score = 0.95
            # Related organisms
            elif (host_organism == "E. coli" and enzyme.source_organism == "P. putida") or \
                 (host_organism == "B. subtilis" and enzyme.source_organism == "C. glutamicum"):
                score = 0.85
            # Eukaryotic to prokaryotic (may need codon optimization)
            elif host_organism in ["E. coli", "B. subtilis", "C. glutamicum", "P. putida"] and \
                 enzyme.source_organism in ["S. cerevisiae", "Arabidopsis thaliana"]:
                score = 0.65
            # Prokaryotic to eukaryotic
            elif host_organism == "S. cerevisiae" and \
                 enzyme.source_organism in ["E. coli", "B. subtilis"]:
                score = 0.70
            
            # Adjust for temperature compatibility
            temp_diff = abs(enzyme.optimal_temp_c - 37)  # Assume 37°C fermentation
            temp_penalty = min(0.2, temp_diff * 0.01)
            score -= temp_penalty
            
            # Adjust for pH compatibility
            ph_diff = abs(enzyme.optimal_ph - 7.2)  # Assume pH 7.2 fermentation
            ph_penalty = min(0.1, ph_diff * 0.05)
            score -= ph_penalty
            
            return max(0, min(1, score))
            
        except Exception as e:
            get_logger().error(f"Host compatibility prediction failed: {type(e).__name__}: {e}")
            return 0.5  # Neutral fallback
    
    def score_esm2_embedding(self, enzyme: EnzymeCandidate) -> float:
        """
        Simulate ESM-2 protein language model embedding score.
        In reality, this would use actual transformer embeddings.
        
        Args:
            enzyme: EnzymeCandidate to score
            
        Returns:
            ESM-2 similarity score (0-1)
        """
        try:
            # Simulate embedding-based scoring
            # Higher for well-characterized enzymes from model organisms
            base_score = 0.6
            
            if enzyme.source_organism in ["E. coli", "S. cerevisiae"]:
                base_score = 0.8
            elif enzyme.source_organism in ["B. subtilis", "C. glutamicum", "P. putida"]:
                base_score = 0.75
            
            # Add some randomness to simulate embedding variation
            embedding_score = base_score + 0.15 * random.random()
            
            # Stability contributes to embedding quality
            embedding_score = 0.7 * embedding_score + 0.3 * enzyme.stability_score
            
            get_logger().debug(f"ESM-2 score for {enzyme.name}: {embedding_score:.3f}")
            return min(1, embedding_score)
            
        except Exception as e:
            get_logger().error(f"ESM-2 scoring failed: {type(e).__name__}: {e}")
            return 0.6  # Neutral fallback
    
    def filter_by_organism_constraints(self, candidates: List[EnzymeCandidate],
                                       organism: str) -> List[EnzymeCandidate]:
        """
        Filter enzyme candidates based on organism-specific constraints.
        
        Args:
            candidates: List of EnzymeCandidates
            organism: Target organism
            
        Returns:
            Filtered list of compatible candidates
        """
        try:
            filtered = []
            
            for candidate in candidates:
                # Temperature constraints
                if organism == "E. coli" and candidate.optimal_temp_c > 42:
                    continue
                if organism == "P. putida" and candidate.optimal_temp_c < 25:
                    continue
                
                # pH constraints
                if organism == "C. glutamicum" and candidate.optimal_ph < 6.5:
                    continue
                
                # Gram stain considerations (simplified)
                gram_positive = organism in ["B. subtilis", "C. glutamicum"]
                if gram_positive and candidate.expression_score < 0.5:
                    continue  # Poor expression in Gram+
                
                filtered.append(candidate)
            
            get_logger().debug(f"Filtered {len(candidates)} → {len(filtered)} candidates for {organism}")
            return filtered if filtered else candidates  # Return original if all filtered out
            
        except Exception as e:
            get_logger().error(f"Organism constraint filtering failed: {type(e).__name__}: {e}")
            return candidates


class CodonOptimizer:
    """
    Optimizes coding sequences for expression in different host organisms.
    """
    
    # Codon usage tables (simplified, % frequency)
    CODON_TABLES = {
        "E. coli": {
            'A': {'GCT': 0.42, 'GCC': 0.32, 'GCA': 0.18, 'GCG': 0.08},
            'R': {'CGT': 0.40, 'CGC': 0.28, 'CGA': 0.12, 'CGG': 0.10, 'AGA': 0.06, 'AGG': 0.04},
            'N': {'AAT': 0.45, 'AAC': 0.55},
            'D': {'GAT': 0.48, 'GAC': 0.52},
            'C': {'TGT': 0.47, 'TGC': 0.53},
            'Q': {'CAA': 0.58, 'CAG': 0.42},
            'E': {'GAA': 0.68, 'GAG': 0.32},
            'G': {'GGT': 0.40, 'GGC': 0.35, 'GGA': 0.15, 'GGG': 0.10},
            'H': {'CAT': 0.43, 'CAC': 0.57},
            'I': {'ATT': 0.42, 'ATC': 0.38, 'ATA': 0.20},
            'L': {'TTA': 0.10, 'TTG': 0.13, 'CTT': 0.13, 'CTC': 0.10, 'CTA': 0.04, 'CTG': 0.50},
            'K': {'AAA': 0.75, 'AAG': 0.25},
            'M': {'ATG': 1.0},
            'F': {'TTT': 0.55, 'TTC': 0.45},
            'P': {'CCT': 0.35, 'CCC': 0.20, 'CCA': 0.25, 'CCG': 0.20},
            'S': {'TCT': 0.20, 'TCC': 0.18, 'TCA': 0.15, 'TCG': 0.12, 'AGT': 0.18, 'AGC': 0.17},
            'T': {'ACT': 0.25, 'ACC': 0.38, 'ACA': 0.22, 'ACG': 0.15},
            'W': {'TGG': 1.0},
            'Y': {'TAT': 0.58, 'TAC': 0.42},
            'V': {'GTT': 0.35, 'GTC': 0.28, 'GTA': 0.18, 'GTG': 0.19},
            '*': {'TAA': 0.65, 'TAG': 0.05, 'TGA': 0.30},
        },
        "S. cerevisiae": {
            'A': {'GCT': 0.25, 'GCC': 0.15, 'GCA': 0.35, 'GCG': 0.25},
            'R': {'CGT': 0.15, 'CGC': 0.10, 'CGA': 0.15, 'CGG': 0.10, 'AGA': 0.25, 'AGG': 0.25},
            'N': {'AAT': 0.52, 'AAC': 0.48},
            'D': {'GAT': 0.55, 'GAC': 0.45},
            'C': {'TGT': 0.52, 'TGC': 0.48},
            'Q': {'CAA': 0.72, 'CAG': 0.28},
            'E': {'GAA': 0.70, 'GAG': 0.30},
            'G': {'GGT': 0.35, 'GGC': 0.20, 'GGA': 0.25, 'GGG': 0.20},
            'H': {'CAT': 0.42, 'CAC': 0.58},
            'I': {'ATT': 0.40, 'ATC': 0.35, 'ATA': 0.25},
            'L': {'TTA': 0.20, 'TTG': 0.20, 'CTT': 0.15, 'CTC': 0.10, 'CTA': 0.05, 'CTG': 0.30},
            'K': {'AAA': 0.80, 'AAG': 0.20},
            'M': {'ATG': 1.0},
            'F': {'TTT': 0.60, 'TTC': 0.40},
            'P': {'CCT': 0.30, 'CCC': 0.18, 'CCA': 0.32, 'CCG': 0.20},
            'S': {'TCT': 0.22, 'TCC': 0.15, 'TCA': 0.25, 'TCG': 0.13, 'AGT': 0.15, 'AGC': 0.10},
            'T': {'ACT': 0.28, 'ACC': 0.22, 'ACA': 0.30, 'ACG': 0.20},
            'W': {'TGG': 1.0},
            'Y': {'TAT': 0.62, 'TAC': 0.38},
            'V': {'GTT': 0.25, 'GTC': 0.18, 'GTA': 0.32, 'GTG': 0.25},
            '*': {'TAA': 0.70, 'TAG': 0.08, 'TGA': 0.22},
        },
        "B. subtilis": {
            'A': {'GCT': 0.30, 'GCC': 0.25, 'GCA': 0.25, 'GCG': 0.20},
            'R': {'CGT': 0.20, 'CGC': 0.18, 'CGA': 0.15, 'CGG': 0.12, 'AGA': 0.18, 'AGG': 0.17},
            'N': {'AAT': 0.48, 'AAC': 0.52},
            'D': {'GAT': 0.52, 'GAC': 0.48},
            'C': {'TGT': 0.45, 'TGC': 0.55},
            'Q': {'CAA': 0.65, 'CAG': 0.35},
            'E': {'GAA': 0.62, 'GAG': 0.38},
            'G': {'GGT': 0.32, 'GGC': 0.28, 'GGA': 0.22, 'GGG': 0.18},
            'H': {'CAT': 0.45, 'CAC': 0.55},
            'I': {'ATT': 0.38, 'ATC': 0.35, 'ATA': 0.27},
            'L': {'TTA': 0.15, 'TTG': 0.18, 'CTT': 0.15, 'CTC': 0.12, 'CTA': 0.08, 'CTG': 0.32},
            'K': {'AAA': 0.72, 'AAG': 0.28},
            'M': {'ATG': 1.0},
            'F': {'TTT': 0.52, 'TTC': 0.48},
            'P': {'CCT': 0.28, 'CCC': 0.22, 'CCA': 0.28, 'CCG': 0.22},
            'S': {'TCT': 0.25, 'TCC': 0.20, 'TCA': 0.20, 'TCG': 0.15, 'AGT': 0.12, 'AGC': 0.08},
            'T': {'ACT': 0.30, 'ACC': 0.28, 'ACA': 0.25, 'ACG': 0.17},
            'W': {'TGG': 1.0},
            'Y': {'TAT': 0.55, 'TAC': 0.45},
            'V': {'GTT': 0.30, 'GTC': 0.25, 'GTA': 0.25, 'GTG': 0.20},
            '*': {'TAA': 0.68, 'TAG': 0.07, 'TGA': 0.25},
        },
        "C. glutamicum": {
            'A': {'GCT': 0.28, 'GCC': 0.30, 'GCA': 0.22, 'GCG': 0.20},
            'R': {'CGT': 0.18, 'CGC': 0.22, 'CGA': 0.12, 'CGG': 0.15, 'AGA': 0.18, 'AGG': 0.15},
            'N': {'AAT': 0.45, 'AAC': 0.55},
            'D': {'GAT': 0.50, 'GAC': 0.50},
            'C': {'TGT': 0.42, 'TGC': 0.58},
            'Q': {'CAA': 0.60, 'CAG': 0.40},
            'E': {'GAA': 0.58, 'GAG': 0.42},
            'G': {'GGT': 0.30, 'GGC': 0.32, 'GGA': 0.20, 'GGG': 0.18},
            'H': {'CAT': 0.42, 'CAC': 0.58},
            'I': {'ATT': 0.35, 'ATC': 0.40, 'ATA': 0.25},
            'L': {'TTA': 0.12, 'TTG': 0.15, 'CTT': 0.12, 'CTC': 0.15, 'CTA': 0.08, 'CTG': 0.38},
            'K': {'AAA': 0.68, 'AAG': 0.32},
            'M': {'ATG': 1.0},
            'F': {'TTT': 0.48, 'TTC': 0.52},
            'P': {'CCT': 0.25, 'CCC': 0.25, 'CCA': 0.25, 'CCG': 0.25},
            'S': {'TCT': 0.22, 'TCC': 0.22, 'TCA': 0.18, 'TCG': 0.15, 'AGT': 0.12, 'AGC': 0.11},
            'T': {'ACT': 0.28, 'ACC': 0.32, 'ACA': 0.22, 'ACG': 0.18},
            'W': {'TGG': 1.0},
            'Y': {'TAT': 0.50, 'TAC': 0.50},
            'V': {'GTT': 0.28, 'GTC': 0.30, 'GTA': 0.22, 'GTG': 0.20},
            '*': {'TAA': 0.62, 'TAG': 0.08, 'TGA': 0.30},
        },
        "P. putida": {
            'A': {'GCT': 0.35, 'GCC': 0.30, 'GCA': 0.20, 'GCG': 0.15},
            'R': {'CGT': 0.25, 'CGC': 0.22, 'CGA': 0.15, 'CGG': 0.12, 'AGA': 0.15, 'AGG': 0.11},
            'N': {'AAT': 0.50, 'AAC': 0.50},
            'D': {'GAT': 0.52, 'GAC': 0.48},
            'C': {'TGT': 0.50, 'TGC': 0.50},
            'Q': {'CAA': 0.62, 'CAG': 0.38},
            'E': {'GAA': 0.65, 'GAG': 0.35},
            'G': {'GGT': 0.35, 'GGC': 0.30, 'GGA': 0.20, 'GGG': 0.15},
            'H': {'CAT': 0.45, 'CAC': 0.55},
            'I': {'ATT': 0.40, 'ATC': 0.38, 'ATA': 0.22},
            'L': {'TTA': 0.12, 'TTG': 0.15, 'CTT': 0.15, 'CTC': 0.13, 'CTA': 0.05, 'CTG': 0.40},
            'K': {'AAA': 0.70, 'AAG': 0.30},
            'M': {'ATG': 1.0},
            'F': {'TTT': 0.55, 'TTC': 0.45},
            'P': {'CCT': 0.30, 'CCC': 0.22, 'CCA': 0.28, 'CCG': 0.20},
            'S': {'TCT': 0.25, 'TCC': 0.20, 'TCA': 0.18, 'TCG': 0.12, 'AGT': 0.15, 'AGC': 0.10},
            'T': {'ACT': 0.30, 'ACC': 0.32, 'ACA': 0.22, 'ACG': 0.16},
            'W': {'TGG': 1.0},
            'Y': {'TAT': 0.58, 'TAC': 0.42},
            'V': {'GTT': 0.32, 'GTC': 0.28, 'GTA': 0.22, 'GTG': 0.18},
            '*': {'TAA': 0.65, 'TAG': 0.06, 'TGA': 0.29},
        },
    }
    
    # Restriction sites to avoid
    RESTRICTION_SITES = {
        "BamHI": "GGATCC",
        "EcoRI": "GAATTC",
        "HindIII": "AAGCTT",
        "XhoI": "CTCGAG",
        "NdeI": "CATATG",
        "NcoI": "CCATGG",
    }
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        get_logger().info("CodonOptimizer initialized")
    
    def optimize_sequence(self, protein_sequence: str, organism: str) -> str:
        """
        Optimize a protein sequence for expression in a specific organism.
        
        Args:
            protein_sequence: Amino acid sequence (single-letter codes)
            organism: Target organism name
            
        Returns:
            Optimized DNA sequence
        """
        try:
            codon_table = self.CODON_TABLES.get(organism, self.CODON_TABLES["E. coli"])
            dna_sequence = []
            
            for aa in protein_sequence.upper():
                if aa in codon_table:
                    # Select most frequent codon
                    codons = codon_table[aa]
                    best_codon = max(codons.keys(), key=lambda c: codons[c])
                    dna_sequence.append(best_codon)
                elif aa == '*':  # Stop codon
                    codons = codon_table.get('*', {'TAA': 0.65, 'TAG': 0.05, 'TGA': 0.30})
                    best_codon = max(codons.keys(), key=lambda c: codons[c])
                    dna_sequence.append(best_codon)
                else:
                    # Unknown amino acid, use NNN
                    dna_sequence.append('NNN')
                    get_logger().warning(f"Unknown amino acid '{aa}', using NNN")
            
            optimized = ''.join(dna_sequence)
            get_logger().debug(f"Optimized sequence length: {len(optimized)} bp for {organism}")
            return optimized
            
        except Exception as e:
            get_logger().error(f"Sequence optimization failed: {type(e).__name__}: {e}")
            raise ModelInferenceError(f"Codon optimization failed: {e}")
    
    def calculate_cai(self, dna_sequence: str, organism: str) -> float:
        """
        Calculate Codon Adaptation Index (CAI) for a DNA sequence.
        
        Args:
            dna_sequence: DNA sequence
            organism: Reference organism for codon usage
            
        Returns:
            CAI value (0-1, higher is better)
        """
        try:
            codon_table = self.CODON_TABLES.get(organism, self.CODON_TABLES["E. coli"])
            
            # Get reference set (most frequent codons)
            reference = {}
            for aa, codons in codon_table.items():
                if codons:
                    max_freq = max(codons.values())
                    for codon, freq in codons.items():
                        reference[codon] = freq / max_freq if max_freq > 0 else 0
            
            # Calculate CAI
            codons = [dna_sequence[i:i+3] for i in range(0, len(dna_sequence)-2, 3)]
            valid_codons = [c for c in codons if c in reference and 'N' not in c]
            
            if not valid_codons:
                return 0.5  # Neutral fallback
            
            product = 1.0
            for codon in valid_codons:
                product *= reference.get(codon, 0.5)
            
            cai = product ** (1 / len(valid_codons))
            get_logger().debug(f"CAI for {organism}: {cai:.3f}")
            return cai
            
        except Exception as e:
            get_logger().error(f"CAI calculation failed: {type(e).__name__}: {e}")
            return 0.5  # Neutral fallback
    
    def check_restriction_sites(self, dna_sequence: str) -> Dict[str, List[int]]:
        """
        Check for restriction enzyme sites in a DNA sequence.
        
        Args:
            dna_sequence: DNA sequence to check
            
        Returns:
            Dict mapping restriction enzyme names to list of positions
        """
        try:
            found_sites = {}
            dna_upper = dna_sequence.upper()
            
            for enzyme, site in self.RESTRICTION_SITES.items():
                positions = []
                pos = 0
                while True:
                    pos = dna_upper.find(site, pos)
                    if pos == -1:
                        break
                    positions.append(pos)
                    pos += 1
                if positions:
                    found_sites[enzyme] = positions
            
            if found_sites:
                get_logger().warning(f"Found restriction sites: {found_sites}")
            else:
                get_logger().debug("No restriction sites found")
            
            return found_sites
            
        except Exception as e:
            get_logger().error(f"Restriction site check failed: {type(e).__name__}: {e}")
            return {}


class PromoterRBSDesigner:
    """
    Designs promoters and RBS sequences for optimal gene expression.
    """
    
    # Promoter libraries per organism
    PROMOTER_LIBRARY = {
        "E. coli": [
            {"name": "Ptac", "strength": 0.9, "constitutive": False, "inducer": "IPTG"},
            {"name": "Plac", "strength": 0.5, "constitutive": False, "inducer": "IPTG"},
            {"name": "Ptrc", "strength": 0.85, "constitutive": False, "inducer": "IPTG"},
            {"name": "PBAD", "strength": 0.7, "constitutive": False, "inducer": "arabinose"},
            {"name": "PT7", "strength": 0.95, "constitutive": False, "inducer": "IPTG (T7 polymerase)"},
            {"name": "J23100", "strength": 0.8, "constitutive": True, "inducer": None},
            {"name": "J23119", "strength": 0.3, "constitutive": True, "inducer": None},
        ],
        "S. cerevisiae": [
            {"name": "PGAL1", "strength": 0.95, "constitutive": False, "inducer": "galactose"},
            {"name": "PGAL10", "strength": 0.85, "constitutive": False, "inducer": "galactose"},
            {"name": "PTEF1", "strength": 0.9, "constitutive": True, "inducer": None},
            {"name": "PGPD", "strength": 0.8, "constitutive": True, "inducer": None},
            {"name": "PCYC1", "strength": 0.4, "constitutive": True, "inducer": None},
            {"name": "PADH1", "strength": 0.7, "constitutive": True, "inducer": None},
        ],
        "B. subtilis": [
            {"name": "Pveg", "strength": 0.85, "constitutive": True, "inducer": None},
            {"name": "PaprE", "strength": 0.9, "constitutive": False, "inducer": "xylose"},
            {"name": "PxylA", "strength": 0.8, "constitutive": False, "inducer": "xylose"},
            {"name": "Pspank", "strength": 0.75, "constitutive": False, "inducer": "IPTG"},
            {"name": "Phyperspank", "strength": 0.85, "constitutive": False, "inducer": "IPTG"},
        ],
        "C. glutamicum": [
            {"name": "Ptac", "strength": 0.8, "constitutive": False, "inducer": "IPTG"},
            {"name": "Psod", "strength": 0.7, "constitutive": True, "inducer": None},
            {"name": "Pgap", "strength": 0.85, "constitutive": True, "inducer": None},
            {"name": "Ptuf", "strength": 0.9, "constitutive": True, "inducer": None},
        ],
        "P. putida": [
            {"name": "Ptac", "strength": 0.8, "constitutive": False, "inducer": "IPTG"},
            {"name": "Pbad", "strength": 0.7, "constitutive": False, "inducer": "arabinose"},
            {"name": "Pm", "strength": 0.85, "constitutive": False, "inducer": "benzoate"},
            {"name": "Pr", "strength": 0.75, "constitutive": False, "inducer": "benzoate"},
            {"name": "Ptet", "strength": 0.65, "constitutive": False, "inducer": "aTc"},
        ],
    }
    
    # RBS sequences with strengths
    RBS_LIBRARY = {
        "E. coli": [
            {"name": "RBS_strong", "sequence": "AGGAGG", "strength": 0.95},
            {"name": "RBS_medium", "sequence": "AGGAGA", "strength": 0.7},
            {"name": "RBS_weak", "sequence": "AGGA", "strength": 0.4},
        ],
        "S. cerevisiae": [
            {"name": "Kozak_strong", "sequence": "AAAAAC", "strength": 0.9},
            {"name": "Kozak_medium", "sequence": "AAAATC", "strength": 0.6},
        ],
        "B. subtilis": [
            {"name": "RBS_Bs_strong", "sequence": "AGGAGGU", "strength": 0.9},
            {"name": "RBS_Bs_medium", "sequence": "AGGAGG", "strength": 0.65},
        ],
        "C. glutamicum": [
            {"name": "RBS_Cg_strong", "sequence": "AGGAGG", "strength": 0.85},
            {"name": "RBS_Cg_medium", "sequence": "AGGA", "strength": 0.55},
        ],
        "P. putida": [
            {"name": "RBS_Pp_strong", "sequence": "AGGAGG", "strength": 0.88},
            {"name": "RBS_Pp_medium", "sequence": "AGGAGA", "strength": 0.6},
        ],
    }
    
    # Terminator sequences
    TERMINATORS = {
        "E. coli": "BBa_B0015",
        "S. cerevisiae": "CYC1_terminator",
        "B. subtilis": "amyE_terminator",
        "C. glutamicum": "rrnB_terminator",
        "P. putida": "BBa_B0015",
    }
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        get_logger().info("PromoterRBSDesigner initialized")
    
    def design_expression_cassette(self, gene_name: str, optimized_sequence: str,
                                   organism: str, desired_strength: float = 0.8) -> ExpressionCassette:
        """
        Design a complete expression cassette for a gene.
        
        Args:
            gene_name: Name of the gene
            optimized_sequence: Codon-optimized DNA sequence
            organism: Target organism
            desired_strength: Desired expression strength (0-1)
            
        Returns:
            ExpressionCassette with complete design
        """
        try:
            # Select promoter
            promoters = self.PROMOTER_LIBRARY.get(organism, self.PROMOTER_LIBRARY["E. coli"])
            best_promoter = min(promoters, key=lambda p: abs(p["strength"] - desired_strength))
            
            # Select RBS
            rbs_list = self.RBS_LIBRARY.get(organism, self.RBS_LIBRARY["E. coli"])
            best_rbs = min(rbs_list, key=lambda r: abs(r["strength"] - desired_strength))
            
            # Get terminator
            terminator = self.TERMINATORS.get(organism, "BBa_B0015")
            
            # Calculate CAI
            codon_optimizer = CodonOptimizer(self.config)
            cai = codon_optimizer.calculate_cai(optimized_sequence, organism)
            
            # Predict TPM (Transcripts Per Million)
            predicted_tpm = self.predict_expression_level(
                best_promoter["strength"], 
                best_rbs["strength"],
                cai
            )
            
            cassette = ExpressionCassette(
                cassette_id=f"CASSETTE_{uuid.uuid4().hex[:8]}",
                promoter=best_promoter["name"],
                promoter_strength=best_promoter["strength"],
                rbs=best_rbs["name"],
                rbs_strength=best_rbs["strength"],
                gene_name=gene_name,
                optimized_sequence=optimized_sequence,
                terminator=terminator,
                predicted_tpm=predicted_tpm,
                cai=cai
            )
            
            get_logger().info(f"Designed cassette {cassette.cassette_id} for {gene_name}")
            return cassette
            
        except Exception as e:
            get_logger().error(f"Cassette design failed for {gene_name}: {type(e).__name__}: {e}")
            raise ModelInferenceError(f"Expression cassette design failed: {e}")
    
    def predict_expression_level(self, promoter_strength: float, rbs_strength: float,
                                 cai: float) -> float:
        """
        Predict expression level in TPM (Transcripts Per Million).
        
        Args:
            promoter_strength: Promoter strength (0-1)
            rbs_strength: RBS strength (0-1)
            cai: Codon Adaptation Index (0-1)
            
        Returns:
            Predicted TPM value
        """
        try:
            # Simple model: TPM = base × promoter × RBS × CAI
            base_tpm = 1000  # Baseline expression
            predicted_tpm = base_tpm * promoter_strength * rbs_strength * cai
            
            # Add some biological noise
            noise = random.gauss(0, 0.1)
            predicted_tpm *= (1 + noise)
            
            predicted_tpm = max(0, min(10000, predicted_tpm))  # Clamp to realistic range
            
            get_logger().debug(f"Predicted TPM: {predicted_tpm:.1f} (promoter={promoter_strength}, RBS={rbs_strength}, CAI={cai})")
            return predicted_tpm
            
        except Exception as e:
            get_logger().error(f"Expression prediction failed: {type(e).__name__}: {e}")
            return 500.0  # Neutral fallback


class PathwayAIEngine:
    """
    Main orchestrator for Stage 2: integrates retrosynthesis, enzyme selection,
    codon optimization, and expression design.
    """
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.retrosynthesis_engine = RetrosynthesisEngine(config)
        self.enzyme_selector = EnzymeSelector(config)
        self.codon_optimizer = CodonOptimizer(config)
        self.promoter_designer = PromoterRBSDesigner(config)
        get_logger().info("PathwayAIEngine initialized")
    
    def run_stage_2(self, stage_1_output: Dict[str, Any]) -> Dict[str, Any]:
        """
        Execute Stage 2 pipeline: pathway prediction and design.
        
        Args:
            stage_1_output: Validated Stage 1 output JSON
            
        Returns:
            Stage 2 output JSON
        """
        import time
        start_time = time.time()
        
        try:
            # Validate input
            get_logger().info(f"=== STAGE 2 START === pipeline_id={stage_1_output.get('pipeline_id')}")
            get_logger().debug(f"Input JSON received: pipeline_id={stage_1_output.get('pipeline_id')}")
            
            validate_stage_output(stage_1_output, 'stage_1_output')
            get_logger().info("Schema validation: PASSED for stage_1_input_schema")
            
            # Extract information from Stage 1
            organism_info = stage_1_output["organism"]
            molecule_info = stage_1_output["target_molecule"]
            genomic_data = stage_1_output["genomic_data"]
            
            organism_name = organism_info["name"]
            target_name = molecule_info["name"]
            target_smiles = molecule_info.get("smiles", "")
            
            get_logger().info(f"Running pathway prediction for {organism_name} → {target_name}")
            
            # Define starting metabolites based on central metabolism
            starting_metabolites = [
                "glucose", "glucose-6-phosphate", "pyruvate", "acetyl-CoA",
                "PEP", "oxaloacetate", "alpha-ketoglutarate", "chorismate"
            ]
            
            # Run retrosynthesis
            pathways_raw = self.retrosynthesis_engine.run_mcts(
                target_molecule=target_name,
                starting_metabolites=starting_metabolites,
                max_iterations=100
            )
            
            get_logger().info(f"Retrosynthesis found {len(pathways_raw)} raw pathway candidates")
            
            # Rank pathways
            ranked_pathways = self.retrosynthesis_engine.rank_pathways(
                pathways=pathways_raw,
                organism=organism_name
            )
            
            # Process pathways - keep track of raw pathways alongside ranked candidates
            pathway_candidates_json = []
            all_gene_modifications = {
                "knockouts": [],
                "overexpressions": [],
                "heterologous_insertions": []
            }
            codon_optimized_sequences = {}
            
            # Use index-based iteration to maintain correspondence
            for idx, pathway in enumerate(ranked_pathways[:5]):  # Top 5 pathways
                steps_json = []
                
                # Get the corresponding raw pathway by index
                if idx < len(pathways_raw):
                    raw_pathway = pathways_raw[idx]
                else:
                    raw_pathway = []
                
                for rule_idx, rule in enumerate(raw_pathway):
                    # Find enzymes for this reaction
                    enzyme_candidates = self.enzyme_selector.query_brenda_simulated(
                        enzyme_name=rule.name,
                        ec_number=rule.ec_number
                    )
                    
                    # Filter by organism
                    filtered_enzymes = self.enzyme_selector.filter_by_organism_constraints(
                        enzyme_candidates, organism_name
                    )
                    
                    # Select best enzyme
                    best_enzyme = max(filtered_enzymes, key=lambda e: (
                        self.enzyme_selector.predict_host_compatibility(e, organism_name) * 0.4 +
                        self.enzyme_selector.score_esm2_embedding(e) * 0.3 +
                        e.kcat_per_sec * 0.2 +
                        e.stability_score * 0.1
                    ))
                    
                    # Determine if heterologous
                    is_heterologous = best_enzyme.source_organism != organism_name
                    
                    # Generate simulated protein sequence for codon optimization
                    protein_seq = "M" + "".join(random.choice("ARNDCEQGHILKMFPSTWYV") for _ in range(best_enzyme.sequence_length - 1))
                    
                    # Codon optimize
                    optimized_dna = self.codon_optimizer.optimize_sequence(protein_seq, organism_name)
                    cai = self.codon_optimizer.calculate_cai(optimized_dna, organism_name)
                    
                    # Check restriction sites
                    restriction_sites = self.codon_optimizer.check_restriction_sites(optimized_dna)
                    if restriction_sites:
                        get_logger().warning(f"Restriction sites found in {best_enzyme.gene_name}: {list(restriction_sites.keys())}")
                    
                    # Store codon-optimized sequence
                    codon_optimized_sequences[best_enzyme.gene_name] = optimized_dna
                    
                    # Create step JSON
                    step = {
                        "step_number": rule_idx + 1,
                        "reaction_id": rule.rule_id,
                        "enzyme_name": best_enzyme.name,
                        "gene_name": best_enzyme.gene_name,
                        "ec_number": best_enzyme.ec_number,
                        "substrate": rule.substrate_pattern,
                        "product": rule.product_pattern,
                        "delta_g_kj_per_mol": rule.delta_g_kj_per_mol,
                        "kcat_per_sec": best_enzyme.kcat_per_sec,
                        "km_mm": best_enzyme.km_mm,
                        "is_heterologous": is_heterologous,
                        "source_organism": best_enzyme.source_organism
                    }
                    steps_json.append(step)
                    
                    # Track gene modifications
                    if is_heterologous:
                        all_gene_modifications["heterologous_insertions"].append(best_enzyme.gene_name)
                    else:
                        all_gene_modifications["overexpressions"].append(best_enzyme.gene_name)
                
                # Design expression cassettes for each gene
                for gene_name, seq in codon_optimized_sequences.items():
                    cassette = self.promoter_designer.design_expression_cassette(
                        gene_name=gene_name,
                        optimized_sequence=seq,
                        organism=organism_name,
                        desired_strength=0.8
                    )
                    get_logger().debug(f"Designed cassette for {gene_name}: TPM={cassette.predicted_tpm:.1f}")
                
                pathway_json = {
                    "pathway_id": pathway.pathway_id,
                    "rank": pathway.rank,
                    "pathway_name": pathway.pathway_name,
                    "steps": steps_json,
                    "total_steps": pathway.total_steps,
                    "predicted_yield_mol_per_mol": round(pathway.predicted_yield_mol_per_mol, 4),
                    "thermodynamic_feasibility_score": round(pathway.thermodynamic_feasibility_score, 4),
                    "gnn_viability_score": round(pathway.gnn_viability_score, 4),
                    "host_compatibility_score": round(pathway.host_compatibility_score, 4)
                }
                pathway_candidates_json.append(pathway_json)
            
            # Suggest knockouts based on target molecule
            suggested_knockouts = self._suggest_knockouts(target_name, organism_name)
            all_gene_modifications["knockouts"].extend(suggested_knockouts)
            
            # Build Stage 2 output
            stage_2_output = {
                "pipeline_id": stage_1_output["pipeline_id"],
                "stage_1_output": stage_1_output,
                "pathway_candidates": pathway_candidates_json,
                "gene_modifications": all_gene_modifications,
                "codon_optimized_sequences": codon_optimized_sequences,
                "stage_2_status": "PASS"
            }
            
            # Validate output
            validate_stage_output(stage_2_output, 'stage_2_output')
            get_logger().info("Schema validation: PASSED for stage_2_output_schema")
            
            get_logger().info(f"Stage 2 output: {len(pathway_candidates_json)} pathways, "
                       f"{len(all_gene_modifications['knockouts'])} knockouts, "
                       f"{len(all_gene_modifications['heterologous_insertions'])} insertions")
            
            # Save stage summary
            import time
            duration = time.time() - start_time
            summary = {
                "pipeline_id": stage_1_output["pipeline_id"],
                "stage": 2,
                "status": "PASS",
                "duration_seconds": round(duration, 3),
                "organism": organism_name,
                "target_molecule": target_name,
                "pathway_candidates_count": len(pathway_candidates_json),
                "knockouts_count": len(all_gene_modifications['knockouts']),
                "overexpressions_count": len(all_gene_modifications['overexpressions']),
                "heterologous_insertions_count": len(all_gene_modifications['heterologous_insertions']),
                "codon_optimized_sequences_count": len(codon_optimized_sequences)
            }
            get_logger().save_stage_summary(summary)
            
            get_logger().log_stage_complete(
                duration_seconds=duration,
                status="PASS",
                output_summary=summary
            )
            
            return stage_2_output
            
        except Exception as e:
            get_logger().error(f"Exception in run_stage_2: {type(e).__name__}: {e}")
            get_logger().error(f"Input JSON that caused error: {json.dumps(stage_1_output, indent=2)}")
            
            # Attempt fallback
            get_logger().warning("Attempting fallback: returning minimal pathway design")
            fallback_output = self._fallback_stage_2(stage_1_output)
            
            # Save error summary
            import time
            summary = {
                "pipeline_id": stage_1_output.get("pipeline_id", "unknown"),
                "stage": 2,
                "status": "FAIL",
                "error_type": type(e).__name__,
                "error_message": str(e)
            }
            get_logger().save_stage_summary(summary)
            
            return fallback_output
    
    def _suggest_knockouts(self, target_molecule: str, organism: str) -> List[str]:
        """Suggest gene knockouts to improve production."""
        knockout_suggestions = {
            "lycopene": ["idi", "ispA"],  # Redirect FPP from native pathways
            "vanillin": ["pheA", "tyrA"],  # Block competing aromatic pathways
            "lysine": ["ldhA", "adhE", "frdA"],  # Reduce byproduct formation
            "riboflavin": ["ribG"],  # Feedback resistance
            "pha": ["fadR", "arcA"],  # Enhance precursor supply
        }
        return knockout_suggestions.get(target_molecule.lower(), [])
    
    def _fallback_stage_2(self, stage_1_output: Dict[str, Any]) -> Dict[str, Any]:
        """Fallback Stage 2 output when main pipeline fails."""
        get_logger().warning("Using fallback Stage 2 output")
        
        organism_name = stage_1_output["organism"]["name"]
        target_name = stage_1_output["target_molecule"]["name"]
        
        return {
            "pipeline_id": stage_1_output["pipeline_id"],
            "stage_1_output": stage_1_output,
            "pathway_candidates": [{
                "pathway_id": f"PATH_FALLBACK_{uuid.uuid4().hex[:8]}",
                "rank": 1,
                "pathway_name": f"Fallback pathway to {target_name}",
                "steps": [],
                "total_steps": 0,
                "predicted_yield_mol_per_mol": 0.1,
                "thermodynamic_feasibility_score": 0.5,
                "gnn_viability_score": 0.5,
                "host_compatibility_score": 0.7
            }],
            "gene_modifications": {
                "knockouts": [],
                "overexpressions": [],
                "heterologous_insertions": []
            },
            "codon_optimized_sequences": {},
            "stage_2_status": "WARN"
        }


def run_stage_2(stage_1_output: Dict[str, Any], config: Optional[PipelineConfig] = None) -> Dict[str, Any]:
    """
    Main entry point for Stage 2 execution.
    
    Args:
        stage_1_output: Validated Stage 1 output JSON
        config: Pipeline configuration (optional, uses defaults)
        
    Returns:
        Stage 2 output JSON
    """
    if config is None:
        config = PipelineConfig()
    
    engine = PathwayAIEngine(config)
    return engine.run_stage_2(stage_1_output)


if __name__ == "__main__":
    # Test Stage 2 with mock Stage 1 output
    import json
    from datetime import datetime
    
    print("Testing Stage 2: Pathway Prediction + AI Engine")
    print("=" * 60)
    
    # Create mock Stage 1 output
    mock_stage_1 = {
        "pipeline_id": str(uuid.uuid4()),
        "timestamp": datetime.now().isoformat(),
        "organism": {
            "name": "E. coli",
            "strain": "K-12 MG1655",
            "gem_model_id": "iJO1366",
            "doubling_time_min": 20.0,
            "optimal_ph": 7.2,
            "optimal_temp_c": 37.0,
            "gram_stain": "negative"
        },
        "target_molecule": {
            "name": "lycopene",
            "smiles": "CC(C=CCCC(C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC1)C=CCCC(C)=CC=CC2=C(C)CC2",
            "chebi_id": "CHEBI:15947",
            "target_titer_g_per_l": 5.0,
            "target_yield_mol_per_mol": 0.3
        },
        "genomic_data": {
            "total_genes": 4145,
            "essential_genes": ["dnaA", "dnaN", "gyrA", "gyrB", "rpoB", "rpoC"],
            "available_promoters": ["Ptac", "Plac", "Ptrc", "PBAD", "PT7"],
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
    
    try:
        # Run Stage 2
        stage_2_output = run_stage_2(mock_stage_1)
        
        print("\nStage 2 Output Summary:")
        print(f"  Pipeline ID: {stage_2_output['pipeline_id']}")
        print(f"  Status: {stage_2_output['stage_2_status']}")
        print(f"  Pathway candidates: {len(stage_2_output['pathway_candidates'])}")
        print(f"  Knockouts: {len(stage_2_output['gene_modifications']['knockouts'])}")
        print(f"  Overexpressions: {len(stage_2_output['gene_modifications']['overexpressions'])}")
        print(f"  Heterologous insertions: {len(stage_2_output['gene_modifications']['heterologous_insertions'])}")
        print(f"  Codon-optimized sequences: {len(stage_2_output['codon_optimized_sequences'])}")
        
        if stage_2_output['pathway_candidates']:
            best_pathway = stage_2_output['pathway_candidates'][0]
            print(f"\nBest pathway: {best_pathway['pathway_name']}")
            print(f"  Steps: {best_pathway['total_steps']}")
            print(f"  Predicted yield: {best_pathway['predicted_yield_mol_per_mol']:.3f}")
            print(f"  Thermodynamic score: {best_pathway['thermodynamic_feasibility_score']:.3f}")
        
        print("\n✓ Stage 2 test completed successfully!")
        
    except Exception as e:
        print(f"\n✗ Stage 2 test failed: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
