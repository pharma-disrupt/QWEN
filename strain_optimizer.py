"""
Strain Optimizer - Genetic modification optimization for metabolic engineering.

This module provides strain design optimization using multi-objective optimization
to identify optimal knockouts and overexpression targets.
"""

import logging
import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any
from scipy.optimize import differential_evolution, minimize
import json

from logger_setup import setup_logger
from exceptions import PipelineError
from pipeline_config import PipelineConfig
from fba_engine import FBAEngine, FBAModel, FBAResult

logger = setup_logger("STAGE_3", "strain_optimizer")


@dataclass
class GeneticModification:
    """Dataclass representing a genetic modification."""
    
    gene_name: str
    modification_type: str  # knockout, overexpression, insertion
    expected_effect: str
    confidence_score: float
    source_organism: Optional[str] = None
    ec_number: Optional[str] = None
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "gene_name": self.gene_name,
            "modification_type": self.modification_type,
            "expected_effect": self.expected_effect,
            "confidence_score": self.confidence_score,
            "source_organism": self.source_organism,
            "ec_number": self.ec_number
        }


@dataclass
class StrainDesign:
    """Dataclass representing a complete strain design."""
    
    design_id: str
    knockouts: List[str]
    overexpressions: List[str]
    insertions: List[str]
    predicted_titer_g_per_l: float
    predicted_productivity_g_per_l_per_h: float
    metabolic_burden_score: float
    genetic_stability_score: float
    pareto_rank: int = 1
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "design_id": self.design_id,
            "knockouts": self.knockouts,
            "overexpressions": self.overexpressions,
            "insertions": self.insertions,
            "predicted_titer_g_per_l": self.predicted_titer_g_per_l,
            "predicted_productivity_g_per_l_per_h": self.predicted_productivity_g_per_l_per_h,
            "metabolic_burden_score": self.metabolic_burden_score,
            "genetic_stability_score": self.genetic_stability_score,
            "pareto_rank": self.pareto_rank
        }


class StrainOptimizer:
    """
    Strain Optimization Engine for identifying optimal genetic modifications.
    
    Uses multi-objective optimization (NSGA-III inspired) to find Pareto-optimal
    strain designs balancing production, growth, and stability.
    """
    
    # Gene essentiality scores by organism (simulated)
    ESSENTIAL_GENES = {
        "ecoli": ["dnaA", "rpoB", "gyrA", "ftsZ", "murA", "fabI"],
        "scerevisiae": ["ACT1", "TUB1", "CDC28", "RPB1", "CDC6"],
        "bsubtilis": ["dnaA", "sigA", "ftsZ", "murA", "accC"],
        "cglutamicum": ["dnaA", "ftsZ", "murA", "fabI", "accD"],
        "pputida": ["dnaA", "rpoD", "ftsZ", "murA", "fabI"]
    }
    
    # Common knockout targets for production enhancement
    KNOWN_KNOCKOUT_TARGETS = {
        "lycopene": ["pgi", "pfkA", "pykF", "ldhA", "adhE", "ackA"],
        "vanillin": ["fadR", "fadE", "fadL", "fadD", "atoDA"],
        "l_lysine": ["ldhA", "mdh", "pqo", "odgA", "hom"],
        "riboflavin": ["ribR", "purH", "guaB", "add"],
        "pha": ["fadR", "fadE", "fadL", "phaZ"]
    }
    
    # Overexpression targets
    OVEREXPRESSION_TARGETS = {
        "lycopene": ["dxs", "idi", "crtE", "crtB", "crtI"],
        "vanillin": ["aroG", "ppsA", "tktA", "talB", "vdh"],
        "l_lysine": ["dapA", "ddh", "lysC", "asd", "hom"],
        "riboflavin": ["ribA", "ribB", "ribC", "ribD", "ribE"],
        "pha": ["phaA", "phaB", "phaC"]
    }
    
    def __init__(self, config: Optional[PipelineConfig] = None):
        """
        Initialize Strain Optimizer.
        
        Args:
            config: Pipeline configuration object
        """
        self.config = config or PipelineConfig()
        self.logger = setup_logger("STAGE_3", "strain_optimizer")
        self.fba_engine = FBAEngine(config)
        
    def run_optknock_simulated(
        self,
        pathway_steps: List[Dict[str, Any]],
        organism: str,
        target_molecule: str,
        max_knockouts: int = 5
    ) -> List[str]:
        """
        Simulate OptKnock algorithm to identify knockout targets.
        
        Args:
            pathway_steps: Pathway reaction steps
            organism: Organism identifier
            target_molecule: Target product name
            max_knockouts: Maximum number of knockouts to suggest
            
        Returns:
            List of gene names to knockout
        """
        self.logger.info(f"Running OptKnock simulation for {target_molecule} in {organism}")
        self.logger.debug(f"Pathway steps: {len(pathway_steps)}, max knockouts: {max_knockouts}")
        
        # Get known knockout targets for this molecule
        known_targets = self.KNOWN_KNOCKOUT_TARGETS.get(target_molecule.lower(), [])
        essential = self.ESSENTIAL_GENES.get(organism, [])
        
        # Filter out essential genes
        candidate_knockouts = [g for g in known_targets if g not in essential]
        
        # Score each knockout based on expected benefit
        scored_knockouts = []
        for gene in candidate_knockouts:
            score = self._score_knockout(gene, organism, target_molecule)
            scored_knockouts.append((gene, score))
        
        # Sort by score and take top N
        scored_knockouts.sort(key=lambda x: x[1], reverse=True)
        selected = [g for g, s in scored_knockouts[:max_knockouts]]
        
        self.logger.info(f"Selected {len(selected)} knockout targets: {selected}")
        return selected
    
    def _score_knockout(self, gene: str, organism: str, target: str) -> float:
        """
        Score a potential knockout target.
        
        Args:
            gene: Gene name
            organism: Organism identifier
            target: Target molecule
            
        Returns:
            Score between 0 and 1 (higher = better target)
        """
        # Base score from literature knowledge
        base_scores = {
            "pgi": 0.85, "pfkA": 0.7, "pykF": 0.65,
            "ldhA": 0.9, "adhE": 0.8, "ackA": 0.75,
            "fadR": 0.8, "fadE": 0.7, "fadL": 0.65,
            "mdh": 0.75, "pqo": 0.7, "hom": 0.85
        }
        
        base = base_scores.get(gene, 0.5)
        
        # Adjust for organism
        organism_factor = {
            "ecoli": 1.0,
            "scerevisiae": 0.9,
            "bsubtilis": 0.85,
            "cglutamicum": 0.95,
            "pputida": 0.88
        }.get(organism, 0.9)
        
        return base * organism_factor
    
    def run_nsga3_optimization(
        self,
        pathway_steps: List[Dict[str, Any]],
        organism: str,
        target_molecule: str,
        population_size: int = 50,
        generations: int = 30
    ) -> List[StrainDesign]:
        """
        Run NSGA-III inspired multi-objective optimization.
        
        Objectives:
        1. Maximize product titer
        2. Maximize growth rate
        3. Minimize metabolic burden
        4. Maximize genetic stability
        
        Args:
            pathway_steps: Pathway reaction steps
            organism: Organism identifier
            target_molecule: Target product
            population_size: Population size for evolution
            generations: Number of generations
            
        Returns:
            List of Pareto-optimal strain designs
        """
        self.logger.info(f"Running NSGA-III optimization: {population_size} individuals, {generations} generations")
        
        # Get candidate modifications
        essential = self.ESSENTIAL_GENES.get(organism, [])
        ko_candidates = [g for g in self.KNOWN_KNOCKOUT_TARGETS.get(target_molecule.lower(), []) 
                        if g not in essential]
        oe_candidates = self.OVEREXPRESSION_TARGETS.get(target_molecule.lower(), [])
        
        # Encode solutions as binary vectors
        n_ko = len(ko_candidates)
        n_oe = len(oe_candidates)
        n_vars = n_ko + n_oe
        
        if n_vars == 0:
            self.logger.warning("No candidate modifications found, returning default design")
            return [self._create_default_design(organism, target_molecule)]
        
        # Initialize population
        population = np.random.randint(0, 2, size=(population_size, n_vars))
        
        # Evolution loop
        best_designs = []
        
        for gen in range(generations):
            # Evaluate objectives
            objectives = []
            for individual in population:
                obj = self._evaluate_objectives(individual, ko_candidates, oe_candidates, 
                                               organism, target_molecule, pathway_steps)
                objectives.append(obj)
            
            # Non-dominated sorting (simplified)
            ranks = self._non_dominated_sort(objectives)
            
            # Selection and variation
            if gen < generations - 1:
                population = self._evolve_population(population, ranks, objectives)
        
        # Extract Pareto front
        pareto_front_indices = [i for i, r in enumerate(ranks) if r == 1]
        
        for idx in pareto_front_indices[:10]:  # Top 10 designs
            individual = population[idx]
            design = self._individual_to_design(
                individual, ko_candidates, oe_candidates, 
                organism, target_molecule, pathway_steps, idx
            )
            best_designs.append(design)
        
        self.logger.info(f"Found {len(best_designs)} Pareto-optimal designs")
        return best_designs
    
    def _evaluate_objectives(
        self,
        individual: np.ndarray,
        ko_candidates: List[str],
        oe_candidates: List[str],
        organism: str,
        target: str,
        pathway_steps: List[Dict[str, Any]]
    ) -> np.ndarray:
        """
        Evaluate multiple objectives for an individual.
        
        Returns array of 4 objectives (to be minimized, so negate maximization objectives).
        """
        # Decode individual
        ko_mask = individual[:len(ko_candidates)]
        oe_mask = individual[len(ko_candidates):]
        
        knockouts = [ko_candidates[i] for i in range(len(ko_candidates)) if ko_mask[i] == 1]
        overexpressions = [oe_candidates[i] for i in range(len(oe_candidates)) if oe_mask[i] == 1]
        
        # Objective 1: Product titer (negate for minimization)
        titer = self._predict_titer(knockouts, overexpressions, organism, target)
        obj1 = -titer / 10.0  # Normalize
        
        # Objective 2: Growth rate (negate)
        growth = self._predict_growth(knockouts, organism)
        obj2 = -growth
        
        # Objective 3: Metabolic burden (minimize)
        burden = self.score_metabolic_burden(knockouts, overexpressions, organism)
        obj3 = burden
        
        # Objective 4: Genetic instability (minimize)
        instability = 1.0 - self.predict_genetic_stability(knockouts, overexpressions, organism)
        obj4 = instability
        
        return np.array([obj1, obj2, obj3, obj4])
    
    def _predict_titer(
        self, 
        knockouts: List[str], 
        overexpressions: List[str],
        organism: str, 
        target: str
    ) -> float:
        """Predict product titer based on modifications."""
        # Base titer
        base_titers = {
            "lycopene": 2.0,
            "vanillin": 1.5,
            "l_lysine": 5.0,
            "riboflavin": 3.0,
            "pha": 4.0
        }
        base = base_titers.get(target.lower(), 1.0)
        
        # Knockout bonus
        ko_bonus = len(knockouts) * 0.3
        
        # Overexpression bonus
        oe_bonus = len(overexpressions) * 0.4
        
        # Synergy bonus
        synergy = min(len(knockouts), len(overexpressions)) * 0.2
        
        predicted = base + ko_bonus + oe_bonus + synergy
        
        # Cap at realistic maximum
        max_titers = {
            "lycopene": 8.0,
            "vanillin": 5.0,
            "l_lysine": 20.0,
            "riboflavin": 10.0,
            "pha": 15.0
        }
        return min(predicted, max_titers.get(target.lower(), 10.0))
    
    def _predict_growth(self, knockouts: List[str], organism: str) -> float:
        """Predict growth rate after knockouts."""
        # Base growth rate
        base_rates = {
            "ecoli": 0.9,
            "scerevisiae": 0.45,
            "bsubtilis": 0.7,
            "cglutamicum": 0.55,
            "pputida": 0.6
        }
        base = base_rates.get(organism, 0.6)
        
        # Essential genes cause severe growth defect
        essential = self.ESSENTIAL_GENES.get(organism, [])
        essential_hits = len([g for g in knockouts if g in essential])
        
        # Growth reduction per knockout
        reduction = len(knockouts) * 0.05 + essential_hits * 0.3
        
        return max(0.1, base - reduction)
    
    def _non_dominated_sort(self, objectives: List[np.ndarray]) -> List[int]:
        """
        Perform non-dominated sorting (simplified NSGA-III).
        
        Returns rank for each individual (1 = Pareto front).
        """
        n = len(objectives)
        ranks = [1] * n
        
        for i in range(n):
            for j in range(n):
                if i != j:
                    # Check if j dominates i
                    if self._dominates(objectives[j], objectives[i]):
                        ranks[i] = max(ranks[i], 2)
                        break
        
        return ranks
    
    def _dominates(self, obj1: np.ndarray, obj2: np.ndarray) -> bool:
        """Check if obj1 dominates obj2 (all values <= and at least one <)."""
        all_better_or_equal = np.all(obj1 <= obj2)
        any_better = np.any(obj1 < obj2)
        return all_better_or_equal and any_better
    
    def _evolve_population(
        self,
        population: np.ndarray,
        ranks: List[int],
        objectives: List[np.ndarray]
    ) -> np.ndarray:
        """Evolve population using tournament selection and crossover."""
        n = len(population)
        new_population = np.zeros_like(population)
        
        for i in range(n):
            # Tournament selection
            candidates = np.random.choice(n, size=3, replace=False)
            best = min(candidates, key=lambda c: ranks[c])
            parent1 = population[best]
            
            candidates = np.random.choice(n, size=3, replace=False)
            best = min(candidates, key=lambda c: ranks[c])
            parent2 = population[best]
            
            # Crossover
            point = np.random.randint(1, len(parent1) - 1)
            child = np.concatenate([parent1[:point], parent2[point:]])
            
            # Mutation
            mutation_rate = 0.1
            for j in range(len(child)):
                if np.random.random() < mutation_rate:
                    child[j] = 1 - child[j]
            
            new_population[i] = child
        
        return new_population
    
    def _individual_to_design(
        self,
        individual: np.ndarray,
        ko_candidates: List[str],
        oe_candidates: List[str],
        organism: str,
        target: str,
        pathway_steps: List[Dict[str, Any]],
        design_idx: int
    ) -> StrainDesign:
        """Convert genetic encoding to StrainDesign."""
        ko_mask = individual[:len(ko_candidates)]
        oe_mask = individual[len(ko_candidates):]
        
        knockouts = [ko_candidates[i] for i in range(len(ko_candidates)) if ko_mask[i] == 1]
        overexpressions = [oe_candidates[i] for i in range(len(oe_candidates)) if oe_mask[i] == 1]
        
        titer = self._predict_titer(knockouts, overexpressions, organism, target)
        productivity = titer / 24.0  # Assume 24h fermentation
        burden = self.score_metabolic_burden(knockouts, overexpressions, organism)
        stability = self.predict_genetic_stability(knockouts, overexpressions, organism)
        
        return StrainDesign(
            design_id=f"DESIGN_{organism}_{target}_{design_idx:03d}",
            knockouts=knockouts,
            overexpressions=overexpressions,
            insertions=[],
            predicted_titer_g_per_l=titer,
            predicted_productivity_g_per_l_per_h=productivity,
            metabolic_burden_score=burden,
            genetic_stability_score=stability,
            pareto_rank=1
        )
    
    def _create_default_design(self, organism: str, target: str) -> StrainDesign:
        """Create a default strain design when no candidates found."""
        return StrainDesign(
            design_id=f"DESIGN_{organism}_{target}_DEFAULT",
            knockouts=[],
            overexpressions=self.OVEREXPRESSION_TARGETS.get(target.lower(), [])[:3],
            insertions=[],
            predicted_titer_g_per_l=1.0,
            predicted_productivity_g_per_l_per_h=0.04,
            metabolic_burden_score=0.3,
            genetic_stability_score=0.9,
            pareto_rank=1
        )
    
    def score_metabolic_burden(
        self,
        knockouts: List[str],
        overexpressions: List[str],
        organism: str
    ) -> float:
        """
        Calculate metabolic burden score.
        
        Args:
            knockouts: List of knocked out genes
            overexpressions: List of overexpressed genes
            organism: Organism identifier
            
        Returns:
            Burden score 0-1 (lower is better)
        """
        self.logger.debug(f"Scoring metabolic burden: {len(knockouts)} KO, {len(overexpressions)} OE")
        
        # Base burden from heterologous expression
        base_burden = 0.1
        
        # Each knockout adds some burden
        ko_burden = len(knockouts) * 0.05
        
        # Each overexpression adds more burden
        oe_burden = len(overexpressions) * 0.08
        
        # Organism-specific factors
        organism_factors = {
            "ecoli": 1.0,
            "scerevisiae": 0.9,
            "bsubtilis": 0.95,
            "cglutamicum": 0.92,
            "pputida": 0.88
        }
        
        total_burden = (base_burden + ko_burden + oe_burden) * organism_factors.get(organism, 1.0)
        
        return min(1.0, total_burden)
    
    def predict_genetic_stability(
        self,
        knockouts: List[str],
        overexpressions: List[str],
        organism: str
    ) -> float:
        """
        Predict genetic stability of the strain design.
        
        Args:
            knockouts: List of knocked out genes
            overexpressions: List of overexpressed genes
            organism: Organism identifier
            
        Returns:
            Stability score 0-1 (higher is better)
        """
        self.logger.debug(f"Predicting genetic stability: {len(knockouts)} KO, {len(overexpressions)} OE")
        
        # Base stability
        base_stability = 0.95
        
        # Essential gene hits reduce stability dramatically
        essential = self.ESSENTIAL_GENES.get(organism, [])
        essential_hits = len([g for g in knockouts if g in essential])
        essential_penalty = essential_hits * 0.3
        
        # Each modification reduces stability slightly
        mod_penalty = (len(knockouts) + len(overexpressions)) * 0.02
        
        # Plasmid vs chromosomal (assume chromosomal for stability)
        integration_bonus = 0.05
        
        stability = base_stability - essential_penalty - mod_penalty + integration_bonus
        
        return max(0.1, min(1.0, stability))
    
    def rank_strain_designs(
        self,
        designs: List[StrainDesign],
        weights: Optional[Dict[str, float]] = None
    ) -> List[StrainDesign]:
        """
        Rank strain designs by weighted score.
        
        Args:
            designs: List of strain designs
            weights: Optional weights for objectives
            
        Returns:
            Sorted list of designs (best first)
        """
        if weights is None:
            weights = {
                "titer": 0.4,
                "productivity": 0.3,
                "burden": 0.15,
                "stability": 0.15
            }
        
        self.logger.info(f"Ranking {len(designs)} strain designs")
        
        scored_designs = []
        for design in designs:
            # Normalize scores
            titer_score = min(design.predicted_titer_g_per_l / 10.0, 1.0)
            prod_score = min(design.predicted_productivity_g_per_l_per_h / 0.5, 1.0)
            burden_score = 1.0 - design.metabolic_burden_score
            stability_score = design.genetic_stability_score
            
            total_score = (
                weights["titer"] * titer_score +
                weights["productivity"] * prod_score +
                weights["burden"] * burden_score +
                weights["stability"] * stability_score
            )
            
            scored_designs.append((design, total_score))
        
        scored_designs.sort(key=lambda x: x[1], reverse=True)
        
        # Update pareto ranks
        for i, (design, _) in enumerate(scored_designs):
            design.pareto_rank = i + 1
        
        self.logger.info(f"Top design: {scored_designs[0][0].design_id} (score={scored_designs[0][1]:.3f})")
        
        return [d for d, _ in scored_designs]


if __name__ == "__main__":
    print("=== Testing Strain Optimizer ===")
    
    optimizer = StrainOptimizer()
    
    # Test pathway
    test_pathway = [
        {"step_number": 1, "reaction_id": "R1", "substrate": "glucose", "product": "g6p"},
        {"step_number": 2, "reaction_id": "R2", "substrate": "g6p", "product": "pyruvate"},
        {"step_number": 3, "reaction_id": "R3", "substrate": "pyruvate", "product": "lycopene"}
    ]
    
    # Test OptKnock
    knockouts = optimizer.run_optknock_simulated(test_pathway, "ecoli", "lycopene")
    print(f"\nOptKnock results: {knockouts}")
    
    # Test NSGA-III
    designs = optimizer.run_nsga3_optimization(test_pathway, "ecoli", "lycopene", 
                                               population_size=20, generations=10)
    print(f"\nFound {len(designs)} Pareto-optimal designs:")
    for d in designs[:3]:
        print(f"  {d.design_id}: titer={d.predicted_titer_g_per_l:.2f}, burden={d.metabolic_burden_score:.2f}")
    
    # Test ranking
    ranked = optimizer.rank_strain_designs(designs)
    print(f"\nTop ranked design: {ranked[0].to_dict()}")
    
    # Test burden scoring
    burden = optimizer.score_metabolic_burden(["ldhA", "adhE"], ["dxs", "crtE"], "ecoli")
    print(f"\nMetabolic burden score: {burden:.3f}")
    
    # Test stability prediction
    stability = optimizer.predict_genetic_stability(["ldhA"], ["dxs"], "ecoli")
    print(f"Genetic stability score: {stability:.3f}")
    
    print("\n=== Strain Optimizer Test Complete ===")
