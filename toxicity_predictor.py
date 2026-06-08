"""
Toxicity Predictor - Predict toxicity of pathway intermediates.

This module provides toxicity assessment for metabolic pathway intermediates
to identify potential bottlenecks and risks in strain design.
"""

import logging
import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any
import json

from logger_setup import setup_logger
from exceptions import PipelineError
from pipeline_config import PipelineConfig

logger = setup_logger("STAGE_3", "toxicity_predictor")


@dataclass
class ToxicityAssessment:
    """Dataclass representing toxicity assessment results."""
    
    metabolite_name: str
    overall_toxicity_score: float
    ros_score: float
    membrane_disruption_score: float
    enzyme_inhibition_score: float
    risk_level: str  # LOW, MEDIUM, HIGH
    flagged: bool
    recommendations: List[str] = field(default_factory=list)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "metabolite_name": self.metabolite_name,
            "overall_toxicity_score": self.overall_toxicity_score,
            "ros_score": self.ros_score,
            "membrane_disruption_score": self.membrane_disruption_score,
            "enzyme_inhibition_score": self.enzyme_inhibition_score,
            "risk_level": self.risk_level,
            "flagged": self.flagged,
            "recommendations": self.recommendations
        }


class ToxicityPredictor:
    """
    Toxicity Prediction Engine for metabolic pathway intermediates.
    
    Assesses multiple toxicity mechanisms:
    - Reactive oxygen species (ROS) generation
    - Membrane disruption
    - Enzyme inhibition
    - General cellular stress
    """
    
    # Known toxic intermediates (simulated database)
    KNOWN_TOXIC_METABOLITES = {
        "acetaldehyde": {"ros": 0.7, "membrane": 0.5, "enzyme": 0.6},
        "formaldehyde": {"ros": 0.9, "membrane": 0.8, "enzyme": 0.9},
        "methylglyoxal": {"ros": 0.8, "membrane": 0.4, "enzyme": 0.7},
        "farnesyl_pyrophosphate": {"ros": 0.3, "membrane": 0.6, "enzyme": 0.4},
        "geranylgeranyl_pyrophosphate": {"ros": 0.3, "membrane": 0.7, "enzyme": 0.4},
        "protocatechuate": {"ros": 0.5, "membrane": 0.3, "enzyme": 0.4},
        "catechol": {"ros": 0.8, "membrane": 0.4, "enzyme": 0.6},
        "vanillin": {"ros": 0.4, "membrane": 0.5, "enzyme": 0.3},
        "p-coumarate": {"ros": 0.3, "membrane": 0.4, "enzyme": 0.3},
        "caffeate": {"ros": 0.4, "membrane": 0.4, "enzyme": 0.3},
        "lysine": {"ros": 0.1, "membrane": 0.1, "enzyme": 0.1},
        "glutamate": {"ros": 0.1, "membrane": 0.1, "enzyme": 0.1},
        "lycopene": {"ros": 0.2, "membrane": 0.3, "enzyme": 0.1},
        "riboflavin": {"ros": 0.2, "membrane": 0.2, "enzyme": 0.1},
        "acetate": {"ros": 0.3, "membrane": 0.2, "enzyme": 0.3},
        "ethanol": {"ros": 0.4, "membrane": 0.6, "enzyme": 0.4},
        "butyrate": {"ros": 0.3, "membrane": 0.4, "enzyme": 0.3},
        "pyruvate": {"ros": 0.2, "membrane": 0.1, "enzyme": 0.2},
        "glucose": {"ros": 0.1, "membrane": 0.1, "enzyme": 0.1},
        "g6p": {"ros": 0.1, "membrane": 0.1, "enzyme": 0.1},
        "f6p": {"ros": 0.1, "membrane": 0.1, "enzyme": 0.1}
    }
    
    # Organism-specific toxicity thresholds
    ORGANISM_THRESHOLDS = {
        "ecoli": {
            "ros_tolerance": 0.7,
            "membrane_tolerance": 0.6,
            "enzyme_tolerance": 0.7,
            "overall_tolerance": 0.5
        },
        "scerevisiae": {
            "ros_tolerance": 0.6,
            "membrane_tolerance": 0.7,
            "enzyme_tolerance": 0.6,
            "overall_tolerance": 0.5
        },
        "bsubtilis": {
            "ros_tolerance": 0.65,
            "membrane_tolerance": 0.55,
            "enzyme_tolerance": 0.65,
            "overall_tolerance": 0.45
        },
        "cglutamicum": {
            "ros_tolerance": 0.7,
            "membrane_tolerance": 0.6,
            "enzyme_tolerance": 0.7,
            "overall_tolerance": 0.5
        },
        "pputida": {
            "ros_tolerance": 0.75,
            "membrane_tolerance": 0.65,
            "enzyme_tolerance": 0.75,
            "overall_tolerance": 0.55
        }
    }
    
    def __init__(self, config: Optional[PipelineConfig] = None):
        """
        Initialize Toxicity Predictor.
        
        Args:
            config: Pipeline configuration object
        """
        self.config = config or PipelineConfig()
        self.logger = setup_logger("STAGE_3", "toxicity_predictor")
        
    def predict_intermediate_toxicity(
        self,
        metabolite_name: str,
        organism: str,
        concentration_mm: Optional[float] = None
    ) -> ToxicityAssessment:
        """
        Predict toxicity of a single metabolite.
        
        Args:
            metabolite_name: Name of the metabolite
            organism: Host organism identifier
            concentration_mm: Optional concentration in mM
            
        Returns:
            ToxicityAssessment with detailed scores
        """
        self.logger.debug(f"Predicting toxicity for {metabolite_name} in {organism}")
        
        # Get base toxicity scores
        if metabolite_name.lower() in self.KNOWN_TOXIC_METABOLITES:
            base_scores = self.KNOWN_TOXIC_METABOLITES[metabolite_name.lower()]
            ros_score = base_scores["ros"]
            membrane_score = base_scores["membrane"]
            enzyme_score = base_scores["enzyme"]
        else:
            # Predict based on chemical properties (simulated)
            ros_score, membrane_score, enzyme_score = self._predict_from_structure(metabolite_name)
        
        # Adjust for concentration if provided
        if concentration_mm is not None:
            conc_factor = min(concentration_mm / 10.0, 2.0)  # Cap at 2x
            ros_score *= (0.5 + 0.5 * conc_factor)
            membrane_score *= (0.5 + 0.5 * conc_factor)
            enzyme_score *= (0.5 + 0.5 * conc_factor)
        
        # Clamp to [0, 1]
        ros_score = min(1.0, max(0.0, ros_score))
        membrane_score = min(1.0, max(0.0, membrane_score))
        enzyme_score = min(1.0, max(0.0, enzyme_score))
        
        # Calculate overall score (weighted average)
        overall_score = 0.4 * ros_score + 0.35 * membrane_score + 0.25 * enzyme_score
        
        # Determine risk level based on organism thresholds
        thresholds = self.ORGANISM_THRESHOLDS.get(organism, self.ORGANISM_THRESHOLDS["ecoli"])
        
        if overall_score > thresholds["overall_tolerance"]:
            risk_level = "HIGH"
            flagged = True
        elif overall_score > thresholds["overall_tolerance"] * 0.7:
            risk_level = "MEDIUM"
            flagged = False
        else:
            risk_level = "LOW"
            flagged = False
        
        # Generate recommendations
        recommendations = self._generate_recommendations(
            metabolite_name, ros_score, membrane_score, enzyme_score, risk_level
        )
        
        assessment = ToxicityAssessment(
            metabolite_name=metabolite_name,
            overall_toxicity_score=overall_score,
            ros_score=ros_score,
            membrane_disruption_score=membrane_score,
            enzyme_inhibition_score=enzyme_score,
            risk_level=risk_level,
            flagged=flagged,
            recommendations=recommendations
        )
        
        self.logger.debug(f"Toxicity assessment: {assessment.to_dict()}")
        return assessment
    
    def _predict_from_structure(self, metabolite_name: str) -> tuple:
        """
        Predict toxicity from chemical structure (simulated).
        
        In a real implementation, this would use RDKit or similar.
        
        Args:
            metabolite_name: Name of metabolite
            
        Returns:
            Tuple of (ros_score, membrane_score, enzyme_score)
        """
        # Simple heuristics based on name patterns
        ros_score = 0.2
        membrane_score = 0.2
        enzyme_score = 0.2
        
        name_lower = metabolite_name.lower()
        
        # Aldehydes tend to be reactive
        if "aldehyde" in name_lower or "al" in name_lower[-2:]:
            ros_score += 0.3
            enzyme_score += 0.2
        
        # Long hydrophobic chains disrupt membranes
        if "farnesyl" in name_lower or "geranyl" in name_lower or "phytyl" in name_lower:
            membrane_score += 0.4
        
        # Phenolic compounds can be toxic
        if "catechol" in name_lower or "phenol" in name_lower or "vanill" in name_lower:
            ros_score += 0.2
            membrane_score += 0.1
        
        # Organic acids
        if "ate" in name_lower[-3:] or "ic acid" in name_lower:
            enzyme_score += 0.1
        
        # Amino acids are generally safe
        if any(aa in name_lower for aa in ["lysine", "glutamate", "threonine", "tryptophan"]):
            ros_score = 0.1
            membrane_score = 0.1
            enzyme_score = 0.1
        
        return (
            min(1.0, ros_score),
            min(1.0, membrane_score),
            min(1.0, enzyme_score)
        )
    
    def _generate_recommendations(
        self,
        metabolite_name: str,
        ros_score: float,
        membrane_score: float,
        enzyme_score: float,
        risk_level: str
    ) -> List[str]:
        """Generate mitigation recommendations based on toxicity profile."""
        recommendations = []
        
        if risk_level == "HIGH":
            recommendations.append("Consider dynamic regulation to limit accumulation")
            recommendations.append("Implement efflux pump overexpression")
        
        if ros_score > 0.6:
            recommendations.append("Overexpress antioxidant enzymes (SOD, catalase)")
            recommendations.append("Consider NADPH balancing strategies")
        
        if membrane_score > 0.6:
            recommendations.append("Modify membrane composition (cyclopropane fatty acids)")
            recommendations.append("Consider two-phase fermentation with organic overlay")
        
        if enzyme_score > 0.6:
            recommendations.append("Use feedback-resistant enzyme variants")
            recommendations.append("Consider enzyme compartmentalization")
        
        return recommendations
    
    def calculate_ros_score(
        self,
        metabolites: List[str],
        organism: str
    ) -> Dict[str, float]:
        """
        Calculate ROS generation scores for multiple metabolites.
        
        Args:
            metabolites: List of metabolite names
            organism: Host organism
            
        Returns:
            Dict mapping metabolite -> ROS score
        """
        self.logger.debug(f"Calculating ROS scores for {len(metabolites)} metabolites")
        
        ros_scores = {}
        for met in metabolites:
            assessment = self.predict_intermediate_toxicity(met, organism)
            ros_scores[met] = assessment.ros_score
        
        return ros_scores
    
    def check_membrane_disruption(
        self,
        metabolites: List[str],
        organism: str
    ) -> Dict[str, float]:
        """
        Check membrane disruption potential for metabolites.
        
        Args:
            metabolites: List of metabolite names
            organism: Host organism
            
        Returns:
            Dict mapping metabolite -> membrane disruption score
        """
        self.logger.debug(f"Checking membrane disruption for {len(metabolites)} metabolites")
        
        membrane_scores = {}
        for met in metabolites:
            assessment = self.predict_intermediate_toxicity(met, organism)
            membrane_scores[met] = assessment.membrane_disruption_score
        
        return membrane_scores
    
    def assess_pathway_toxicity(
        self,
        pathway_steps: List[Dict[str, Any]],
        organism: str
    ) -> Dict[str, Any]:
        """
        Assess overall toxicity risk for an entire pathway.
        
        Args:
            pathway_steps: List of pathway step dictionaries
            organism: Host organism
            
        Returns:
            Dict with overall toxicity assessment
        """
        self.logger.info(f"Assessing pathway toxicity for {len(pathway_steps)} steps in {organism}")
        
        # Extract all metabolites
        metabolites = set()
        for step in pathway_steps:
            if step.get("substrate"):
                metabolites.add(step["substrate"])
            if step.get("product"):
                metabolites.add(step["product"])
        
        # Assess each metabolite
        assessments = {}
        flagged_intermediates = []
        max_toxicity = 0.0
        
        for met in metabolites:
            assessment = self.predict_intermediate_toxicity(met, organism)
            assessments[met] = assessment.to_dict()
            
            if assessment.flagged:
                flagged_intermediates.append(met)
            
            max_toxicity = max(max_toxicity, assessment.overall_toxicity_score)
        
        # Determine overall risk
        if max_toxicity > 0.7 or len(flagged_intermediates) > 2:
            overall_risk = "HIGH"
        elif max_toxicity > 0.4 or len(flagged_intermediates) > 0:
            overall_risk = "MEDIUM"
        else:
            overall_risk = "LOW"
        
        result = {
            "intermediate_toxicity_scores": {m: a["overall_toxicity_score"] for m, a in assessments.items()},
            "overall_toxicity_risk": overall_risk,
            "flagged_intermediates": flagged_intermediates,
            "detailed_assessments": assessments,
            "max_toxicity_score": max_toxicity,
            "num_flagged": len(flagged_intermediates)
        }
        
        self.logger.info(f"Pathway toxicity assessment: risk={overall_risk}, flagged={len(flagged_intermediates)}")
        return result
    
    def organism_toxicity_thresholds(self, organism: str) -> Dict[str, float]:
        """
        Get toxicity thresholds for a specific organism.
        
        Args:
            organism: Organism identifier
            
        Returns:
            Dict with threshold values
        """
        return self.ORGANISM_THRESHOLDS.get(organism, self.ORGANISM_THRESHOLDS["ecoli"])


if __name__ == "__main__":
    print("=== Testing Toxicity Predictor ===")
    
    predictor = ToxicityPredictor()
    
    # Test individual metabolite toxicity
    test_metabolites = ["lycopene", "acetaldehyde", "vanillin", "glucose", "farnesyl_pyrophosphate"]
    
    print("\nIndividual metabolite assessments:")
    for met in test_metabolites:
        assessment = predictor.predict_intermediate_toxicity(met, "ecoli")
        print(f"  {met}: score={assessment.overall_toxicity_score:.2f}, risk={assessment.risk_level}")
    
    # Test pathway toxicity
    test_pathway = [
        {"step_number": 1, "substrate": "glucose", "product": "g6p"},
        {"step_number": 2, "substrate": "g6p", "product": "pyruvate"},
        {"step_number": 3, "substrate": "pyruvate", "product": "acetyl_coa"},
        {"step_number": 4, "substrate": "acetyl_coa", "product": "lycopene"}
    ]
    
    print("\nPathway toxicity assessment:")
    pathway_result = predictor.assess_pathway_toxicity(test_pathway, "ecoli")
    print(f"  Overall risk: {pathway_result['overall_toxicity_risk']}")
    print(f"  Flagged intermediates: {pathway_result['flagged_intermediates']}")
    print(f"  Max toxicity score: {pathway_result['max_toxicity_score']:.2f}")
    
    # Test organism thresholds
    print("\nOrganism toxicity thresholds:")
    for org in ["ecoli", "scerevisiae", "pputida"]:
        thresholds = predictor.organism_toxicity_thresholds(org)
        print(f"  {org}: overall_tolerance={thresholds['overall_tolerance']}")
    
    print("\n=== Toxicity Predictor Test Complete ===")
