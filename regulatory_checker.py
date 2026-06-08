"""
Regulatory Compliance Checker for Industrial Metabolic Pathway Pipeline.
Assesses biosafety, allergenicity, toxicity, and GRAS status.
"""

import logging
import random
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional
from datetime import datetime

# Internal imports
from pipeline_config import PipelineConfig, OrganismType, MoleculeType
from exceptions import PipelineError
import logging

def get_logger(name: str):
    """Get a standard Python logger for the given module name."""
    return logging.getLogger(name)

logger = get_logger(__name__)

@dataclass
class RegulatoryAssessment:
    """Complete regulatory assessment results."""
    biosafety_level: str  # BSL-1, BSL-2, BSL-3
    gras_status: bool
    allergenicity_risk: str  # LOW, MEDIUM, HIGH
    toxicity_risk: str  # LOW, MEDIUM, HIGH
    antibiotic_marker_required: bool
    compliance_score: float  # 0-100
    flagged_issues: List[str] = field(default_factory=list)
    recommendations: List[str] = field(default_factory=list)

class RegulatoryChecker:
    """
    Simulates regulatory compliance assessment for biotechnology products.
    Based on FDA, EMA, and international biosafety guidelines.
    """
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.logger = get_logger(__name__)
        
        # Biosafety Level classification by organism
        self.BIOSAFETY_LEVELS = {
            OrganismType.ECOLI: "BSL-1",  # K-12 strains are BSL-1
            OrganismType.SCEREVISIAE: "BSL-1",  # Baker's yeast is BSL-1
            OrganismType.BSUBTILIS: "BSL-1",  # Generally regarded as safe
            OrganismType.CGLUTAMICUM: "BSL-1",  # Food-grade organism
            OrganismType.PPUTIDA: "BSL-1"  # Environmental strain, low pathogenicity
        }
        
        # GRAS (Generally Recognized As Safe) status
        self.GRAS_STATUS = {
            OrganismType.ECOLI: False,  # Some strains used in production but not fully GRAS
            OrganismType.SCEREVISIAE: True,  # Bread/brewer's yeast is GRAS
            OrganismType.BSUBTILIS: True,  # Used in natto production, GRAS
            OrganismType.CGLUTAMICUM: True,  # Amino acid production, GRAS
            OrganismType.PPUTIDA: False  # Environmental, not food-grade
        }
        
        # Known allergenicity patterns (simplified simulation)
        self.ALLERGEN_PATTERNS = {
            'high_risk': ['tropomyosin', 'casein', 'gliadin', 'arachin'],
            'medium_risk': ['chitinase', 'thaumatin', 'defensin'],
            'low_risk': []  # Most microbial metabolites
        }
        
        # Toxicity thresholds (simulated ProTox-like assessment)
        self.TOXICITY_THRESHOLDS = {
            'ld50_oral_rat_mg_kg': {
                'LOW': 2000,    # > 2000 mg/kg = low toxicity
                'MEDIUM': 300,  # 300-2000 = medium
                'HIGH': 0       # < 300 = high toxicity
            }
        }

    def classify_biosafety_level(self, organism: OrganismType, genetic_modifications: int = 0) -> str:
        """
        Determine biosafety level based on organism and modifications.
        """
        base_bsl = self.BIOSAFETY_LEVELS.get(organism, "BSL-2")
        
        # Increase BSL if extensive genetic modification
        if genetic_modifications > 10:
            if base_bsl == "BSL-1":
                return "BSL-2"
        elif genetic_modifications > 5:
            self.logger.warning(f"Organism {organism.value} has {genetic_modifications} modifications - monitor closely")
            
        return base_bsl

    def check_gras_status(self, organism: OrganismType, product_type: str) -> Dict[str, Any]:
        """
        Check GRAS status and regulatory pathway.
        """
        is_gras = self.GRAS_STATUS.get(organism, False)
        
        regulatory_pathway = "Food Additive Petition" if not is_gras else "GRAS Notification"
        
        notes = []
        if is_gras:
            notes.append(f"{organism.value} has established GRAS status for food applications")
        else:
            notes.append(f"{organism.value} requires full safety assessment")
            
        if product_type.lower() in ['amino_acid', 'vitamin', 'flavor']:
            notes.append("Product class may qualify for expedited review")
            
        return {
            "is_gras": is_gras,
            "organism_status": "GRAS" if is_gras else "Non-GRAS",
            "regulatory_pathway": regulatory_pathway,
            "notes": notes
        }

    def predict_allergenicity(self, molecule_name: str, gene_sequences: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Simulate allergenicity prediction (like AllerTop/AllergenFP).
        """
        self.logger.debug(f"Predicting allergenicity for {molecule_name}")
        
        # Check if molecule name matches known allergens
        mol_lower = molecule_name.lower()
        
        risk_level = "LOW"
        confidence = 0.85
        details = []
        
        # Small molecule metabolites are typically non-allergenic
        if any(pattern in mol_lower for pattern in ['acid', 'ine', 'ol', 'one']):
            risk_level = "LOW"
            details.append("Small molecule metabolite - low allergenic potential")
        else:
            # Check for protein-like names
            if any(pattern in mol_lower for pattern in ['protein', 'peptide', 'enzyme']):
                risk_level = "MEDIUM"
                confidence = 0.70
                details.append("Protein/peptide structure detected - requires sequence analysis")
                
                # Simulate sequence-based prediction if sequences provided
                if gene_sequences:
                    # In real implementation: BLAST against allergen databases
                    # Here: random simulation
                    if random.random() < 0.1:  # 10% chance of flagging
                        risk_level = "HIGH"
                        details.append("Sequence similarity to known allergens detected")
        
        return {
            "risk_level": risk_level,
            "confidence": confidence,
            "details": details,
            "recommendation": "No additional testing required" if risk_level == "LOW" else "In vitro IgE binding assay recommended"
        }

    def predict_toxicity(self, molecule_name: str, smiles: Optional[str] = None) -> Dict[str, Any]:
        """
        Simulate toxicity prediction (like ProTox-3).
        Returns predicted LD50 and toxicity class.
        """
        self.logger.debug(f"Predicting toxicity for {molecule_name}")
        
        # Simulated LD50 values based on molecule type
        # Real implementation would use QSAR models
        toxicities = {
            'lycopene': {'ld50': 5000, 'class': 'LOW', 'notes': 'Natural carotenoid, very safe'},
            'vanillin': {'ld50': 2500, 'class': 'LOW', 'notes': 'Common flavor compound'},
            'l_lysine': {'ld50': 5000, 'class': 'LOW', 'notes': 'Essential amino acid'},
            'l_glutamate': {'ld50': 5000, 'class': 'LOW', 'notes': 'Common amino acid'},
            'riboflavin': {'ld50': 5000, 'class': 'LOW', 'notes': 'Vitamin B2, essential nutrient'},
            'pha': {'ld50': 5000, 'class': 'LOW', 'notes': 'Biodegradable polymer, biocompatible'},
            'artemisinic_acid': {'ld50': 1500, 'class': 'MEDIUM', 'notes': 'Bioactive compound, monitor dosage'}
        }
        
        mol_key = molecule_name.lower().replace(" ", "_").replace("-", "_")
        
        if mol_key in toxicities:
            result = toxicities[mol_key]
        else:
            # Default prediction for unknown molecules
            result = {
                'ld50': 2000,
                'class': 'LOW',
                'notes': 'Predicted low toxicity based on structural analogs'
            }
        
        return {
            "predicted_ld50_oral_rat_mg_kg": result['ld50'],
            "toxicity_class": result['class'],
            "epa_toxicity_category": "IV" if result['class'] == 'LOW' else ("III" if result['class'] == 'MEDIUM' else "II"),
            "notes": result['notes'],
            "recommendation": "Standard safety testing" if result['class'] == 'LOW' else "Extended toxicology studies recommended"
        }

    def check_antibiotic_marker_removal(self, organism: OrganismType, 
                                        intended_use: str = "food") -> Dict[str, Any]:
        """
        Assess requirements for antibiotic marker removal.
        """
        self.logger.info(f"Checking antibiotic marker requirements for {intended_use} use")
        
        requirements = {
            'food': {
                'removal_required': True,
                'rationale': 'FDA/EFSA guidelines prohibit antibiotic resistance markers in food-producing organisms',
                'alternatives': ['Markerless deletion systems', 'Auxotrophic complementation', 'Site-specific recombination']
            },
            'pharma': {
                'removal_required': True,
                'rationale': 'EMA/FDA require removal or justification for pharmaceutical production strains',
                'alternatives': ['Markerless systems', 'Well-characterized markers with safety data']
            },
            'industrial_non_food': {
                'removal_required': False,
                'rationale': 'May be acceptable for contained industrial processes with proper waste treatment',
                'alternatives': ['Standard selection markers acceptable with containment']
            }
        }
        
        req = requirements.get(intended_use.lower(), requirements['industrial_non_food'])
        
        return {
            "removal_required": req['removal_required'],
            "regulatory_rationale": req['rationale'],
            "recommended_alternatives": req['alternatives'],
            "compliance_note": f"For {organism.value} in {intended_use} applications"
        }

    def generate_compliance_report(self, 
                                   organism: OrganismType,
                                   molecule_name: str,
                                   stage_4_json: Dict[str, Any],
                                   num_genetic_modifications: int = 5,
                                   intended_use: str = "food") -> RegulatoryAssessment:
        """
        Generate comprehensive regulatory compliance report.
        """
        self.logger.info("=== GENERATING REGULATORY COMPLIANCE REPORT ===")
        
        flagged_issues = []
        recommendations = []
        
        # 1. Biosafety Level
        bsl = self.classify_biosafety_level(organism, num_genetic_modifications)
        self.logger.info(f"Biosafety Level: {bsl}")
        
        # 2. GRAS Status
        gras_info = self.check_gras_status(organism, molecule_name)
        if not gras_info['is_gras']:
            flagged_issues.append(f"Organism {organism.value} is not GRAS - full safety assessment required")
            recommendations.append(gras_info['regulatory_pathway'])
        
        # 3. Allergenicity
        allergy_result = self.predict_allergenicity(molecule_name)
        if allergy_result['risk_level'] in ['MEDIUM', 'HIGH']:
            flagged_issues.append(f"Allergenicity concern: {allergy_result['details']}")
            recommendations.append(allergy_result['recommendation'])
        
        # 4. Toxicity
        tox_result = self.predict_toxicity(molecule_name)
        if tox_result['toxicity_class'] in ['MEDIUM', 'HIGH']:
            flagged_issues.append(f"Toxicity concern: LD50 = {tox_result['predicted_ld50_oral_rat_mg_kg']} mg/kg")
            recommendations.append(tox_result['recommendation'])
        
        # 5. Antibiotic Markers
        abx_result = self.check_antibiotic_marker_removal(organism, intended_use)
        if abx_result['removal_required']:
            recommendations.append(f"Implement antibiotic marker removal: {abx_result['recommended_alternatives'][0]}")
        
        # Calculate compliance score (0-100)
        score = 100
        if not gras_info['is_gras']:
            score -= 20
        if allergy_result['risk_level'] == 'MEDIUM':
            score -= 15
        elif allergy_result['risk_level'] == 'HIGH':
            score -= 30
        if tox_result['toxicity_class'] == 'MEDIUM':
            score -= 15
        elif tox_result['toxicity_class'] == 'HIGH':
            score -= 30
        if bsl == "BSL-2":
            score -= 10
            
        score = max(0, score)
        
        # Determine overall risks
        overall_allergy = allergy_result['risk_level']
        overall_tox = tox_result['toxicity_class']
        
        assessment = RegulatoryAssessment(
            biosafety_level=bsl,
            gras_status=gras_info['is_gras'],
            allergenicity_risk=overall_allergy,
            toxicity_risk=overall_tox,
            antibiotic_marker_required=abx_result['removal_required'],
            compliance_score=score,
            flagged_issues=flagged_issues,
            recommendations=recommendations
        )
        
        self.logger.info(f"Compliance Score: {score}/100")
        self.logger.info(f"Flagged Issues: {len(flagged_issues)}")
        
        return assessment

    def run_stage_5_regulatory(self, stage_4_json: Dict[str, Any]) -> Dict[str, Any]:
        """
        Execute Stage 5 Regulatory Assessment logic.
        Input: Stage 4 JSON
        Output: Regulatory compliance JSON
        """
        self.logger.info("=== STAGE 5 REGULATORY ASSESSMENT START ===")
        
        try:
            # Extract info
            org_name = stage_4_json['stage_3_output']['stage_1_output']['organism']['name']
            mol_name = stage_4_json['stage_3_output']['stage_1_output']['target_molecule']['name']
            organism_enum = OrganismType(org_name.replace(" ", "_").replace("-", "_"))
            
            # Count genetic modifications from stage 3
            gene_mods = stage_4_json['stage_3_output'].get('gene_modifications', {})
            num_mods = (len(gene_mods.get('knockouts', [])) + 
                       len(gene_mods.get('overexpressions', [])) + 
                       len(gene_mods.get('heterologous_insertions', [])))
            
            # Generate report
            assessment = self.generate_compliance_report(
                organism=organism_enum,
                molecule_name=mol_name,
                stage_4_json=stage_4_json,
                num_genetic_modifications=num_mods,
                intended_use="food"  # Default assumption
            )
            
            # Format output
            result = {
                "biosafety_assessment": {
                    "biosafety_level": assessment.biosafety_level,
                    "containment_requirements": "Basic microbiological practices" if assessment.biosafety_level == "BSL-1" else "Enhanced containment required"
                },
                "gras_status": {
                    "is_gras": assessment.gras_status,
                    "regulatory_pathway": "GRAS Notification" if assessment.gras_status else "Food Additive Petition / NDA"
                },
                "safety_assessments": {
                    "allergenicity": {
                        "risk_level": assessment.allergenicity_risk,
                        "testing_required": assessment.allergenicity_risk != "LOW"
                    },
                    "toxicity": {
                        "risk_level": assessment.toxicity_risk,
                        "testing_required": assessment.toxicity_risk != "LOW"
                    }
                },
                "genetic_stability": {
                    "antibiotic_marker_removal_required": assessment.antibiotic_marker_required,
                    "recommended_strategy": "Markerless genome editing" if assessment.antibiotic_marker_required else "Standard selection acceptable"
                },
                "compliance_summary": {
                    "overall_score": assessment.compliance_score,
                    "regulatory_feasibility": "HIGH" if assessment.compliance_score >= 70 else ("MEDIUM" if assessment.compliance_score >= 40 else "LOW"),
                    "flagged_issues": assessment.flagged_issues,
                    "critical_recommendations": assessment.recommendations
                },
                "estimated_timeline_months": 12 if assessment.gras_status else 36,
                "estimated_regulatory_cost_usd": 500000 if assessment.gras_status else 2000000
            }
            
            self.logger.info(f"Regulatory assessment complete. Feasibility: {result['compliance_summary']['regulatory_feasibility']}")
            return result
            
        except Exception as e:
            self.logger.error(f"Regulatory assessment failed: {e}", exc_info=True)
            raise PipelineError(f"Stage 5 Regulatory failed: {str(e)}")
