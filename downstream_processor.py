"""
Downstream Processing Module for Industrial Metabolic Pathway Pipeline.
Simulates purification, recovery, and cost estimation for bioproducts.
"""

import logging
import math
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple
from enum import Enum

# Internal imports
from pipeline_config import PipelineConfig, OrganismType, MoleculeType
from exceptions import PipelineError
import logging

def get_logger(name: str):
    """Get a standard Python logger for the given module name."""
    return logging.getLogger(name)

logger = get_logger(__name__)

class ProductLocation(Enum):
    INTRACELLULAR = "intracellular"
    SECRETED = "secreted"
    PERIPLASMIC = "periplasmic"
    EXTRACELLULAR = "extracellular"

class ChromatographyMode(Enum):
    ION_EXCHANGE = "ion_exchange"
    HYDROPHOBIC_INTERACTION = "hydrophobic_interaction"
    AFFINITY = "affinity"
    SIZE_EXCLUSION = "size_exclusion"
    REVERSE_PHASE = "reverse_phase"

@dataclass
class PurificationStep:
    """Represents a single unit operation in downstream processing."""
    step_name: str
    step_type: str  # harvest, lysis, clarification, chromatography, crystallization, etc.
    recovery_percent: float
    purity_increase_factor: float
    cost_per_kg_product: float
    description: str = ""

@dataclass
class DownstreamResult:
    """Complete downstream processing results."""
    product_location: ProductLocation
    purification_train: List[PurificationStep]
    overall_recovery_percent: float
    final_purity_percent: float
    total_cost_per_kg_usd: float
    critical_steps: List[str] = field(default_factory=list)

class DownstreamProcessor:
    """
    Simulates downstream processing for various products and organisms.
    """
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.logger = get_logger(__name__)
        
        # Product location classification based on molecule type
        self.PRODUCT_LOCATIONS = {
            'lycopene': ProductLocation.INTRACELLULAR,
            'vanillin': ProductLocation.SECRETED,
            'l_lysine': ProductLocation.SECRETED,
            'l_glutamate': ProductLocation.SECRETED,
            'artemisinic_acid': ProductLocation.SECRETED,
            'riboflavin': ProductLocation.SECRETED,
            'pha': ProductLocation.INTRACELLULAR,
            'hyaluronic_acid': ProductLocation.EXTRACELLULAR
        }
        
        # Default locations if not specified
        self.DEFAULT_LOCATION = ProductLocation.INTRACELLULAR
        
        # Chromatography selection rules based on product properties
        # (Simplified: based on molecular weight and polarity)
        self.CHROMATOGRAPHY_RULES = {
            'small_molecule_polar': [ChromatographyMode.ION_EXCHANGE, ChromatographyMode.REVERSE_PHASE],
            'small_molecule_nonpolar': [ChromatographyMode.HYDROPHOBIC_INTERACTION, ChromatographyMode.REVERSE_PHASE],
            'large_molecule_protein': [ChromatographyMode.AFFINITY, ChromatographyMode.ION_EXCHANGE, ChromatographyMode.SIZE_EXCLUSION],
            'polymer': [ChromatographyMode.SIZE_EXCLUSION, ChromatographyMode.ION_EXCHANGE]
        }

    def classify_product_location(self, molecule_name: str, organism: OrganismType) -> ProductLocation:
        """Determine where the product accumulates."""
        self.logger.debug(f"Classifying product location for {molecule_name}")
        
        # Check explicit mapping
        mol_key = molecule_name.lower().replace(" ", "_").replace("-", "_")
        if mol_key in self.PRODUCT_LOCATIONS:
            return self.PRODUCT_LOCATIONS[mol_key]
        
        # Heuristics based on molecule name
        if 'acid' in mol_key and 'amino' not in mol_key:
            return ProductLocation.SECRETED
        elif 'amino' in mol_key or 'lysine' in mol_key or 'glutamate' in mol_key:
            return ProductLocation.SECRETED
        elif 'pha' in mol_key or 'polyhydroxy' in mol_key:
            return ProductLocation.INTRACELLULAR
            
        return self.DEFAULT_LOCATION

    def select_chromatography_mode(self, molecule_name: str, molecular_weight: float = 200.0) -> List[ChromatographyMode]:
        """Select appropriate chromatography methods based on product properties."""
        self.logger.debug(f"Selecting chromatography for MW={molecular_weight}")
        
        # Simple MW-based classification
        if molecular_weight < 500:
            if 'polar' in molecule_name.lower() or 'acid' in molecule_name.lower() or 'amino' in molecule_name.lower():
                return self.CHROMATOGRAPHY_RULES['small_molecule_polar']
            else:
                return self.CHROMATOGRAPHY_RULES['small_molecule_nonpolar']
        elif molecular_weight < 50000:
            return self.CHROMATOGRAPHY_RULES['large_molecule_protein']
        else:
            return self.CHROMATOGRAPHY_RULES['polymer']

    def organism_specific_harvest_strategy(self, organism: OrganismType) -> Tuple[str, float, str]:
        """
        Returns: (harvest_method, recovery_efficiency, description)
        """
        strategies = {
            OrganismType.ECOLI: (
                "Centrifugation + High-pressure homogenization",
                0.92,
                "E. coli requires mechanical lysis due to robust cell wall. Use 2 passes at 15,000 psi."
            ),
            OrganismType.SCEREVISIAE: (
                "Centrifugation + Enzymatic lysis (Zymolyase)",
                0.88,
                "Yeast cell wall requires enzymatic digestion. Costlier but gentler than bead beating."
            ),
            OrganismType.BSUBTILIS: (
                "Centrifugation + Bead milling",
                0.90,
                "Bacillus forms spores; ensure vegetative cells only. Bead milling effective for Gram-positive."
            ),
            OrganismType.CGLUTAMICUM: (
                "Centrifugation (product secreted)",
                0.95,
                "C. glutamicum typically secretes amino acids. Simple centrifugation removes biomass."
            ),
            OrganismType.PPUTIDA: (
                "Centrifugation + Solvent extraction (for PHA)",
                0.85,
                "Pseudomonas often accumulates PHA granules. Solvent extraction with chloroform or eco-friendly alternatives."
            )
        }
        
        return strategies.get(organism, ("Centrifugation", 0.90, "Standard harvest protocol"))

    def simulate_purification_train(self, 
                                    product_location: ProductLocation,
                                    chromatography_modes: List[ChromatographyMode],
                                    target_purity: float = 98.0) -> List[PurificationStep]:
        """Design a complete purification train."""
        self.logger.info(f"Designing purification train for {product_location.value}")
        
        steps = []
        current_purity = 10.0  # Start with crude broth (10% purity estimate)
        cumulative_recovery = 1.0
        
        # Step 1: Harvest/Clarification
        if product_location in [ProductLocation.INTRACELLULAR, ProductLocation.PERIPLASMIC]:
            steps.append(PurificationStep(
                step_name="Cell Harvest",
                step_type="centrifugation",
                recovery_percent=95.0,
                purity_increase_factor=1.2,
                cost_per_kg_product=5.0,
                description="Discard supernatant, keep cell pellet"
            ))
            current_purity *= 1.2
            cumulative_recovery *= 0.95
            
            # Lysis step
            steps.append(PurificationStep(
                step_name="Cell Lysis",
                step_type="lysis",
                recovery_percent=90.0,
                purity_increase_factor=1.1,
                cost_per_kg_product=15.0,
                description="Mechanical or enzymatic disruption"
            ))
            current_purity *= 1.1
            cumulative_recovery *= 0.90
            
            # Debris removal
            steps.append(PurificationStep(
                step_name="Debris Clarification",
                step_type="filtration",
                recovery_percent=92.0,
                purity_increase_factor=1.3,
                cost_per_kg_product=8.0,
                description="Remove cell debris by depth filtration"
            ))
            current_purity *= 1.3
            cumulative_recovery *= 0.92
            
        else:  # Secreted/Extracellular
            steps.append(PurificationStep(
                step_name="Biomass Removal",
                step_type="centrifugation",
                recovery_percent=97.0,
                purity_increase_factor=1.5,
                cost_per_kg_product=4.0,
                description="Remove cells, keep supernatant containing product"
            ))
            current_purity *= 1.5
            cumulative_recovery *= 0.97
        
        # Step 2: Capture (First Chromatography)
        if chromatography_modes:
            cap_mode = chromatography_modes[0]
            cap_recovery = 85.0
            cap_purity_factor = 4.0
            
            steps.append(PurificationStep(
                step_name=f"Capture ({cap_mode.value.replace('_', ' ').title()})",
                step_type="chromatography",
                recovery_percent=cap_recovery,
                purity_increase_factor=cap_purity_factor,
                cost_per_kg_product=50.0,
                description=f"Primary capture using {cap_mode.value}"
            ))
            current_purity *= cap_purity_factor
            cumulative_recovery *= (cap_recovery / 100.0)
            
            # Step 3: Polishing (Second Chromatography if needed)
            if current_purity < target_purity and len(chromatography_modes) > 1:
                pol_mode = chromatography_modes[1]
                pol_recovery = 90.0
                pol_purity_factor = 2.5
                
                steps.append(PurificationStep(
                    step_name=f"Polishing ({pol_mode.value.replace('_', ' ').title()})",
                    step_type="chromatography",
                    recovery_percent=pol_recovery,
                    purity_increase_factor=pol_purity_factor,
                    cost_per_kg_product=60.0,
                    description=f"Polishing step using {pol_mode.value}"
                ))
                current_purity *= pol_purity_factor
                cumulative_recovery *= (pol_recovery / 100.0)
        
        # Step 4: Final concentration/crystallization
        if current_purity < target_purity:
            steps.append(PurificationStep(
                step_name="Crystallization/Precipitation",
                step_type="crystallization",
                recovery_percent=88.0,
                purity_increase_factor=1.8,
                cost_per_kg_product=20.0,
                description="Final purification and solid form generation"
            ))
            current_purity *= 1.8
            cumulative_recovery *= 0.88
        
        # Update actual recovery and purity in steps
        running_recovery = 1.0
        running_purity = 10.0
        for step in steps:
            running_recovery *= (step.recovery_percent / 100.0)
            running_purity *= step.purity_increase_factor
            step.recovery_percent = round(running_recovery * 100, 2)  # Cumulative
            # Purity increase factor stays as-is for reference
            
        return steps

    def calculate_cost_per_kg(self, steps: List[PurificationStep], overall_recovery: float) -> float:
        """Calculate total downstream cost per kg of pure product."""
        # Base cost from steps
        base_cost = sum(step.cost_per_kg_product for step in steps)
        
        # Adjust for recovery losses (you need to process more to get 1 kg)
        if overall_recovery > 0:
            true_cost = base_cost / (overall_recovery / 100.0)
        else:
            true_cost = float('inf')
            
        return round(true_cost, 2)

    def run_stage_5_downstream(self, stage_4_json: Dict[str, Any]) -> Dict[str, Any]:
        """
        Execute Stage 5 Downstream Processing logic.
        Input: Stage 4 JSON
        Output: Downstream assessment JSON
        """
        self.logger.info("=== STAGE 5 DOWNSTREAM START ===")
        
        try:
            # Extract info
            org_name = stage_4_json['stage_3_output']['stage_1_output']['organism']['name']
            mol_name = stage_4_json['stage_3_output']['stage_1_output']['target_molecule']['name']
            organism_enum = OrganismType(org_name.replace(" ", "_").replace("-", "_"))
            
            # Classify product
            location = self.classify_product_location(mol_name, organism_enum)
            
            # Get harvest strategy
            harvest_method, harvest_recovery, harvest_desc = self.organism_specific_harvest_strategy(organism_enum)
            
            # Select chromatography (using dummy MW for now)
            mw_estimates = {
                'lycopene': 536.9, 'vanillin': 152.1, 'l_lysine': 146.2,
                'l_glutamate': 147.1, 'riboflavin': 376.4, 'pha': 100000.0
            }
            mw = mw_estimates.get(mol_name.lower(), 200.0)
            chrom_modes = self.select_chromatography_mode(mol_name, mw)
            
            # Design purification train
            steps = self.simulate_purification_train(location, chrom_modes)
            
            # Calculate overall metrics
            if steps:
                overall_recovery = steps[-1].recovery_percent
                final_purity = min(99.5, 10.0 * math.prod(s.purity_increase_factor for s in steps))
            else:
                overall_recovery = 100.0
                final_purity = 10.0
                
            total_cost = self.calculate_cost_per_kg(steps, overall_recovery)
            
            # Identify critical steps (lowest recovery)
            critical_steps = []
            if len(steps) >= 2:
                recoveries = [s.recovery_percent / (steps[i-1].recovery_percent if i > 0 else 100) 
                             for i, s in enumerate(steps)]
                min_rec_idx = recoveries.index(min(recoveries))
                critical_steps.append(steps[min_rec_idx].step_name)
            
            result = {
                "product_location": location.value,
                "harvest_strategy": {
                    "method": harvest_method,
                    "recovery_percent": harvest_recovery * 100,
                    "description": harvest_desc
                },
                "purification_train": [
                    {
                        "step_name": s.step_name,
                        "step_type": s.step_type,
                        "cumulative_recovery_percent": s.recovery_percent,
                        "purity_increase_factor": s.purity_increase_factor,
                        "cost_contribution_usd_per_kg": s.cost_per_kg_product,
                        "description": s.description
                    } for s in steps
                ],
                "overall_recovery_percent": round(overall_recovery, 2),
                "final_purity_percent": round(final_purity, 2),
                "estimated_cost_per_kg_usd": total_cost,
                "critical_steps": critical_steps,
                "economic_feasibility": "HIGH" if total_cost < 100 else ("MEDIUM" if total_cost < 500 else "LOW")
            }
            
            self.logger.info(f"Downstream design complete. Recovery: {overall_recovery:.1f}%, Cost: ${total_cost}/kg")
            return result
            
        except Exception as e:
            self.logger.error(f"Downstream processing simulation failed: {e}", exc_info=True)
            raise PipelineError(f"Stage 5 Downstream failed: {str(e)}")
