"""
Scale-Up Engine for Industrial Metabolic Pathway Pipeline.
Predicts performance changes when moving from lab-scale to industrial-scale bioreactors.
"""

import logging
import math
import random
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional
from datetime import datetime

# Internal imports
from pipeline_config import PipelineConfig, OrganismType
from exceptions import PipelineError, FermentationSimulationError
from logger_setup import setup_logger, PipelineLogger
import logging

def get_logger(name: str):
    """Get a standard Python logger for the given module name."""
    return logging.getLogger(name)

logger = get_logger(__name__)

@dataclass
class ScaleLevel:
    """Represents a specific scale of operation."""
    name: str
    volume_liters: float
    height_to_diameter_ratio: float = 2.0
    max_agitation_rpm: float = 1200.0
    max_aeration_vvm: float = 2.0
    power_number_np: float = 5.0  # Standard Rushton turbine
    
    # Calculated fields
    kla_per_hour: float = 0.0
    mixing_time_sec: float = 0.0
    power_per_volume_kw_m3: float = 0.0
    tip_speed_m_sec: float = 0.0

@dataclass
class ScaleUpResult:
    """Results of scaling up a specific process."""
    source_scale: str
    target_scale: str
    volume_ratio: float
    kla_change_factor: float
    mixing_time_change_factor: float
    predicted_yield_loss_percent: float
    predicted_titer_g_per_l: float
    critical_parameters: Dict[str, float] = field(default_factory=dict)
    recommendations: List[str] = field(default_factory=list)

class ScaleUpEngine:
    """
    Engine for predicting scale-up effects on fermentation performance.
    Uses dimensional analysis and empirical correlations.
    """
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.logger = get_logger(__name__)
        
        # Empirical constants for kLa correlation: kLa = K * (P/V)^alpha * Vs^beta
        self.KLA_CONSTANTS = {
            'general': {'K': 0.026, 'alpha': 0.4, 'beta': 0.5},
            'E_coli': {'K': 0.030, 'alpha': 0.45, 'beta': 0.55},  # High oxygen demand
            'S_cerevisiae': {'K': 0.025, 'alpha': 0.35, 'beta': 0.45},
            'B_subtilis': {'K': 0.028, 'alpha': 0.4, 'beta': 0.5},
            'C_glutamicum': {'K': 0.027, 'alpha': 0.4, 'beta': 0.5},
            'P_putida': {'K': 0.029, 'alpha': 0.42, 'beta': 0.52}  # High respiration
        }
        
        # Mixing time constants: tm = C * (V)^(1/3) / (N*D^2) ... simplified
        self.MIXING_CONSTANTS = {
            'general': 15.0,
            'viscous': 25.0  # For filamentous or high-cell-density
        }

    def calculate_reactor_geometry(self, volume_liters: float, h_d_ratio: float = 2.0) -> Dict[str, float]:
        """Calculate reactor dimensions based on volume."""
        # V = pi * (D/2)^2 * H = pi * (D/2)^2 * (H_D * D) = pi/4 * H_D * D^3
        # D = (4 * V / (pi * H_D))^(1/3)
        vol_m3 = volume_liters / 1000.0
        diameter_m = math.pow((4.0 * vol_m3) / (math.pi * h_d_ratio), 1.0/3.0)
        height_m = diameter_m * h_d_ratio
        radius_m = diameter_m / 2.0
        
        return {
            'diameter_m': diameter_m,
            'height_m': height_m,
            'radius_m': radius_m,
            'surface_area_m2': math.pi * diameter_m * height_m + 2 * math.pi * radius_m**2
        }

    def predict_kla(self, scale: ScaleLevel, organism: OrganismType) -> float:
        """
        Predict volumetric mass transfer coefficient (kLa).
        Correlation: kLa = K * (P/V)^alpha * Vs^beta
        """
        self.logger.debug(f"Predicting kLa for {organism} at {scale.name}")
        
        consts = self.KLA_CONSTANTS.get(organism.value, self.KLA_CONSTANTS['general'])
        
        # Superficial gas velocity Vs (m/h) approx proportional to VVM * height
        # Simplified: Vs ~ VVM * 60 (very rough approx for correlation)
        vs = scale.max_aeration_vvm * 60.0 
        
        # P/V is provided or calculated. Here we use a typical value for aerobic fermentation
        # Typical P/V range: 1-5 kW/m3
        p_over_v = 2.0 # kW/m3 baseline
        
        kla = consts['K'] * (p_over_v ** consts['alpha']) * (vs ** consts['beta'])
        
        # Scale effect: kLa often decreases slightly at very large scale due to mixing heterogeneity
        if scale.volume_liters > 10000:
            kla *= 0.85
            
        return round(kla, 3)

    def predict_mixing_time(self, scale: ScaleLevel, viscosity_factor: float = 1.0) -> float:
        """
        Predict mixing time in seconds.
        Correlation: tm ~ V^(1/3) / (N * D^2) simplified to function of Volume and Agitation
        """
        geom = self.calculate_reactor_geometry(scale.volume_liters)
        diameter = geom['diameter_m']
        
        # N in rev/sec
        n_rps = scale.max_agitation_rpm / 60.0
        
        # Simplified correlation: tm = C * (V^(1/3)) / (N * D^2) * viscosity_factor
        # Actually tm is often correlated to number of circulations. 
        # Let's use a simpler empirical relation: tm increases with scale
        base_tm = self.MIXING_CONSTANTS['general'] * (scale.volume_liters ** 0.25)
        
        # Agitation reduces mixing time
        agitation_factor = 1000.0 / max(scale.max_agitation_rpm, 100.0)
        
        tm = base_tm * agitation_factor * viscosity_factor
        
        return round(tm, 2)

    def calculate_power_input(self, scale: ScaleLevel) -> float:
        """Calculate Power per Volume (P/V) in kW/m3."""
        geom = self.calculate_reactor_geometry(scale.volume_liters)
        d = geom['diameter_m']
        n_rps = scale.max_agitation_rpm / 60.0
        rho = 1000.0 # kg/m3
        
        # Power number equation: P = Np * rho * N^3 * D^5
        power_watts = scale.power_number_np * rho * (n_rps ** 3) * (d ** 5)
        power_kw = power_watts / 1000.0
        
        vol_m3 = scale.volume_liters / 1000.0
        p_over_v = power_kw / vol_m3
        
        return round(p_over_v, 3)

    def calculate_tip_speed(self, scale: ScaleLevel) -> float:
        """Calculate impeller tip speed in m/s."""
        geom = self.calculate_reactor_geometry(scale.volume_liters)
        d = geom['diameter_m']
        # Impeller diameter is typically 1/3 of tank diameter
        d_impeller = d / 3.0
        n_rps = scale.max_agitation_rpm / 60.0
        
        tip_speed = math.pi * d_impeller * n_rps
        return round(tip_speed, 3)

    def simulate_scale_up(self, 
                          current_titer: float, 
                          organism: OrganismType, 
                          scales: List[ScaleLevel]) -> List[ScaleUpResult]:
        """
        Simulate cascade from lab to pilot to production.
        """
        self.logger.info(f"Simulating scale-up cascade for {organism.value}")
        results = []
        
        current_titer_val = current_titer
        
        for i in range(len(scales) - 1):
            src = scales[i]
            tgt = scales[i+1]
            
            # Update physics for target scale
            tgt.kla_per_hour = self.predict_kla(tgt, organism)
            tgt.mixing_time_sec = self.predict_mixing_time(tgt)
            tgt.power_per_volume_kw_m3 = self.calculate_power_input(tgt)
            tgt.tip_speed_m_sec = self.calculate_tip_speed(tgt)
            
            # Calculate Yield Loss Factors
            loss_factors = []
            recommendations = []
            
            # 1. Oxygen Transfer Limitation
            if tgt.kla_per_hour < src.kla_per_hour * 0.8:
                loss_factors.append(0.05) # 5% loss
                recommendations.append(f"Increase aeration or O2 enrichment at {tgt.name}")
            
            # 2. Mixing Heterogeneity (Gradient formation)
            if tgt.mixing_time_sec > src.mixing_time_sec * 2.0:
                loss_factors.append(0.03) # 3% loss
                recommendations.append(f"Optimize feed point location to avoid substrate gradients in {tgt.name}")
            
            # 3. Shear Stress (Tip speed)
            if tgt.tip_speed_m_sec > 5.0: # High shear threshold
                loss_factors.append(0.02)
                recommendations.append(f"Monitor cell viability due to high tip speed ({tgt.tip_speed_m_sec:.2f} m/s)")
            
            # 4. Organism Specifics
            if organism == OrganismType.ECOLI:
                if tgt.power_per_volume_kw_m3 < 1.0:
                    loss_factors.append(0.04) # Acetate accumulation risk due to poor mixing
                    recommendations.append("Watch for acetate overflow due to local oxygen starvation")
            elif organism == OrganismType.SCEREVISIAE:
                 if tgt.mixing_time_sec > 60:
                    loss_factors.append(0.05) # Crabtree effect exacerbation
                    recommendations.append("Strict fed-batch control required to avoid ethanol formation")
            
            total_loss = sum(loss_factors)
            predicted_titer = current_titer_val * (1.0 - total_loss)
            
            result = ScaleUpResult(
                source_scale=src.name,
                target_scale=tgt.name,
                volume_ratio=tgt.volume_liters / src.volume_liters,
                kla_change_factor=tgt.kla_per_hour / max(src.kla_per_hour, 0.001),
                mixing_time_change_factor=tgt.mixing_time_sec / max(src.mixing_time_sec, 0.001),
                predicted_yield_loss_percent=round(total_loss * 100, 2),
                predicted_titer_g_per_l=round(predicted_titer, 3),
                critical_parameters={
                    "kLa_1h": tgt.kla_per_hour,
                    "mixing_time_s": tgt.mixing_time_sec,
                    "power_per_volume_kw_m3": tgt.power_per_volume_kw_m3,
                    "tip_speed_m_s": tgt.tip_speed_m_sec
                },
                recommendations=recommendations
            )
            
            results.append(result)
            current_titer_val = predicted_titer # Cascade effect
            
            self.logger.debug(f"Scale-up step {src.name}->{tgt.name}: Loss={total_loss*100:.1f}%, New Titer={predicted_titer:.2f}")
            
        return results

    def run_stage_5_scaleup(self, stage_4_json: Dict[str, Any]) -> Dict[str, Any]:
        """
        Execute Stage 5 Scale-Up logic.
        Input: Stage 4 JSON (Fermentation Simulation results)
        Output: Scale-up assessment JSON
        """
        self.logger.info("=== STAGE 5 SCALE-UP START ===")
        
        try:
            # Extract data from Stage 4
            org_name = stage_4_json['stage_3_output']['stage_1_output']['organism']['name']
            organism_enum = OrganismType(org_name.replace(" ", "_").replace("-", "_"))
            
            lab_titer = stage_4_json['fermentation_simulation']['final_titer_g_per_l']
            
            # Define standard scale-up cascade
            scales = [
                ScaleLevel(name="Lab_Scale", volume_liters=2.0, max_agitation_rpm=1000, max_aeration_vvm=1.5),
                ScaleLevel(name="Pilot_Plant", volume_liters=200.0, max_agitation_rpm=600, max_aeration_vvm=1.2),
                ScaleLevel(name="Production", volume_liters=20000.0, max_agitation_rpm=300, max_aeration_vvm=0.8)
            ]
            
            # Initialize physics for lab scale
            scales[0].kla_per_hour = self.predict_kla(scales[0], organism_enum)
            scales[0].mixing_time_sec = self.predict_mixing_time(scales[0])
            
            results = self.simulate_scale_up(lab_titer, organism_enum, scales)
            
            # Format output
            output_data = {
                "scale_up_cascade": [
                    {
                        "transition": f"{r.source_scale} -> {r.target_scale}",
                        "volume_ratio": r.volume_ratio,
                        "predicted_titer_g_per_l": r.predicted_titer_g_per_l,
                        "yield_loss_percent": r.predicted_yield_loss_percent,
                        "critical_parameters": r.critical_parameters,
                        "recommendations": r.recommendations
                    } for r in results
                ],
                "final_production_titer_g_per_l": results[-1].predicted_titer_g_per_l if results else lab_titer,
                "scale_up_risk": "HIGH" if any(r.predicted_yield_loss_percent > 10 for r in results) else "LOW"
            }
            
            self.logger.info(f"Scale-up simulation complete. Final production titer: {output_data['final_production_titer_g_per_l']} g/L")
            return output_data
            
        except Exception as e:
            self.logger.error(f"Scale-up simulation failed: {e}", exc_info=True)
            raise FermentationSimulationError(f"Stage 5 Scale-Up failed: {str(e)}")
