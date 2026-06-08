"""
Fermentation Simulator Module - Bioprocess Simulation with ODE Kinetics
Stage 4: DBTL Loop + Fermentation Simulation

This module implements fermentation process simulation using Monod kinetics
and organism-specific metabolic behaviors.
"""

import logging
import numpy as np
from scipy.integrate import solve_ivp
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Any, Callable
from datetime import datetime
import pandas as pd

from pipeline_config import PipelineConfig
from logger_setup import setup_logger, PipelineLogger
from exceptions import PipelineError, FermentationSimulationError


def get_pipeline_logger(stage: int = 4) -> logging.Logger:
    """Helper to get a logger for a specific stage."""
    # Try to get existing logger or create new one
    from logger_setup import get_existing_logger
    logger_obj = get_existing_logger("pipeline-default")
    if logger_obj:
        return logger_obj.logger
    
    # Fallback: create new logger
    pipeline_logger = setup_logger(
        pipeline_id="pipeline-default",
        stage=stage
    )
    return pipeline_logger.logger


@dataclass
class FermentationState:
    """Current state of the fermentation process."""
    time_hours: float
    biomass_g_per_l: float
    substrate_g_per_l: float
    product_g_per_l: float
    volume_l: float
    dissolved_oxygen_percent: float = 100.0
    ph: float = 7.0
    temperature_c: float = 37.0
    acetate_g_per_l: float = 0.0  # For E. coli
    ethanol_g_per_l: float = 0.0  # For yeast
    sporulation_fraction: float = 0.0  # For B. subtilis
    pha_fraction: float = 0.0  # For P. putida


@dataclass
class FermentationParameters:
    """Kinetic parameters for fermentation model."""
    mu_max: float  # Maximum specific growth rate (1/h)
    Ks: float  # Substrate saturation constant (g/L)
    Yxs: float  # Biomass yield on substrate (g/g)
    Yps: float  # Product yield on substrate (g/g)
    Kd: float  # Death/decay rate (1/h)
    maintenance_coefficient: float  # g substrate / g biomass / h
    oxygen_yield: float  # g biomass / g O2
    critical_do: float  # Critical DO below which growth limited


class FermentationSimulator:
    """
    Simulates fermentation processes using ODE-based kinetic models.
    
    Supports batch, fed-batch, and continuous modes with organism-specific
    metabolic behaviors.
    """
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.logger = get_pipeline_logger(stage=4)
        self.state: Optional[FermentationState] = None
        self.params: Optional[FermentationParameters] = None
        self.time_series: List[FermentationState] = []
        
    def initialize(
        self,
        organism: str,
        initial_biomass: float = 0.5,
        initial_substrate: float = 20.0,
        initial_product: float = 0.0,
        initial_volume: float = 2.0
    ) -> None:
        """
        Initialize fermentation state with organism-specific parameters.
        
        Args:
            organism: Name of the host organism
            initial_biomass: Initial biomass concentration (g/L)
            initial_substrate: Initial substrate concentration (g/L)
            initial_product: Initial product concentration (g/L)
            initial_volume: Initial working volume (L)
        """
        self.logger.info(f"Initializing fermentation for {organism}")
        
        # Get organism-specific parameters
        params = self._get_organism_parameters(organism)
        self.params = params
        
        # Initialize state
        self.state = FermentationState(
            time_hours=0.0,
            biomass_g_per_l=initial_biomass,
            substrate_g_per_l=initial_substrate,
            product_g_per_l=initial_product,
            volume_l=initial_volume,
            dissolved_oxygen_percent=100.0,
            ph=self._get_optimal_ph(organism),
            temperature_c=self._get_optimal_temp(organism)
        )
        
        self.time_series = [self.state]
        self.logger.debug(f"Initial state: X={initial_biomass:.2f}, S={initial_substrate:.2f}")
    
    def _get_organism_parameters(self, organism: str) -> FermentationParameters:
        """Get kinetic parameters for specific organism."""
        organism_lower = organism.lower()
        
        # Default parameters (E. coli-like)
        params = {
            "mu_max": 0.8,
            "Ks": 0.1,
            "Yxs": 0.5,
            "Yps": 0.3,
            "Kd": 0.01,
            "maintenance_coefficient": 0.05,
            "oxygen_yield": 1.4,
            "critical_do": 10.0
        }
        
        # Organism-specific adjustments
        if "coli" in organism_lower:
            params.update({"mu_max": 0.9, "Yxs": 0.45, "Ks": 0.05})
        elif "cerevisiae" in organism_lower or "yeast" in organism_lower:
            params.update({"mu_max": 0.4, "Yxs": 0.5, "Ks": 0.2})
        elif "subtilis" in organism_lower:
            params.update({"mu_max": 0.6, "Yxs": 0.4, "Ks": 0.15})
        elif "glutamicum" in organism_lower:
            params.update({"mu_max": 0.35, "Yxs": 0.45, "Ks": 0.1})
        elif "putida" in organism_lower:
            params.update({"mu_max": 0.5, "Yxs": 0.42, "Ks": 0.08})
        
        return FermentationParameters(**params)
    
    def _get_optimal_ph(self, organism: str) -> float:
        """Get optimal pH for organism."""
        organism_lower = organism.lower()
        if "coli" in organism_lower:
            return 7.0
        elif "cerevisiae" in organism_lower:
            return 5.5
        elif "subtilis" in organism_lower:
            return 6.8
        elif "glutamicum" in organism_lower:
            return 7.2
        elif "putida" in organism_lower:
            return 7.0
        return 7.0
    
    def _get_optimal_temp(self, organism: str) -> float:
        """Get optimal temperature for organism."""
        organism_lower = organism.lower()
        if "coli" in organism_lower:
            return 37.0
        elif "cerevisiae" in organism_lower:
            return 30.0
        elif "subtilis" in organism_lower:
            return 37.0
        elif "glutamicum" in organism_lower:
            return 32.0
        elif "putida" in organism_lower:
            return 30.0
        return 37.0
    
    def monod_kinetics(
        self,
        substrate: float,
        inhibitor: float = 0.0,
        do_limitation: float = 1.0
    ) -> float:
        """
        Calculate specific growth rate using Monod equation with limitations.
        
        μ = μmax × S/(Ks + S) × inhibitor_term × DO_term
        
        Args:
            substrate: Substrate concentration (g/L)
            inhibitor: Inhibitor concentration (e.g., acetate, ethanol)
            do_limitation: DO limitation factor (0-1)
            
        Returns:
            Specific growth rate (1/h)
        """
        if self.params is None:
            raise FermentationSimulationError("Simulator not initialized")
        
        # Basic Monod term
        if substrate + self.params.Ks > 0:
            mu = self.params.mu_max * substrate / (self.params.Ks + substrate)
        else:
            mu = 0.0
        
        # Substrate inhibition (Haldane kinetics)
        if inhibitor > 0:
            Ki = 10.0  # Inhibition constant
            inhibitor_term = 1.0 / (1.0 + inhibitor / Ki)
            mu *= inhibitor_term
        
        # Oxygen limitation
        mu *= do_limitation
        
        return max(0.0, mu)
    
    def acetate_overflow_term(self, acetate: float, growth_rate: float) -> float:
        """
        Calculate acetate production term for E. coli (overflow metabolism).
        
        Acetate is produced when growth rate exceeds critical threshold.
        
        Args:
            acetate: Current acetate concentration
            growth_rate: Current specific growth rate
            
        Returns:
            Acetate production rate (g/L/h)
        """
        critical_mu = 0.4  # Critical growth rate for overflow
        
        if growth_rate > critical_mu:
            # Overflow metabolism active
            excess_mu = growth_rate - critical_mu
            yield_ac = 0.3  # g acetate / g biomass
            production = excess_mu * self.state.biomass_g_per_l * yield_ac
            
            # Acetate consumption at low glucose
            if self.state.substrate_g_per_l < 1.0:
                consumption_rate = 0.1 * acetate
                production -= consumption_rate
            
            return max(0.0, production)
        else:
            # Acetate consumption
            return -0.05 * acetate
    
    def crabtree_effect_term(self, glucose: float, growth_rate: float) -> float:
        """
        Calculate ethanol production term for S. cerevisiae (Crabtree effect).
        
        Ethanol is produced even under aerobic conditions at high glucose.
        
        Args:
            glucose: Glucose concentration
            growth_rate: Specific growth rate
            
        Returns:
            Ethanol production rate (g/L/h)
        """
        critical_glucose = 0.5  # g/L
        
        if glucose > critical_glucose:
            # Crabtree effect active
            excess_glucose = glucose - critical_glucose
            yield_etoh = 0.4  # g ethanol / g glucose
            production = growth_rate * self.state.biomass_g_per_l * yield_etoh
            return production
        else:
            # Ethanol consumption
            return -0.1 * self.state.ethanol_g_per_l
    
    def sporulation_trigger_term(
        self,
        nutrient_level: float,
        biomass: float,
        time_hours: float
    ) -> float:
        """
        Calculate sporulation fraction for B. subtilis.
        
        Sporulation triggered by nutrient depletion and high cell density.
        
        Args:
            nutrient_level: Current nutrient level
            biomass: Current biomass
            time_hours: Time in culture
            
        Returns:
            Rate of change in sporulation fraction
        """
        # Sporulation triggers
        nutrient_starvation = max(0.0, 1.0 - nutrient_level / 5.0)
        density_signal = min(1.0, biomass / 50.0)
        time_signal = min(1.0, max(0.0, (time_hours - 12) / 12))
        
        # Combined trigger
        sporulation_signal = nutrient_starvation * density_signal * time_signal
        
        # Sporulation rate
        max_sporulation_rate = 0.1
        current_fraction = self.state.sporulation_fraction if self.state else 0.0
        
        rate = max_sporulation_rate * sporulation_signal * (1.0 - current_fraction)
        return rate
    
    def amino_acid_secretion_term(
        self,
        biomass: float,
        product: float,
        growth_rate: float
    ) -> float:
        """
        Calculate amino acid secretion for C. glutamicum.
        
        C. glutamicum naturally secretes amino acids due to membrane changes.
        
        Args:
            biomass: Current biomass
            product: Current product concentration
            growth_rate: Specific growth rate
            
        Returns:
            Product secretion rate (g/L/h)
        """
        # Growth-associated production
        growth_assoc = growth_rate * biomass * 0.1
        
        # Non-growth associated production (major for C. glutamicum)
        non_growth_assoc = 0.05 * biomass
        
        # Product inhibition
        if product > 50.0:
            inhibition = 1.0 / (1.0 + (product - 50.0) / 20.0)
        else:
            inhibition = 1.0
        
        return (growth_assoc + non_growth_assoc) * inhibition
    
    def pha_accumulation_term(
        self,
        biomass: float,
        substrate: float,
        nitrogen_limited: bool
    ) -> float:
        """
        Calculate PHA accumulation for P. putida.
        
        PHA accumulates under nitrogen limitation with excess carbon.
        
        Args:
            biomass: Current biomass
            substrate: Substrate concentration
            nitrogen_limited: Whether nitrogen is limiting
            
        Returns:
            PHA accumulation rate (fraction/h)
        """
        if not nitrogen_limited:
            return 0.0
        
        # PHA synthesis under N-limitation
        if substrate > 5.0:
            max_pha_fraction = 0.8  # Max 80% CDW as PHA
            current_pha = self.state.pha_fraction if self.state else 0.0
            
            accumulation_rate = 0.05 * (max_pha_fraction - current_pha)
            return accumulation_rate
        else:
            # PHA degradation when substrate depleted
            return -0.02 * self.state.pha_fraction
    
    def _ode_system(
        self,
        t: float,
        y: np.ndarray,
        feeding_rate: float = 0.0,
        feed_concentration: float = 500.0
    ) -> np.ndarray:
        """
        Define the ODE system for fermentation dynamics.
        
        State variables: [biomass, substrate, product, volume]
        
        Args:
            t: Time (h)
            y: State vector
            feeding_rate: Feed rate (L/h)
            feed_concentration: Feed substrate concentration (g/L)
            
        Returns:
            Derivatives [dX/dt, dS/dt, dP/dt, dV/dt]
        """
        X, S, P, V = y
        
        # Update temporary state for helper functions
        if self.state:
            self.state.substrate_g_per_l = S
            self.state.biomass_g_per_l = X
        
        # Calculate growth rate
        inhibitor = 0.0
        if self.state and hasattr(self.state, 'acetate_g_per_l'):
            inhibitor = self.state.acetate_g_per_l
        
        do_limitation = 1.0 if (self.state and self.state.dissolved_oxygen_percent > 20) else 0.5
        
        mu = self.monod_kinetics(S, inhibitor, do_limitation)
        
        # Biomass balance: dX/dt = μX - KdX - D(X - X_in)
        dX_dt = mu * X - self.params.Kd * X
        
        # Substrate balance: dS/dt = -1/Yxs * μX - m*X + D(S_in - S)
        substrate_consumption = -(1.0 / self.params.Yxs) * mu * X
        maintenance = -self.params.maintenance_coefficient * X
        feed_term_S = feeding_rate / V * (feed_concentration - S) if V > 0 else 0
        dS_dt = substrate_consumption + maintenance + feed_term_S
        
        # Product balance: dP/dt = Yps * μX + beta*X - D*P
        growth_assoc_prod = self.params.Yps * mu * X
        non_growth_assoc = 0.02 * X  # Simple non-growth term
        feed_term_P = -feeding_rate / V * P if V > 0 else 0
        dP_dt = growth_assoc_prod + non_growth_assoc + feed_term_P
        
        # Volume balance: dV/dt = F
        dV_dt = feeding_rate
        
        # Ensure non-negative values
        dX_dt = max(dX_dt, -X / 0.1)  # Can't decrease faster than washout
        dS_dt = max(dS_dt, -S / 0.1)
        dP_dt = max(dP_dt, -P / 0.1)
        
        return np.array([dX_dt, dS_dt, dP_dt, dV_dt])
    
    def run_ode_simulation(
        self,
        duration_hours: float = 72.0,
        mode: str = "fed-batch",
        n_points: int = 100,
        feeding_strategy: Optional[Callable] = None
    ) -> pd.DataFrame:
        """
        Run ODE simulation for fermentation process.
        
        Args:
            duration_hours: Total simulation time
            mode: "batch", "fed-batch", or "continuous"
            n_points: Number of output time points
            feeding_strategy: Function returning feed rate at time t
            
        Returns:
            DataFrame with time-series data
        """
        self.logger.info(
            f"Running ODE simulation: {duration_hours}h, mode={mode}"
        )
        
        if self.state is None or self.params is None:
            raise FermentationSimulationError("Simulator not initialized")
        
        # Set up feeding strategy
        if feeding_strategy is None:
            if mode == "batch":
                feeding_strategy = lambda t: 0.0
            elif mode == "fed-batch":
                feeding_strategy = self._default_fed_batch_strategy
            else:  # continuous
                feeding_strategy = lambda t: 0.1  # Constant dilution
        
        # Initial conditions
        y0 = np.array([
            self.state.biomass_g_per_l,
            self.state.substrate_g_per_l,
            self.state.product_g_per_l,
            self.state.volume_l
        ])
        
        # Time span
        t_span = (0, duration_hours)
        t_eval = np.linspace(0, duration_hours, n_points)
        
        # Solve ODE
        try:
            solution = solve_ivp(
                fun=lambda t, y: self._ode_system(t, y, feeding_strategy(t)),
                t_span=t_span,
                y0=y0,
                t_eval=t_eval,
                method='RK45',
                rtol=1e-6,
                atol=1e-8
            )
            
            ode_converged = solution.success
            
        except Exception as e:
            self.logger.error(f"ODE solver failed: {e}")
            ode_converged = False
            # Create fallback linear profile
            solution = type('obj', (object,), {
                't': t_eval,
                'y': np.zeros((4, len(t_eval)))
            })()
            solution.y[0] = np.linspace(y0[0], y0[0] * 1.5, len(t_eval))
            solution.y[1] = np.linspace(y0[1], 0, len(t_eval))
            solution.y[2] = np.linspace(y0[2], y0[2] + 5, len(t_eval))
            solution.y[3] = np.full(len(t_eval), y0[3])
        
        # Build time-series DataFrame
        df_data = {
            'time_hours': solution.t,
            'biomass_g_per_l': solution.y[0],
            'substrate_g_per_l': solution.y[1],
            'product_g_per_l': solution.y[2],
            'volume_l': solution.y[3]
        }
        
        # Add organism-specific states
        organism = self.config.organism_name
        if "coli" in organism.lower():
            # Simulate acetate profile
            acetate = np.maximum(0, solution.y[0] * 0.1 - solution.y[1] * 0.05)
            df_data['acetate_g_per_l'] = acetate
            df_data['organism_event'] = 'acetate_overflow'
        elif "cerevisiae" in organism.lower():
            # Simulate ethanol profile
            ethanol = np.maximum(0, solution.y[1] * 0.2)
            df_data['ethanol_g_per_l'] = ethanol
            df_data['organism_event'] = 'crabtree_effect'
        elif "subtilis" in organism.lower():
            # Simulate sporulation
            sporulation = np.clip((solution.t - 12) / 24, 0, 1) * np.exp(-solution.y[1] / 5)
            df_data['sporulation_fraction'] = sporulation
            df_data['organism_event'] = 'sporulation'
        elif "glutamicum" in organism.lower():
            df_data['organism_event'] = 'amino_acid_secretion'
        elif "putida" in organism.lower():
            # Simulate PHA accumulation
            pha = np.clip((20 - solution.y[1]) / 20, 0, 0.8)
            df_data['pha_fraction'] = pha
            df_data['organism_event'] = 'pha_accumulation'
        else:
            df_data['organism_event'] = 'none'
        
        df = pd.DataFrame(df_data)
        
        # Store results
        self.time_series_df = df
        self.ode_converged = ode_converged
        
        # Update final state
        final_idx = -1
        self.state = FermentationState(
            time_hours=float(df['time_hours'].iloc[final_idx]),
            biomass_g_per_l=float(df['biomass_g_per_l'].iloc[final_idx]),
            substrate_g_per_l=float(df['substrate_g_per_l'].iloc[final_idx]),
            product_g_per_l=float(df['product_g_per_l'].iloc[final_idx]),
            volume_l=float(df['volume_l'].iloc[final_idx]),
            dissolved_oxygen_percent=30.0,
            ph=self.state.ph,
            temperature_c=self.state.temperature_c
        )
        
        self.logger.info(
            f"Simulation complete: final_titer={self.state.product_g_per_l:.2f} g/L"
        )
        
        return df
    
    def _default_fed_batch_strategy(self, t: float) -> float:
        """
        Default fed-batch feeding strategy (glucose-stat control).
        
        Args:
            t: Time (h)
            
        Returns:
            Feed rate (L/h)
        """
        # Exponential feeding to maintain constant specific growth rate
        target_mu = 0.3  # 1/h
        feed_conc = 500.0  # g/L
        
        if self.state is None or self.params is None:
            return 0.0
        
        # F = (μ * X * V) / (Yxs * (S_feed - S))
        X = self.state.biomass_g_per_l
        V = self.state.volume_l
        S = self.state.substrate_g_per_l
        
        feed_rate = (target_mu * X * V) / (self.params.Yxs * (feed_conc - S))
        
        # Limit feed rate
        feed_rate = min(feed_rate, 0.5)  # Max 0.5 L/h
        feed_rate = max(feed_rate, 0.0)
        
        return feed_rate
    
    def get_optimal_conditions(self) -> Dict[str, float]:
        """
        Determine optimal fermentation conditions based on simulation.
        
        Returns:
            Dictionary of optimal setpoints
        """
        organism = self.config.organism_name
        
        conditions = {
            "temperature_c": self._get_optimal_temp(organism),
            "ph": self._get_optimal_ph(organism),
            "do_percent_saturation": 30.0,
            "glucose_feed_g_per_l_per_h": 1.0,
            "agitation_rpm": 500,
            "aeration_vvm": 1.0
        }
        
        # Organism-specific adjustments
        if "cerevisiae" in organism.lower():
            conditions["do_percent_saturation"] = 20.0  # Lower OK for yeast
            conditions["agitation_rpm"] = 400
        elif "putida" in organism.lower():
            conditions["do_percent_saturation"] = 40.0  # Higher for Pseudomonas
            conditions["agitation_rpm"] = 600
        
        return conditions
    
    def get_final_state_json(self) -> Dict[str, Any]:
        """
        Get final fermentation state as JSON-serializable dict.
        
        Returns:
            Final state dictionary
        """
        if self.state is None:
            raise FermentationSimulationError("No simulation results available")
        
        return {
            "mode": "fed-batch",
            "duration_hours": round(self.state.time_hours, 2),
            "final_titer_g_per_l": round(self.state.product_g_per_l, 4),
            "final_yield_g_per_g": round(
                self.state.product_g_per_l / max(1.0, 20.0 - self.state.substrate_g_per_l),
                4
            ),
            "final_productivity_g_per_l_per_h": round(
                self.state.product_g_per_l / max(1.0, self.state.time_hours),
                4
            ),
            "ode_convergence": getattr(self, 'ode_converged', True),
            "organism_specific_events": self._get_organism_events()
        }
    
    def _get_organism_events(self) -> List[str]:
        """Get list of organism-specific metabolic events."""
        events = []
        organism = self.config.organism_name
        
        if "coli" in organism.lower():
            if self.state and self.state.acetate_g_per_l > 2.0:
                events.append("acetate_overflow_detected")
        elif "cerevisiae" in organism.lower():
            if self.state and self.state.ethanol_g_per_l > 5.0:
                events.append("crabtree_effect_active")
        elif "subtilis" in organism.lower():
            if self.state and self.state.sporulation_fraction > 0.3:
                events.append("sporulation_initiated")
        elif "putida" in organism.lower():
            if self.state and self.state.pha_fraction > 0.5:
                events.append("pha_accumulation_high")
        
        if not events:
            events.append("normal_metabolism")
        
        return events


def run_fermentation_simulation(
    stage_4_partial: Dict[str, Any],
    config: PipelineConfig
) -> Dict[str, Any]:
    """
    Run fermentation simulation and update Stage 4 output.
    
    Args:
        stage_4_partial: Partial Stage 4 output from DBTL
        config: Pipeline configuration
        
    Returns:
        Updated fermentation simulation results
    """
    logger = get_pipeline_logger(stage=4)
    
    try:
        organism = config.organism_name
        dbtl_cycles = stage_4_partial.get("dbtl_cycles", [])
        
        # Get final titer from DBTL
        if dbtl_cycles:
            initial_titer = dbtl_cycles[-1]["best_titer_g_per_l"]
        else:
            initial_titer = 0.5
        
        # Initialize simulator
        simulator = FermentationSimulator(config)
        simulator.initialize(
            organism=organism,
            initial_biomass=0.5,
            initial_substrate=20.0,
            initial_product=initial_titer * 0.1,  # Start with inoculum
            initial_volume=2.0
        )
        
        # Run simulation
        df = simulator.run_ode_simulation(
            duration_hours=72.0,
            mode="fed-batch"
        )
        
        # Get results
        fermentation_results = simulator.get_final_state_json()
        optimal_conditions = simulator.get_optimal_conditions()
        
        logger.info(
            f"Fermentation simulation complete: "
            f"titer={fermentation_results['final_titer_g_per_l']:.2f} g/L"
        )
        
        return {
            "fermentation_simulation": fermentation_results,
            "optimal_fermentation_conditions": optimal_conditions,
            "time_series_summary": {
                "n_points": len(df),
                "max_biomass_g_per_l": round(df['biomass_g_per_l'].max(), 2),
                "min_substrate_g_per_l": round(df['substrate_g_per_l'].min(), 2)
            }
        }
        
    except Exception as e:
        logger.error(f"Fermentation simulation failed: {e}")
        return {
            "fermentation_simulation": {
                "mode": "batch",
                "duration_hours": 0,
                "final_titer_g_per_l": 0,
                "final_yield_g_per_g": 0,
                "final_productivity_g_per_l_per_h": 0,
                "ode_convergence": False,
                "organism_specific_events": [f"ERROR: {str(e)}"]
            },
            "optimal_fermentation_conditions": {},
            "error": str(e)
        }


if __name__ == "__main__":
    # Test fermentation simulator standalone
    print("Testing Fermentation Simulator")
    print("=" * 60)
    
    config = PipelineConfig()
    
    # Test with different organisms
    for organism in ["E. coli K-12 MG1655", "S. cerevisiae S288C", "B. subtilis 168"]:
        print(f"\n--- Testing {organism} ---")
        
        simulator = FermentationSimulator(config)
        simulator.initialize(
            organism=organism,
            initial_biomass=0.5,
            initial_substrate=20.0
        )
        
        df = simulator.run_ode_simulation(duration_hours=48.0, mode="fed-batch")
        results = simulator.get_final_state_json()
        
        print(f"Final titer: {results['final_titer_g_per_l']:.2f} g/L")
        print(f"Productivity: {results['final_productivity_g_per_l_per_h']:.3f} g/L/h")
        print(f"Events: {results['organism_specific_events']}")
