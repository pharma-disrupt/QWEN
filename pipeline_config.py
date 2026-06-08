"""
Pipeline Configuration Module
Defines all configurable parameters for the metabolic pathway pipeline.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional
from enum import Enum


class OrganismType(Enum):
    """Supported industrial microorganisms."""
    ECOLI = "ecoli"
    SCEREVISIAE = "scerevisiae"
    BSUBTILIS = "bsubtilis"
    CGLUTAMICUM = "cglutamicum"
    PPUTIDA = "pputida"


class MoleculeType(Enum):
    """Target product categories."""
    LYCOPENE = "lycopene"
    VANILLIN = "vanillin"
    ARTEMISINIC_ACID = "artemisinic_acid"
    L_LYSINE = "l_lysine"
    L_GLUTAMATE = "l_glutamate"
    L_THREONINE = "l_threonine"
    PHA = "pha"
    HYALURONIC_ACID = "hyaluronic_acid"
    RIBOFLAVIN = "riboflavin"


@dataclass
class OrganismConfig:
    """Configuration for a specific organism strain."""
    name: str
    strain: str
    gem_model_id: str
    doubling_time_min: float
    optimal_ph: float
    optimal_temp_c: float
    gram_stain: str
    codon_table_id: int
    gc_content_percent: float
    essential_genes: List[str]
    available_promoters: List[str]
    mu_max_per_hour: float
    ks_glucose_mm: float
    maintenance_coefficient: float
    
    @classmethod
    def get_ecoli_config(cls) -> 'OrganismConfig':
        """E. coli K-12 MG1655 configuration."""
        return cls(
            name="Escherichia coli",
            strain="K-12 MG1655",
            gem_model_id="iML1515",
            doubling_time_min=20.0,
            optimal_ph=7.0,
            optimal_temp_c=37.0,
            gram_stain="negative",
            codon_table_id=11,
            gc_content_percent=50.8,
            essential_genes=[
                "dnaA", "dnaN", "gyrA", "gyrB", "rpoA", "rpoB", "rpoC",
                "ftsZ", "murA", "murB", "fabH", "fabD", "accA", "accD",
                "rpsA", "rpsB", "rplA", "rplB", "tufA", "tufB", "tsf",
                "gapA", "pgk", "tpiA", "pykF", "pfkA", "zwf", "gnd",
                "acnA", "acnB", "icd", "sucA", "sucB", "sucC", "sucD",
                "mdh", "fumA", "fumB", "fumC", "nuoA", "nuoB", "cyoA"
            ],
            available_promoters=[
                "J23100", "J23101", "J23102", "J23103", "J23104",
                "J23105", "J23106", "J23107", "J23108", "J23109",
                "Ptac", "Ptrc", "Plac", "PBAD", "Prha", "Ptet"
            ],
            mu_max_per_hour=3.0,
            ks_glucose_mm=0.1,
            maintenance_coefficient=0.05
        )
    
    @classmethod
    def get_ecoli_bl21_config(cls) -> 'OrganismConfig':
        """E. coli BL21(DE3) configuration."""
        config = cls.get_ecoli_config()
        config.strain = "BL21(DE3)"
        config.gem_model_id = "iJO1366"
        config.available_promoters.extend(["PT7", "PT7lacO"])
        return config
    
    @classmethod
    def get_scerevisiae_config(cls) -> 'OrganismConfig':
        """S. cerevisiae S288C configuration."""
        return cls(
            name="Saccharomyces cerevisiae",
            strain="S288C",
            gem_model_id="yeast8",
            doubling_time_min=90.0,
            optimal_ph=5.5,
            optimal_temp_c=30.0,
            gram_stain="positive",
            codon_table_id=12,
            gc_content_percent=38.3,
            essential_genes=[
                "CDC28", "CLN1", "CLN2", "CLB1", "CLB2", "TUB1", "TUB2",
                "ACT1", "MYO1", "RPB1", "RPB2", "RPB3", "RPC1", "RPC2",
                "NOP1", "NSR1", "NUP1", "NUP2", "GSP1", "RAN1",
                "ERG1", "ERG2", "ERG3", "ERG4", "ERG5", "ERG6", "ERG7",
                "ADE1", "ADE2", "ADE3", "HIS1", "HIS3", "HIS4", "HIS5",
                "LEU1", "LEU2", "TRP1", "TRP2", "TRP3", "TRP4", "TRP5"
            ],
            available_promoters=[
                "PGAL1", "PGAL10", "PGPD", "PTEF1", "PACT1", "PADH1",
                "PCYC1", "PPGK1", "PHXT1", "PHXT7", "PTDH3", "PCUP1"
            ],
            mu_max_per_hour=0.47,
            ks_glucose_mm=1.0,
            maintenance_coefficient=0.02
        )
    
    @classmethod
    def get_bsubtilis_config(cls) -> 'OrganismConfig':
        """B. subtilis 168 configuration."""
        return cls(
            name="Bacillus subtilis",
            strain="168",
            gem_model_id="iBSU1103",
            doubling_time_min=25.0,
            optimal_ph=7.0,
            optimal_temp_c=37.0,
            gram_stain="positive",
            codon_table_id=11,
            gc_content_percent=43.5,
            essential_genes=[
                "dnaA", "dnaB", "dnaD", "dnaI", "dnaC", "gyrA", "gyrB",
                "rpoA", "rpoB", "rpoC", "rpoD", "sigA", "ftsZ", "ftsA",
                "murAA", "murAB", "murB", "murC", "murD", "murE", "murF",
                "accC", "accD", "accA", "accB", "fabH", "fabD", "fabF",
                "rpsA", "rpsB", "rpsC", "rplA", "rplB", "rplC", "tuf"
            ],
            available_promoters=[
                "Pveg", "PaprE", "PamyE", "PspoIIG", "Phag", "Pylb",
                "P43", "Pgrac", "Pxyl", "PglnA", "Pipt", "PsrfA"
            ],
            mu_max_per_hour=2.4,
            ks_glucose_mm=0.5,
            maintenance_coefficient=0.04
        )
    
    @classmethod
    def get_cglutamicum_config(cls) -> 'OrganismConfig':
        """C. glutamicum ATCC 13032 configuration."""
        return cls(
            name="Corynebacterium glutamicum",
            strain="ATCC 13032",
            gem_model_id="iCW773",
            doubling_time_min=60.0,
            optimal_ph=7.0,
            optimal_temp_c=30.0,
            gram_stain="positive",
            codon_table_id=11,
            gc_content_percent=53.8,
            essential_genes=[
                "dnaA", "dnaN", "gyrA", "gyrB", "rpoA", "rpoB", "rpoC",
                "ftsZ", "ftsQ", "ftsA", "murA", "murB", "murC", "murD",
                "accD", "accA", "accB", "fabH", "fabD", "fabF", "fabG",
                "rpsA", "rpsB", "rplA", "rplB", "tuf", "tsf", "fusA"
            ],
            available_promoters=[
                "Ptac", "Psod", "Pgap", "Ptuf", "Pffs", "PdapA",
                "Pcg2927", "Pcg1766", "Pcg0933", "Pinducible_tet"
            ],
            mu_max_per_hour=0.69,
            ks_glucose_mm=0.2,
            maintenance_coefficient=0.03
        )
    
    @classmethod
    def get_pputida_config(cls) -> 'OrganismConfig':
        """P. putida KT2440 configuration."""
        return cls(
            name="Pseudomonas putida",
            strain="KT2440",
            gem_model_id="iJN1462",
            doubling_time_min=40.0,
            optimal_ph=7.0,
            optimal_temp_c=30.0,
            gram_stain="negative",
            codon_table_id=11,
            gc_content_percent=61.6,
            essential_genes=[
                "dnaA", "dnaN", "gyrA", "gyrB", "rpoA", "rpoB", "rpoC",
                "rpoD", "sigN", "ftsZ", "ftsA", "murA", "murB", "murC",
                "accA", "accB", "accC", "accD", "fabH", "fabD", "fabF",
                "rpsA", "rpsB", "rpsC", "rplA", "rplB", "tuf", "tsf"
            ],
            available_promoters=[
                "Ptac", "Pbad", "Prha", "Pm", "Ps", "PalkB",
                "Pcat", "Pben", "PxylS", "Ptet", "Plac", "Ptrc"
            ],
            mu_max_per_hour=0.9,
            ks_glucose_mm=0.15,
            maintenance_coefficient=0.06
        )


@dataclass
class MoleculeConfig:
    """Configuration for target molecule production."""
    name: str
    smiles: str
    chebi_id: str
    target_titer_g_per_l: float
    target_yield_mol_per_mol: float
    molecular_weight: float
    category: str
    pathway_type: str
    toxicity_level: str  # LOW, MEDIUM, HIGH
    secretion_type: str  # intracellular, secreted, periplasmic
    
    @classmethod
    def get_lycopene_config(cls) -> 'MoleculeConfig':
        """Lycopene production configuration."""
        return cls(
            name="lycopene",
            smiles="CC1=CCC2(C1)CC=C3C(=CC=C4C(=C3)C(=CC=C5C(=C4)CC=C6C(=C5)CC(=CC=C6)C)C)C2",
            chebi_id="CHEBI:15956",
            target_titer_g_per_l=5.0,
            target_yield_mol_per_mol=0.3,
            molecular_weight=536.87,
            category="carotenoid",
            pathway_type="terpenoid",
            toxicity_level="LOW",
            secretion_type="intracellular"
        )
    
    @classmethod
    def get_vanillin_config(cls) -> 'MoleculeConfig':
        """Vanillin production configuration."""
        return cls(
            name="vanillin",
            smiles="COc1ccc(C=O)cc1O",
            chebi_id="CHEBI:16938",
            target_titer_g_per_l=3.0,
            target_yield_mol_per_mol=0.4,
            molecular_weight=152.15,
            category="phenylpropanoid",
            pathway_type="aromatic",
            toxicity_level="MEDIUM",
            secretion_type="secreted"
        )
    
    @classmethod
    def get_artemisinic_acid_config(cls) -> 'MoleculeConfig':
        """Artemisinic acid production configuration."""
        return cls(
            name="artemisinic acid",
            smiles="CC(=C)C1CC2C(C1(C)C)CC=C(C2)C(=O)O",
            chebi_id="CHEBI:2256",
            target_titer_g_per_l=25.0,
            target_yield_mol_per_mol=0.25,
            molecular_weight=234.33,
            category="sesquiterpene",
            pathway_type="terpenoid",
            toxicity_level="MEDIUM",
            secretion_type="secreted"
        )
    
    @classmethod
    def get_l_lysine_config(cls) -> 'MoleculeConfig':
        """L-lysine production configuration."""
        return cls(
            name="L-lysine",
            smiles="NCCCCC(N)C(=O)O",
            chebi_id="CHEBI:25654",
            target_titer_g_per_l=100.0,
            target_yield_mol_per_mol=0.6,
            molecular_weight=146.19,
            category="amino_acid",
            pathway_type="aspartate_family",
            toxicity_level="LOW",
            secretion_type="secreted"
        )
    
    @classmethod
    def get_l_glutamate_config(cls) -> 'MoleculeConfig':
        """L-glutamate production configuration."""
        return cls(
            name="L-glutamate",
            smiles="NC(CCC(=O)O)C(=O)O",
            chebi_id="CHEBI:29985",
            target_titer_g_per_l=80.0,
            target_yield_mol_per_mol=0.7,
            molecular_weight=147.13,
            category="amino_acid",
            pathway_type="tca_cycle_derived",
            toxicity_level="LOW",
            secretion_type="secreted"
        )
    
    @classmethod
    def get_l_threonine_config(cls) -> 'MoleculeConfig':
        """L-threonine production configuration."""
        return cls(
            name="L-threonine",
            smiles="CC(O)C(N)C(=O)O",
            chebi_id="CHEBI:57926",
            target_titer_g_per_l=60.0,
            target_yield_mol_per_mol=0.5,
            molecular_weight=119.12,
            category="amino_acid",
            pathway_type="aspartate_family",
            toxicity_level="LOW",
            secretion_type="secreted"
        )
    
    @classmethod
    def get_pha_config(cls) -> 'MoleculeConfig':
        """PHA (polyhydroxyalkanoate) production configuration."""
        return cls(
            name="PHA",
            smiles="CC(CC(=O)O)O",  # Simplified monomer unit
            chebi_id="CHEBI:47783",
            target_titer_g_per_l=50.0,
            target_yield_mol_per_mol=0.4,
            molecular_weight=86.09,  # Per monomer
            category="biopolymer",
            pathway_type="fatty_acid_derived",
            toxicity_level="LOW",
            secretion_type="intracellular"
        )
    
    @classmethod
    def get_hyaluronic_acid_config(cls) -> 'MoleculeConfig':
        """Hyaluronic acid production configuration."""
        return cls(
            name="hyaluronic acid",
            smiles="CC(=O)NC1C(CC(C(C1O)O)O)OC(C(C(CO)O)O)C(=O)O",  # Disaccharide unit
            chebi_id="CHEBI:28440",
            target_titer_g_per_l=10.0,
            target_yield_mol_per_mol=0.3,
            molecular_weight=401.38,  # Per disaccharide
            category="glycosaminoglycan",
            pathway_type="sugar_nucleotide",
            toxicity_level="LOW",
            secretion_type="secreted"
        )
    
    @classmethod
    def get_riboflavin_config(cls) -> 'MoleculeConfig':
        """Riboflavin (Vitamin B2) production configuration."""
        return cls(
            name="riboflavin",
            smiles="CC1C(C(C(C(O1)O)O)O)OC2C(=O)Nc3nc(c4c3nc(=O)c(=CN4C)C)O",
            chebi_id="CHEBI:16126",
            target_titer_g_per_l=15.0,
            target_yield_mol_per_mol=0.35,
            molecular_weight=376.36,
            category="vitamin",
            pathway_type="ribulose_5phosphate",
            toxicity_level="LOW",
            secretion_type="secreted"
        )


@dataclass
class LoggingConfig:
    """Logging configuration for the pipeline."""
    log_dir: str = "logs"
    log_level: str = "DEBUG"
    console_output: bool = True
    file_output: bool = True
    log_format: str = "%(asctime)s | %(levelname)-8s | STAGE:%(stage)s | %(funcName)s | %(message)s"
    date_format: str = "%Y-%m-%d %H:%M:%S"
    max_log_size_mb: int = 100
    backup_count: int = 5


@dataclass
class FBAConfig:
    """Flux Balance Analysis configuration."""
    objective_reaction: str = "BIOMASS"
    substrate_uptake_limit: float = 10.0  # mmol/gDW/h
    oxygen_uptake_limit: float = 15.0  # mmol/gDW/h
    atp_maintenance: float = 8.39  # mmol/gDW/h
    ngam: float = 4.18  # Non-growth associated maintenance
    fva_tolerance: float = 1e-6
    max_iterations: int = 1000


@dataclass
class FermentationConfig:
    """Fermentation simulation configuration."""
    mode: str = "fed_batch"  # batch, fed_batch, continuous
    initial_volume_l: float = 1.0
    max_volume_l: float = 2.0
    initial_glucose_g_per_l: float = 20.0
    feed_glucose_g_per_l: float = 500.0
    duration_hours: float = 48.0
    sampling_interval_min: float = 30.0
    temperature_setpoint_c: float = 37.0
    ph_setpoint: float = 7.0
    do_setpoint_percent: float = 30.0
    agitation_rpm: float = 500.0
    aeration_vvm: float = 1.0


@dataclass
class DBTLConfig:
    """Design-Build-Test-Learn cycle configuration."""
    num_cycles: int = 3
    constructs_per_cycle: int = 48
    screening_noise_std: float = 0.1
    bayesian_optimization_candidates: int = 12
    improvement_threshold: float = 1.2  # Minimum fold improvement to continue


@dataclass
class ScaleUpConfig:
    """Scale-up configuration."""
    scales: List[float] = field(default_factory=lambda: [2.0, 200.0, 20000.0])  # Liters
    kla_target_per_hour: float = 100.0
    mixing_time_target_min: float = 2.0
    power_per_volume_kw_per_m3: float = 1.0


@dataclass
class DownstreamConfig:
    """Downstream processing configuration."""
    harvest_method: str = "centrifugation"  # centrifugation, filtration, flocculation
    lysis_method: str = "homogenization"  # homogenization, sonication, enzymatic
    chromatography_modes: List[str] = field(default_factory=lambda: ["ion_exchange", "hydrophobic_interaction", "size_exclusion"])
    target_purity_percent: float = 95.0
    minimum_recovery_percent: float = 70.0


@dataclass
class RegulatoryConfig:
    """Regulatory compliance configuration."""
    biosafety_level_target: int = 1
    require_gras_status: bool = True
    require_antibiotic_free: bool = True
    allergenicity_threshold: float = 0.5
    toxicity_threshold: float = 0.5


@dataclass
class PipelineConfig:
    """Master configuration for the entire pipeline."""
    # Core settings
    organism_type: OrganismType = OrganismType.ECOLI
    molecule_type: MoleculeType = MoleculeType.LYCOPENE
    
    # Sub-configurations
    organism_config: Optional[OrganismConfig] = None
    molecule_config: Optional[MoleculeConfig] = None
    logging_config: LoggingConfig = field(default_factory=LoggingConfig)
    fba_config: FBAConfig = field(default_factory=FBAConfig)
    fermentation_config: FermentationConfig = field(default_factory=FermentationConfig)
    dbtl_config: DBTLConfig = field(default_factory=DBTLConfig)
    scaleup_config: ScaleUpConfig = field(default_factory=ScaleUpConfig)
    downstream_config: DownstreamConfig = field(default_factory=DownstreamConfig)
    regulatory_config: RegulatoryConfig = field(default_factory=RegulatoryConfig)
    
    # Output settings
    output_dir: str = "./pipeline_output"
    random_seed: int = 42
    stop_on_error: bool = False
    
    @classmethod
    def from_args(
        cls,
        organism: str = "ecoli",
        molecule: str = "lycopene",
        cycles: int = 3,
        output_dir: str = "./pipeline_output"
    ) -> 'PipelineConfig':
        """Create configuration from command-line arguments."""
        
        # Map organism string to config
        organism_map = {
            "ecoli": OrganismConfig.get_ecoli_config,
            "ecoli_bl21": OrganismConfig.get_ecoli_bl21_config,
            "scerevisiae": OrganismConfig.get_scerevisiae_config,
            "bsubtilis": OrganismConfig.get_bsubtilis_config,
            "cglutamicum": OrganismConfig.get_cglutamicum_config,
            "pputida": OrganismConfig.get_pputida_config
        }
        
        # Map molecule string to config
        molecule_map = {
            "lycopene": MoleculeConfig.get_lycopene_config,
            "vanillin": MoleculeConfig.get_vanillin_config,
            "artemisinic_acid": MoleculeConfig.get_artemisinic_acid_config,
            "l_lysine": MoleculeConfig.get_l_lysine_config,
            "l_glutamate": MoleculeConfig.get_l_glutamate_config,
            "l_threonine": MoleculeConfig.get_l_threonine_config,
            "pha": MoleculeConfig.get_pha_config,
            "hyaluronic_acid": MoleculeConfig.get_hyaluronic_acid_config,
            "riboflavin": MoleculeConfig.get_riboflavin_config
        }
        
        org_config_func = organism_map.get(organism.lower(), OrganismConfig.get_ecoli_config)
        mol_config_func = molecule_map.get(molecule.lower(), MoleculeConfig.get_lycopene_config)
        
        return cls(
            organism_type=OrganismType(organism.lower()),
            molecule_type=MoleculeType(molecule.lower()),
            organism_config=org_config_func(),
            molecule_config=mol_config_func(),
            dbtl_config=DBTLConfig(num_cycles=cycles),
            output_dir=output_dir
        )
    
    def to_dict(self) -> Dict:
        """Convert configuration to dictionary for JSON serialization."""
        return {
            "organism_type": self.organism_type.value,
            "molecule_type": self.molecule_type.value,
            "output_dir": self.output_dir,
            "random_seed": self.random_seed,
            "dbtl_cycles": self.dbtl_config.num_cycles
        }


if __name__ == "__main__":
    # Test configuration creation
    import json
    
    print("Testing PipelineConfig...")
    
    # Create default config
    config = PipelineConfig.from_args(
        organism="ecoli",
        molecule="lycopene",
        cycles=3
    )
    
    print(f"Organism: {config.organism_config.name} ({config.organism_config.strain})")
    print(f"Target molecule: {config.molecule_config.name}")
    print(f"Target titer: {config.molecule_config.target_titer_g_per_l} g/L")
    print(f"DBTL cycles: {config.dbtl_config.num_cycles}")
    
    # Test all organism configs
    print("\n--- Testing all organism configurations ---")
    for org_func in [
        OrganismConfig.get_ecoli_config,
        OrganismConfig.get_scerevisiae_config,
        OrganismConfig.get_bsubtilis_config,
        OrganismConfig.get_cglutamicum_config,
        OrganismConfig.get_pputida_config
    ]:
        org = org_func()
        print(f"✓ {org.name} {org.strain}: μmax={org.mu_max_per_hour}/h, GC={org.gc_content_percent}%")
    
    # Test all molecule configs
    print("\n--- Testing all molecule configurations ---")
    for mol_func in [
        MoleculeConfig.get_lycopene_config,
        MoleculeConfig.get_vanillin_config,
        MoleculeConfig.get_l_lysine_config,
        MoleculeConfig.get_pha_config,
        MoleculeConfig.get_riboflavin_config
    ]:
        mol = mol_func()
        print(f"✓ {mol.name}: target={mol.target_titer_g_per_l} g/L, MW={mol.molecular_weight}")
    
    print("\n✓ All configuration tests passed!")
