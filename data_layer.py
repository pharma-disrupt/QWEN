"""
Data Layer Module
Provides data loading and management for organisms, molecules, and genomic data.
"""

import json
import uuid
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple
from dataclasses import dataclass, asdict

from pipeline_config import (
    OrganismConfig, MoleculeConfig, PipelineConfig,
    OrganismType, MoleculeType
)
from logger_setup import setup_logger, PipelineLogger
from schema_validator import (
    validate_stage_output, create_sample_pipeline_id, create_timestamp
)
from exceptions import DataIngestionError, PipelineError


@dataclass
class Stage1Output:
    """Stage 1 output data structure."""
    pipeline_id: str
    timestamp: str
    organism: Dict[str, Any]
    target_molecule: Dict[str, Any]
    genomic_data: Dict[str, Any]
    data_quality_report: Dict[str, Any]
    stage_1_status: str
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return asdict(self)
    
    def to_json(self) -> str:
        """Convert to JSON string."""
        return json.dumps(self.to_dict(), indent=2, default=str)


class OrganismDatabase:
    """
    Database of organism configurations.
    
    Provides access to all supported industrial microorganisms.
    """
    
    _organisms: Dict[str, OrganismConfig] = {}
    
    def __init__(self):
        """Initialize the organism database."""
        self._load_all_organisms()
    
    def _load_all_organisms(self) -> None:
        """Load all organism configurations."""
        self._organisms = {
            "ecoli": OrganismConfig.get_ecoli_config(),
            "ecoli_bl21": OrganismConfig.get_ecoli_bl21_config(),
            "scerevisiae": OrganismConfig.get_scerevisiae_config(),
            "bsubtilis": OrganismConfig.get_bsubtilis_config(),
            "cglutamicum": OrganismConfig.get_cglutamicum_config(),
            "pputida": OrganismConfig.get_pputida_config()
        }
    
    def get_organism(self, organism_type: str) -> OrganismConfig:
        """
        Get organism configuration by type.
        
        Args:
            organism_type: String identifier (e.g., "ecoli", "scerevisiae")
        
        Returns:
            OrganismConfig instance
        
        Raises:
            DataIngestionError: If organism not found
        """
        organism_key = organism_type.lower().replace(" ", "_").replace(".", "")
        
        if organism_key not in self._organisms:
            available = list(self._organisms.keys())
            raise DataIngestionError(
                message=f"Organism '{organism_type}' not found. Available: {available}",
                data_source="OrganismDatabase"
            )
        
        return self._organisms[organism_key]
    
    def list_organisms(self) -> List[str]:
        """List all available organism identifiers."""
        return list(self._organisms.keys())
    
    def get_all_configs(self) -> Dict[str, OrganismConfig]:
        """Get all organism configurations."""
        return self._organisms.copy()


class MoleculeDatabase:
    """
    Database of target molecule configurations.
    
    Provides access to all supported bioproducts.
    """
    
    _molecules: Dict[str, MoleculeConfig] = {}
    
    def __init__(self):
        """Initialize the molecule database."""
        self._load_all_molecules()
    
    def _load_all_molecules(self) -> None:
        """Load all molecule configurations."""
        self._molecules = {
            "lycopene": MoleculeConfig.get_lycopene_config(),
            "vanillin": MoleculeConfig.get_vanillin_config(),
            "artemisinic_acid": MoleculeConfig.get_artemisinic_acid_config(),
            "l_lysine": MoleculeConfig.get_l_lysine_config(),
            "l_glutamate": MoleculeConfig.get_l_glutamate_config(),
            "l_threonine": MoleculeConfig.get_l_threonine_config(),
            "pha": MoleculeConfig.get_pha_config(),
            "hyaluronic_acid": MoleculeConfig.get_hyaluronic_acid_config(),
            "riboflavin": MoleculeConfig.get_riboflavin_config()
        }
    
    def get_molecule(self, molecule_name: str) -> MoleculeConfig:
        """
        Get molecule configuration by name.
        
        Args:
            molecule_name: String identifier (e.g., "lycopene", "vanillin")
        
        Returns:
            MoleculeConfig instance
        
        Raises:
            DataIngestionError: If molecule not found
        """
        molecule_key = molecule_name.lower().replace(" ", "_").replace("-", "_")
        
        if molecule_key not in self._molecules:
            available = list(self._molecules.keys())
            raise DataIngestionError(
                message=f"Molecule '{molecule_name}' not found. Available: {available}",
                data_source="MoleculeDatabase"
            )
        
        return self._molecules[molecule_key]
    
    def list_molecules(self) -> List[str]:
        """List all available molecule identifiers."""
        return list(self._molecules.keys())
    
    def get_all_configs(self) -> Dict[str, MoleculeConfig]:
        """Get all molecule configurations."""
        return self._molecules.copy()


class GenomicDataLoader:
    """
    Simulated genomic data loader.
    
    Simulates loading genomic data from external databases like KEGG, UniProt, NCBI.
    In a real implementation, this would connect to actual APIs.
    """
    
    # Simulated gene counts for different organisms
    GENE_COUNTS = {
        "ecoli": 4377,
        "scerevisiae": 6692,
        "bsubtilis": 4100,
        "cglutamicum": 3000,
        "pputida": 5434
    }
    
    # Simulated pathway associations
    PATHWAY_GENES = {
        "terpenoid": ["dxs", "dxr", "ispD", "ispE", "ispF", "ispG", "ispH", "ispA"],
        "aromatic": ["aroA", "aroB", "aroC", "aroD", "aroE", "aroK", "aroL"],
        "aspartate_family": ["asd", "dapA", "dapB", "dapC", "dapD", "dapE", "dapF", "lysA"],
        "tca_cycle_derived": ["gltA", "acnA", "acnB", "icd", "sucA", "sucB", "sucC", "sucD", "sdhA", "sdhB", "fumA", "mdh"],
        "fatty_acid_derived": ["fabA", "fabB", "fabD", "fabF", "fabG", "fabH", "fabI", "fabZ"],
        "ribulose_5phosphate": ["ribA", "ribB", "ribC", "ribD", "ribE", "ribF"]
    }
    
    def __init__(self, organism_config: OrganismConfig):
        """
        Initialize the genomic data loader.
        
        Args:
            organism_config: Configuration for the target organism
        """
        self.organism_config = organism_config
        self.organism_key = organism_config.name.lower().replace(" ", "_").split("_")[0][:10]
    
    def load_genomic_data(self) -> Dict[str, Any]:
        """
        Load simulated genomic data.
        
        Returns:
            Dictionary containing genomic data
        """
        # Determine organism key for gene count lookup
        org_key = None
        for key in self.GENE_COUNTS.keys():
            if key in self.organism_config.strain.lower() or key in self.organism_config.name.lower():
                org_key = key
                break
        
        if org_key is None:
            org_key = "ecoli"  # Default fallback
        
        total_genes = self.GENE_COUNTS.get(org_key, 4000)
        
        return {
            "total_genes": total_genes,
            "essential_genes": self.organism_config.essential_genes[:20],  # Top 20 essential genes
            "available_promoters": self.organism_config.available_promoters,
            "codon_table_id": self.organism_config.codon_table_id,
            "gc_content_percent": self.organism_config.gc_content_percent
        }
    
    def get_pathway_genes(self, pathway_type: str) -> List[str]:
        """
        Get genes associated with a specific pathway type.
        
        Args:
            pathway_type: Type of pathway (e.g., "terpenoid", "aromatic")
        
        Returns:
            List of gene names
        """
        return self.PATHWAY_GENES.get(pathway_type, [])
    
    def simulate_database_query(
        self,
        database: str,
        query_type: str,
        gene_name: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Simulate querying an external database.
        
        Args:
            database: Database name (KEGG, UniProt, NCBI)
            query_type: Type of query (gene_info, sequence, pathway)
            gene_name: Optional gene name to query
        
        Returns:
            Simulated query result
        """
        if database == "KEGG":
            return self._simulate_kegg_query(query_type, gene_name)
        elif database == "UniProt":
            return self._simulate_uniprot_query(query_type, gene_name)
        elif database == "NCBI":
            return self._simulate_ncbi_query(query_type, gene_name)
        else:
            return {"error": f"Unknown database: {database}"}
    
    def _simulate_kegg_query(self, query_type: str, gene_name: Optional[str]) -> Dict[str, Any]:
        """Simulate KEGG database query."""
        if query_type == "pathway":
            return {
                "pathways": [
                    {"id": "map00900", "name": "Terpenoid backbone biosynthesis"},
                    {"id": "map00400", "name": "Phenylalanine, tyrosine and tryptophan biosynthesis"},
                    {"id": "map00300", "name": "Lysine biosynthesis"}
                ]
            }
        elif query_type == "gene" and gene_name:
            return {
                "gene_id": f"b{hash(gene_name) % 10000:04d}",
                "name": gene_name,
                "organism": self.organism_config.strain,
                "pathway": "metabolism"
            }
        return {}
    
    def _simulate_uniprot_query(self, query_type: str, gene_name: Optional[str]) -> Dict[str, Any]:
        """Simulate UniProt database query."""
        if gene_name:
            return {
                "accession": f"P{hash(gene_name) % 100000:05d}",
                "entry_name": gene_name.upper() + "_ECOLI",
                "protein_name": f"{gene_name} protein",
                "organism": self.organism_config.name
            }
        return {}
    
    def _simulate_ncbi_query(self, query_type: str, gene_name: Optional[str]) -> Dict[str, Any]:
        """Simulate NCBI database query."""
        if gene_name:
            return {
                "gene_id": hash(gene_name) % 1000000,
                "symbol": gene_name,
                "description": f"{gene_name} gene",
                "chromosome": "1",
                "genomic_region": f"{hash(gene_name) % 1000000}-{hash(gene_name) % 1000000 + 1000}"
            }
        return {}


class DataQualityChecker:
    """
    Data quality checker for genomic and metabolic data.
    
    Validates completeness and consistency of loaded data.
    """
    
    def __init__(self, logger: Optional[PipelineLogger] = None):
        """
        Initialize the data quality checker.
        
        Args:
            logger: Optional pipeline logger
        """
        self.logger = logger
        self.warnings: List[str] = []
        self.errors: List[str] = []
    
    def check_data_quality(
        self,
        organism_data: Dict[str, Any],
        molecule_data: Dict[str, Any],
        genomic_data: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Check overall data quality.
        
        Args:
            organism_data: Organism configuration data
            molecule_data: Molecule configuration data
            genomic_data: Genomic data
        
        Returns:
            Quality report dictionary
        """
        self.warnings = []
        self.errors = []
        
        # Check organism data
        self._check_organism_data(organism_data)
        
        # Check molecule data
        self._check_molecule_data(molecule_data)
        
        # Check genomic data
        self._check_genomic_data(genomic_data)
        
        # Calculate completeness score
        completeness_score = self._calculate_completeness()
        
        # Determine status
        if self.errors:
            status = "FAIL"
        elif self.warnings:
            status = "WARN"
        else:
            status = "PASS"
        
        report = {
            "completeness_score": completeness_score,
            "warnings": self.warnings,
            "errors": self.errors
        }
        
        if self.logger:
            self.logger.debug(f"Data quality report: {report}")
        
        return report
    
    def _check_organism_data(self, data: Dict[str, Any]) -> None:
        """Check organism data quality."""
        required_fields = ["name", "strain", "gem_model_id", "doubling_time_min"]
        
        for field in required_fields:
            if field not in data or not data[field]:
                self.errors.append(f"Missing required organism field: {field}")
        
        # Check for realistic values
        if "doubling_time_min" in data:
            dt = data["doubling_time_min"]
            if dt < 10 or dt > 300:
                self.warnings.append(f"Unusual doubling time: {dt} min")
        
        if "optimal_ph" in data:
            ph = data["optimal_ph"]
            if ph < 3 or ph > 10:
                self.warnings.append(f"Extreme optimal pH: {ph}")
        
        if "optimal_temp_c" in data:
            temp = data["optimal_temp_c"]
            if temp < 15 or temp > 50:
                self.warnings.append(f"Extreme optimal temperature: {temp}°C")
    
    def _check_molecule_data(self, data: Dict[str, Any]) -> None:
        """Check molecule data quality."""
        required_fields = ["name", "smiles", "chebi_id", "target_titer_g_per_l"]
        
        for field in required_fields:
            if field not in data or not data[field]:
                self.errors.append(f"Missing required molecule field: {field}")
        
        # Check SMILES format (basic check)
        if "smiles" in data:
            smiles = data["smiles"]
            if not any(c in smiles for c in "CCOHN"):
                self.warnings.append(f"SMILES may be invalid: {smiles[:50]}")
        
        # Check target titer realism
        if "target_titer_g_per_l" in data:
            titer = data["target_titer_g_per_l"]
            if titer > 200:
                self.warnings.append(f"Very high target titer: {titer} g/L")
            elif titer < 0.1:
                self.warnings.append(f"Very low target titer: {titer} g/L")
    
    def _check_genomic_data(self, data: Dict[str, Any]) -> None:
        """Check genomic data quality."""
        required_fields = ["total_genes", "essential_genes", "available_promoters"]
        
        for field in required_fields:
            if field not in data or not data[field]:
                self.errors.append(f"Missing required genomic field: {field}")
        
        # Check gene count
        if "total_genes" in data:
            count = data["total_genes"]
            if count < 500 or count > 10000:
                self.warnings.append(f"Unusual gene count: {count}")
        
        # Check essential genes list
        if "essential_genes" in data:
            if len(data["essential_genes"]) < 10:
                self.warnings.append(f"Few essential genes listed: {len(data['essential_genes'])}")
        
        # Check promoter availability
        if "available_promoters" in data:
            if len(data["available_promoters"]) < 3:
                self.warnings.append(f"Limited promoter options: {len(data['available_promoters'])}")
    
    def _calculate_completeness(self) -> float:
        """Calculate data completeness score."""
        total_checks = 10  # Base number of checks
        passed_checks = total_checks - len(self.errors) - (len(self.warnings) * 0.3)
        
        score = max(0.0, min(1.0, passed_checks / total_checks))
        return round(score, 2)


def run_stage_1(
    organism: str,
    molecule: str,
    log_dir: str = "logs"
) -> Dict[str, Any]:
    """
    Execute Stage 1: Core Infrastructure + Data Layer.
    
    This function orchestrates all Stage 1 components:
    1. Load organism configuration
    2. Load molecule configuration
    3. Load genomic data
    4. Perform data quality checks
    5. Generate Stage 1 output JSON
    
    Args:
        organism: Organism name (e.g., 'ecoli', 'scerevisiae')
        molecule: Target molecule (e.g., 'lycopene', 'vanillin')
        log_dir: Directory for log files
    
    Returns:
        Stage 1 output dictionary matching STAGE_1_OUTPUT_SCHEMA
    
    Raises:
        DataIngestionError: If data loading fails
        PipelineError: If validation fails
    """
    import time
    start_time = time.time()
    
    # Generate pipeline ID
    pipeline_id = str(uuid.uuid4())
    
    # Create config from arguments
    config = PipelineConfig.from_args(organism=organism, molecule=molecule)
    
    # Set up logger
    logger = setup_logger(pipeline_id, stage=1, log_dir=log_dir)
    
    try:
        logger.log_stage_start(
            organism=config.organism_config.name,
            molecule=config.molecule_config.name,
            config=config.to_dict()
        )
        
        # Initialize databases
        logger.info("Initializing organism and molecule databases...")
        organism_db = OrganismDatabase()
        molecule_db = MoleculeDatabase()
        
        # Load organism data
        logger.info(f"Loading organism data for {config.organism_config.strain}...")
        organism_data = {
            "name": config.organism_config.name,
            "strain": config.organism_config.strain,
            "gem_model_id": config.organism_config.gem_model_id,
            "doubling_time_min": config.organism_config.doubling_time_min,
            "optimal_ph": config.organism_config.optimal_ph,
            "optimal_temp_c": config.organism_config.optimal_temp_c,
            "gram_stain": config.organism_config.gram_stain,
            "codon_table_id": config.organism_config.codon_table_id,
            "gc_content_percent": config.organism_config.gc_content_percent
        }
        logger.debug(f"Organism data loaded: {organism_data}")
        
        # Load molecule data
        logger.info(f"Loading molecule data for {config.molecule_config.name}...")
        molecule_data = {
            "name": config.molecule_config.name,
            "smiles": config.molecule_config.smiles,
            "chebi_id": config.molecule_config.chebi_id,
            "target_titer_g_per_l": config.molecule_config.target_titer_g_per_l,
            "target_yield_mol_per_mol": config.molecule_config.target_yield_mol_per_mol,
            "molecular_weight": config.molecule_config.molecular_weight,
            "category": config.molecule_config.category,
            "toxicity_level": config.molecule_config.toxicity_level,
            "secretion_type": config.molecule_config.secretion_type
        }
        logger.debug(f"Molecule data loaded: {molecule_data}")
        
        # Load genomic data
        logger.info("Loading genomic data...")
        genomic_loader = GenomicDataLoader(config.organism_config)
        genomic_data = genomic_loader.load_genomic_data()
        logger.debug(f"Genomic data loaded: {genomic_data}")
        
        # Perform data quality checks
        logger.info("Performing data quality checks...")
        quality_checker = DataQualityChecker(logger)
        quality_report = quality_checker.check_data_quality(
            organism_data=organism_data,
            molecule_data=molecule_data,
            genomic_data=genomic_data
        )
        logger.info(f"Data quality score: {quality_report['completeness_score']}")
        
        # Determine stage status
        if quality_report["errors"]:
            stage_status = "FAIL"
        elif quality_report["warnings"]:
            stage_status = "WARN"
        else:
            stage_status = "PASS"
        
        # Create Stage 1 output
        output = Stage1Output(
            pipeline_id=pipeline_id,
            timestamp=create_timestamp(),
            organism=organism_data,
            target_molecule=molecule_data,
            genomic_data=genomic_data,
            data_quality_report=quality_report,
            stage_1_status=stage_status
        )
        
        output_dict = output.to_dict()
        
        # Validate output against schema
        logger.info("Validating output against Stage 1 schema...")
        is_valid, errors = validate_stage_output(
            output_dict,
            "stage_1_output",
            raise_on_error=False
        )
        
        if not is_valid:
            logger.warning(f"Schema validation warnings: {errors}")
        
        # Log JSON contract
        logger.log_json_contract(output_dict, "stage_1_output", "output")
        
        # Calculate duration
        duration = time.time() - start_time
        
        # Log completion
        logger.log_stage_complete(
            status=stage_status,
            duration_seconds=duration,
            output_summary={
                "pipeline_id": pipeline_id,
                "organism": organism_data["name"],
                "molecule": molecule_data["name"],
                "quality_score": quality_report["completeness_score"]
            }
        )
        
        # Save stage summary
        summary = {
            "pipeline_id": pipeline_id,
            "stage": 1,
            "status": stage_status,
            "duration_seconds": round(duration, 3),
            "organism": organism_data,
            "molecule": molecule_data,
            "quality_report": quality_report
        }
        logger.save_stage_summary(summary)
        
        logger.info(f"Stage 1 complete. Output saved to logs/stage_1_summary.json")
        
        return output_dict
        
    except DataIngestionError as e:
        logger.log_error_summary(
            exception=e,
            function_name="run_stage_1",
            input_data=config.to_dict(),
            recovery_action=None
        )
        raise
    
    except Exception as e:
        error_msg = f"Unexpected error in Stage 1: {str(e)}"
        logger.error(error_msg)
        raise PipelineError(
            message=error_msg,
            stage=1,
            function_name="run_stage_1",
            input_json=config.to_dict(),
            original_exception=e
        )


if __name__ == "__main__":
    # Test Stage 1 execution
    print("Testing Stage 1: Core Infrastructure + Data Layer\n")
    
    # Create test configuration
    config = PipelineConfig.from_args(
        organism="ecoli",
        molecule="lycopene",
        cycles=3
    )
    
    # Run Stage 1
    try:
        output = run_stage_1(config)
        
        print("\n✓ Stage 1 completed successfully!")
        print(f"Pipeline ID: {output['pipeline_id']}")
        print(f"Organism: {output['organism']['name']} ({output['organism']['strain']})")
        print(f"Target molecule: {output['target_molecule']['name']}")
        print(f"Data quality score: {output['data_quality_report']['completeness_score']}")
        print(f"Stage status: {output['stage_1_status']}")
        
        if output['data_quality_report']['warnings']:
            print(f"\nWarnings: {output['data_quality_report']['warnings']}")
        if output['data_quality_report']['errors']:
            print(f"\nErrors: {output['data_quality_report']['errors']}")
        
    except PipelineError as e:
        print(f"\n✗ Stage 1 failed: {e.message}")
        print(f"Error details: {e.to_dict()}")
