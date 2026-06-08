"""
Tests for Stage 1: Core Infrastructure + Data Layer
"""
import pytest
import json
from datetime import datetime
from pathlib import Path

from pipeline_config import PipelineConfig, OrganismConfig
from exceptions import PipelineError, DataIngestionError, SchemaValidationError
from schema_validator import validate_stage_output, STAGE_1_OUTPUT_SCHEMA
from data_layer import OrganismDatabase, MoleculeDatabase, run_stage_1


class TestPipelineConfig:
    """Test PipelineConfig dataclass"""
    
    def test_default_config(self):
        """Test default configuration values"""
        config = PipelineConfig()
        assert config.organism_name == "ecoli"
        assert config.dbtl_cycles == 3
        assert config.fermentation_mode == "fed-batch"
    
    def test_custom_config(self):
        """Test custom configuration"""
        config = PipelineConfig(
            organism_name="scerevisiae",
            dbtl_cycles=5,
            fermentation_mode="batch"
        )
        assert config.organism_name == "scerevisiae"
        assert config.dbtl_cycles == 5


class TestOrganismDatabase:
    """Test OrganismDatabase class"""
    
    def test_get_ecoli_config(self):
        """Test E. coli configuration retrieval"""
        db = OrganismDatabase()
        config = db.get_organism_config("ecoli")
        
        assert config is not None
        assert config.name == "Escherichia coli"
        assert config.doubling_time_min < 30  # E. coli doubles in ~20 min
    
    def test_get_all_organisms(self):
        """Test retrieval of all supported organisms"""
        db = OrganismDatabase()
        organisms = ["ecoli", "scerevisiae", "bsubtilis", "cglutamicum", "pputida"]
        
        for org in organisms:
            config = db.get_organism_config(org)
            assert config is not None, f"Failed to get config for {org}"


class TestMoleculeDatabase:
    """Test MoleculeDatabase class"""
    
    def test_get_lycopene_info(self):
        """Test lycopene molecule information"""
        db = MoleculeDatabase()
        molecule = db.get_molecule_info("lycopene")
        
        assert molecule is not None
        assert molecule.name == "Lycopene"
        assert len(molecule.smiles) > 0
    
    def test_get_all_molecules(self):
        """Test retrieval of all supported molecules"""
        db = MoleculeDatabase()
        molecules = ["lycopene", "vanillin", "l_lysine", "riboflavin", "pha"]
        
        for mol in molecules:
            info = db.get_molecule_info(mol)
            assert info is not None, f"Failed to get info for {mol}"


class TestSchemaValidation:
    """Test JSON schema validation"""
    
    def test_valid_stage1_output(self):
        """Test validation of valid Stage 1 output"""
        valid_data = {
            "pipeline_id": "test-uuid-123",
            "timestamp": datetime.now().isoformat(),
            "organism": {
                "name": "Escherichia coli",
                "strain": "K-12 MG1655",
                "gem_model_id": "iJO1366",
                "doubling_time_min": 20.0,
                "optimal_ph": 7.0,
                "optimal_temp_c": 37.0,
                "gram_stain": "negative"
            },
            "target_molecule": {
                "name": "Lycopene",
                "smiles": "CC(C)=CCC=C(C)C=CC=C(C)C=CC=C(C)C=CC=C(C)C=CC=C(C)C=CC=C(C)C=CC=C(C)C=CC=C(C)C=C",
                "chebi_id": "CHEBI:15956",
                "target_titer_g_per_l": 5.0,
                "target_yield_mol_per_mol": 0.5
            },
            "genomic_data": {
                "total_genes": 4145,
                "essential_genes": ["dnaA", "rpoB"],
                "available_promoters": ["lac", "tet", "ara"],
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
        
        result = validate_stage_output(valid_data, STAGE_1_OUTPUT_SCHEMA)
        assert result["valid"] is True
    
    def test_invalid_stage1_output_missing_field(self):
        """Test validation fails for missing required field"""
        invalid_data = {
            "pipeline_id": "test-uuid",
            # Missing timestamp and other required fields
        }
        
        result = validate_stage_output(invalid_data, STAGE_1_OUTPUT_SCHEMA)
        assert result["valid"] is False


class TestRunStage1:
    """Test run_stage_1 function"""
    
    def test_run_stage1_ecoli_lycopene(self):
        """Test Stage 1 execution for E. coli + Lycopene"""
        result = run_stage_1(
            organism_name="ecoli",
            target_molecule="lycopene",
            pipeline_id="test-pipeline-001"
        )
        
        assert result is not None
        assert result["pipeline_id"] == "test-pipeline-001"
        assert result["organism"]["name"] == "Escherichia coli"
        assert result["target_molecule"]["name"] == "Lycopene"
        assert result["stage_1_status"] in ["PASS", "WARN", "FAIL"]
    
    def test_run_stage1_all_organisms(self):
        """Test Stage 1 for all supported organisms"""
        organisms = ["ecoli", "scerevisiae", "bsubtilis", "cglutamicum", "pputida"]
        
        for org in organisms:
            result = run_stage_1(
                organism_name=org,
                target_molecule="lycopene",
                pipeline_id=f"test-{org}"
            )
            
            assert result is not None
            assert result["stage_1_status"] in ["PASS", "WARN", "FAIL"]


class TestExceptions:
    """Test custom exception classes"""
    
    def test_pipeline_error(self):
        """Test PipelineError exception"""
        error = PipelineError("Test error message")
        assert str(error) == "Test error message"
    
    def test_data_ingestion_error(self):
        """Test DataIngestionError exception"""
        error = DataIngestionError("Failed to load data", {"key": "value"})
        assert "Failed to load data" in str(error)
        assert error.input_json == {"key": "value"}
    
    def test_schema_validation_error(self):
        """Test SchemaValidationError exception"""
        error = SchemaValidationError("Invalid schema", ["error1", "error2"])
        assert "Invalid schema" in str(error)
        assert len(error.validation_errors) == 2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
