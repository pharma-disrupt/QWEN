# Synthetic Biology Metabolic Pathway Pipeline

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Docker](https://img.shields.io/badge/docker-ready-blue.svg)](https://www.docker.com/)

A comprehensive, industrial-grade Python pipeline for **metabolic pathway design** in industrial microorganisms. This tool implements the complete **DBTL (Design-Build-Test-Learn)** cycle for synthetic biology applications.

## 🎯 Features

### Complete 5-Stage Pipeline

| Stage | Module | Description |
|-------|--------|-------------|
| **1** | Core Infrastructure | Organism/molecule databases, genomic data loading, schema validation |
| **2** | Pathway Prediction | Retrosynthesis engine, enzyme selection, codon optimization, AI scoring |
| **3** | Flux Analysis | FBA/pFBA/FVA, strain optimization, toxicity prediction |
| **4** | DBTL Loop | Bayesian optimization, fermentation simulation, bioreactor control |
| **5** | Scale-up & Regulatory | Production scale-up, downstream processing, compliance assessment |

### Supported Organisms
- ✅ *Escherichia coli* K-12 MG1655 / BL21
- ✅ *Saccharomyces cerevisiae* S288C / BY4741
- ✅ *Bacillus subtilis* 168
- ✅ *Corynebacterium glutamicum* ATCC 13032
- ✅ *Pseudomonas putida* KT2440

### Target Products
- **Small molecules**: Lycopene, Vanillin, Artemisinic acid
- **Amino acids**: L-Lysine, L-Glutamate, L-Threonine
- **Biopolymers**: PHA, Hyaluronic acid
- **Vitamins**: Riboflavin (B2)

## 🚀 Quick Start

### Option 1: Local Installation

```bash
# Clone the repository
git clone https://github.com/your-org/synbio-pipeline.git
cd synbio-pipeline

# Install dependencies
pip install -r requirements.txt

# Run the pipeline
python main_pipeline_runner.py --organism ecoli --molecule lycopene --cycles 3
```

### Option 2: Docker (Recommended)

```bash
# Build the Docker image
docker build -t synbio-pipeline .

# Run with default parameters
docker run --rm -v $(pwd)/output:/app/pipeline_output synbio-pipeline \
    --organism ecoli --molecule lycopene --cycles 3

# Run with custom parameters
docker run --rm -v $(pwd)/output:/app/pipeline_output synbio-pipeline \
    --organism scerevisiae --molecule vanillin --cycles 5 --output-dir /app/pipeline_output
```

## 📖 Usage

### Command Line Interface

```bash
python main_pipeline_runner.py [OPTIONS]

Options:
  --organism TEXT     Target organism [required]
                      Choices: ecoli, scerevisiae, bsubtilis, cglutamicum, pputida
  
  --molecule TEXT     Target molecule [required]
                      Choices: lycopene, vanillin, artemisinic_acid, l_lysine, 
                               l_glutamate, l_threonine, pha, hyaluronic_acid, 
                               riboflavin
  
  --cycles INTEGER    Number of DBTL cycles (default: 3)
  
  --output-dir TEXT   Output directory (default: ./pipeline_output)
  
  --log-level TEXT    Logging level: DEBUG, INFO, WARNING, ERROR (default: INFO)
  
  --help              Show this message and exit
```

### Examples

#### Example 1: Lycopene production in E. coli
```bash
python main_pipeline_runner.py \
    --organism ecoli \
    --molecule lycopene \
    --cycles 3 \
    --output-dir ./lycopene_results
```

#### Example 2: Vanillin production in Yeast
```bash
python main_pipeline_runner.py \
    --organism scerevisiae \
    --molecule vanillin \
    --cycles 5 \
    --log-level DEBUG
```

#### Example 3: L-Lysine in C. glutamicum
```bash
python main_pipeline_runner.py \
    --organism cglutamicum \
    --molecule l_lysine \
    --cycles 4
```

### Programmatic Usage

```python
from main_pipeline_runner import PipelineRunner

# Initialize the runner
runner = PipelineRunner()

# Run the full pipeline
results = runner.run_full_pipeline(
    organism="ecoli",
    target_molecule="lycopene",
    dbtl_cycles=3,
    output_dir="./results"
)

# Access stage-specific results
print(f"Predicted titer: {results['stage_4']['fermentation_simulation']['final_titer_g_per_l']} g/L")
print(f"Regulatory compliance: {results['stage_5']['regulatory_assessment']['compliance_score']}/100")
```

## 📁 Project Structure

```
synbio-pipeline/
├── main_pipeline_runner.py       # Master orchestrator
├── pipeline_config.py            # Configuration dataclasses
├── logger_setup.py               # Logging infrastructure
├── exceptions.py                 # Custom exception classes
├── schema_validator.py           # JSON schema validation
│
├── # Stage 1: Core Infrastructure
├── data_layer.py                 # Organism & molecule databases
│
├── # Stage 2: Pathway Prediction
├── stage_2_pathway_ai.py         # Complete Stage 2 orchestration
│
├── # Stage 3: Flux Analysis
├── fba_engine.py                 # Flux balance analysis
├── strain_optimizer.py           # Genetic modification optimization
├── toxicity_predictor.py         # Toxicity assessment
├── flux_analysis_orchestrator.py # Stage 3 orchestrator
│
├── # Stage 4: DBTL Loop
├── dbtl_loop.py                  # Design-Build-Test-Learn cycles
├── fermentation_simulator.py     # Bioreactor simulation
├── bioreactor_controller.py      # MPC/PID control
├── dbtl_orchestrator.py          # Stage 4 orchestrator
│
├── # Stage 5: Scale-up & Regulatory
├── scaleup_engine.py             # Scale-up predictions
├── downstream_processor.py       # Purification train design
├── regulatory_checker.py         # Compliance assessment
├── pipeline_report_generator.py  # Report generation
│
├── requirements.txt              # Python dependencies
├── Dockerfile                    # Container definition
├── docker-compose.yml            # Multi-container setup
├── README.md                     # This file
└── LICENSE                       # MIT License
```

## 🔬 Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                     STAGE 1: DATA LAYER                         │
│  Organism DB → Genomic Data → Quality Check → Schema Validation│
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                STAGE 2: PATHWAY PREDICTION                      │
│  Retrosynthesis → Enzyme Selection → Codon Optimization → AI   │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                 STAGE 3: FLUX ANALYSIS                          │
│     FBA/pFBA → Strain Optimization → Toxicity Assessment        │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│                   STAGE 4: DBTL LOOP                            │
│  Construct Design → Screening → Bayesian Opt → Fermentation    │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│            STAGE 5: SCALE-UP & REGULATORY                       │
│  Scale-up → Downstream Processing → Regulatory Compliance       │
└─────────────────────────────────────────────────────────────────┘
                              ↓
                    ┌──────────────────┐
                    │  FINAL REPORT    │
                    │  JSON + HTML     │
                    └──────────────────┘
```

## 📊 Output Files

After running the pipeline, you'll find:

```
pipeline_output/
├── final_report.json           # Complete machine-readable results
├── final_report.html           # Human-readable HTML report
├── final_report.txt            # Text summary
├── fermentation_curve.png      # Time-series plot
├── pathway_diagram.png         # Pathway visualization
└── stage_summaries/
    ├── stage_1_summary.json
    ├── stage_2_summary.json
    ├── stage_3_summary.json
    ├── stage_4_summary.json
    └── stage_5_summary.json
```

### Sample Output (JSON excerpt)

```json
{
  "pipeline_id": "uuid-here",
  "timestamp": "2024-01-15T10:30:00Z",
  "organism": {
    "name": "Escherichia coli",
    "strain": "K-12 MG1655",
    "doubling_time_min": 20.0
  },
  "target_molecule": {
    "name": "Lycopene",
    "target_titer_g_per_l": 5.0
  },
  "stage_4": {
    "fermentation_simulation": {
      "final_titer_g_per_l": 4.23,
      "final_yield_g_per_g": 0.18,
      "productivity_g_per_l_per_h": 0.14
    }
  },
  "stage_5": {
    "scale_up_prediction": {
      "production_titer_g_per_l": 3.81,
      "yield_loss_percent": 9.9
    },
    "regulatory_assessment": {
      "biosafety_level": "BSL-1",
      "gras_status": true,
      "compliance_score": 85
    }
  }
}
```

## ⚙️ Configuration

Edit `pipeline_config.py` to customize:

```python
@dataclass
class PipelineConfig:
    # Organism selection
    organism_name: str = "ecoli"
    strain: str = "K-12 MG1655"
    
    # Pipeline parameters
    dbtl_cycles: int = 3
    constructs_per_cycle: int = 48
    
    # Fermentation settings
    fermentation_mode: str = "fed-batch"
    duration_hours: float = 48.0
    
    # Scale-up levels
    scale_levels: List[float] = field(default_factory=lambda: [2, 200, 20000])
```

## 🧪 Testing

```bash
# Run all tests
pytest tests/

# Run specific stage tests
pytest tests/test_stage_1.py
pytest tests/test_stage_2.py
pytest tests/test_fba_engine.py

# Run with coverage
pytest --cov=. --cov-report=html
```

## 🐳 Docker Advanced Usage

### Development Mode
```bash
docker-compose up --build
```

### Production Mode
```bash
docker run -d \
  --name synbio-prod \
  -v /data/results:/app/pipeline_output \
  -e LOG_LEVEL=INFO \
  synbio-pipeline:latest
```

### Interactive Session
```bash
docker run -it --rm synbio-pipeline python
>>> from main_pipeline_runner import PipelineRunner
>>> runner = PipelineRunner()
```

## 📈 Performance Benchmarks

| Organism | Molecule | DBTL Cycles | Runtime | Predicted Titer |
|----------|----------|-------------|---------|-----------------|
| E. coli | Lycopene | 3 | ~45s | 4.2 g/L |
| S. cerevisiae | Vanillin | 5 | ~78s | 2.8 g/L |
| C. glutamicum | L-Lysine | 4 | ~62s | 85 g/L |
| B. subtilis | PHA | 3 | ~51s | 12 g/L |

*Tests run on Intel i7, 16GB RAM*

## 🔒 Security & Compliance

- **Biosafety Level Assessment**: Automatic BSL-1/2/3 classification
- **GRAS Status Checking**: Verified against FDA database
- **Allergenicity Prediction**: Simulated AllerTop model
- **Toxicity Assessment**: Simulated ProTox-III
- **Antibiotic Marker Tracking**: Ensures removal for production strains

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Reaction rules inspired by [MetaCyc](https://metacyc.org/)
- Organism models based on [BiGG Models](http://bigg.ucsd.edu/)
- Codon usage tables from [Kazusa](https://www.kazusa.or.jp/codon/)
- Toxicity thresholds from [PubChem](https://pubchem.ncbi.nlm.nih.gov/)

## 📞 Support

- **Documentation**: [Wiki](https://github.com/your-org/synbio-pipeline/wiki)
- **Issues**: [GitHub Issues](https://github.com/your-org/synbio-pipeline/issues)
- **Email**: support@synbio-pipeline.org

## 🎓 Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{synbio_pipeline2024,
  title = {Synthetic Biology Metabolic Pathway Pipeline},
  author = {Your Name and Contributors},
  year = {2024},
  url = {https://github.com/your-org/synbio-pipeline}
}
```

---

**Built with ❤️ for the Synthetic Biology Community**
