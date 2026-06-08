# Synthetic Biology Pipeline - Makefile

.PHONY: help install test lint type-check docker-build docker-run clean all

## Help
help: ## Show this help message
	@echo "Usage: make [target]"
	@echo ""
	@echo "Available targets:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  %-15s %s\n", $$1, $$2}'

## Installation
install: ## Install all dependencies
	pip install -r requirements.txt

install-dev: ## Install development dependencies
	pip install -r requirements.txt
	pip install jupyterlab ipywidgets plotly seaborn

## Testing
test: ## Run all tests
	pytest tests/ -v

test-cov: ## Run tests with coverage
	pytest tests/ -v --cov=. --cov-report=html --cov-report=term-missing

test-single: ## Run a single test file (usage: make test-single TEST_FILE=tests/test_fba.py)
	pytest $(TEST_FILE) -v

## Code Quality
lint: ## Run flake8 linter
	flake8 *.py --count --select=E9,F63,F7,F82 --show-source --statistics
	flake8 *.py --count --exit-zero --max-complexity=10 --max-line-length=127 --statistics

type-check: ## Run mypy type checker
	mypy *.py --ignore-missing-imports

format: ## Format code with black
	black *.py

code-quality: lint type-check ## Run all code quality checks

## Docker
docker-build: ## Build Docker image
	docker build -t synbio-pipeline:latest .

docker-build-no-cache: ## Build Docker image without cache
	docker build --no-cache -t synbio-pipeline:latest .

docker-run: ## Run pipeline in Docker (default: ecoli + lycopene)
	docker run --rm -v $(PWD)/pipeline_output:/app/pipeline_output synbio-pipeline:latest \
		--organism ecoli --molecule lycopene --cycles 3

docker-run-custom: ## Run custom pipeline in Docker (usage: make docker-run-custom ORG=ecoli MOL=lycopene)
	docker run --rm -v $(PWD)/pipeline_output:/app/pipeline_output synbio-pipeline:latest \
		--organism $(ORG) --molecule $(MOL) --cycles $(CYCLES)

docker-compose-up: ## Start all services with docker-compose
	docker-compose up --build

docker-compose-up-detached: ## Start all services in detached mode
	docker-compose up -d --build

docker-compose-down: ## Stop all services
	docker-compose down

docker-compose-full: ## Start full stack (with PostgreSQL and Redis)
	docker-compose --profile full up --build

docker-jupyter: ## Start Jupyter notebook environment
	docker-compose --profile dev up jupyter

docker-clean: ## Remove all Docker containers and images
	docker system prune -af

## Pipeline Execution
run-ecoli-lycopene: ## Run pipeline for E. coli + Lycopene
	python main_pipeline_runner.py --organism ecoli --molecule lycopene --cycles 3

run-yeast-vanillin: ## Run pipeline for S. cerevisiae + Vanillin
	python main_pipeline_runner.py --organism scerevisiae --molecule vanillin --cycles 5

run-coryne-lysine: ## Run pipeline for C. glutamicum + L-Lysine
	python main_pipeline_runner.py --organism cglutamicum --molecule l_lysine --cycles 4

run-all-organisms: ## Run pipeline for all organisms with lycopene
	python main_pipeline_runner.py --organism ecoli --molecule lycopene --cycles 2
	python main_pipeline_runner.py --organism scerevisiae --molecule lycopene --cycles 2
	python main_pipeline_runner.py --organism bsubtilis --molecule lycopene --cycles 2
	python main_pipeline_runner.py --organism cglutamicum --molecule lycopene --cycles 2
	python main_pipeline_runner.py --organism pputida --molecule lycopene --cycles 2

## Cleanup
clean: ## Remove build artifacts and output
	rm -rf __pycache__/
	rm -rf *.egg-info
	rm -rf .pytest_cache/
	rm -rf htmlcov/
	rm -rf .coverage
	rm -rf pipeline_output/
	rm -rf logs/
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

clean-all: clean ## Remove everything including downloaded data
	rm -rf notebooks/*.ipynb_checkpoints

## Development
notebook: ## Start Jupyter notebook locally
	jupyter notebook --notebook-dir=./notebooks

dev: install-dev docker-compose-up-detached ## Setup development environment

## CI/CD
ci: code-quality test ## Run CI checks locally

## All
all: install lint test docker-build ## Run complete build process
