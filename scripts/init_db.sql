-- Initialize database for storing pipeline results
-- Run automatically when PostgreSQL container starts

CREATE TABLE IF NOT EXISTS pipeline_runs (
    id SERIAL PRIMARY KEY,
    pipeline_id VARCHAR(255) UNIQUE NOT NULL,
    timestamp TIMESTAMPTZ DEFAULT CURRENT_TIMESTAMP,
    organism_name VARCHAR(100) NOT NULL,
    organism_strain VARCHAR(100),
    target_molecule VARCHAR(100) NOT NULL,
    overall_status VARCHAR(20),
    stage_1_status VARCHAR(20),
    stage_2_status VARCHAR(20),
    stage_3_status VARCHAR(20),
    stage_4_status VARCHAR(20),
    stage_5_status VARCHAR(20)
);

CREATE TABLE IF NOT EXISTS fermentation_results (
    id SERIAL PRIMARY KEY,
    pipeline_id VARCHAR(255) REFERENCES pipeline_runs(pipeline_id),
    final_titer_g_per_l DECIMAL(10, 4),
    final_yield_g_per_g DECIMAL(10, 4),
    final_productivity_g_per_l_per_h DECIMAL(10, 4),
    duration_hours DECIMAL(10, 2),
    mode VARCHAR(50)
);

CREATE TABLE IF NOT EXISTS scaleup_predictions (
    id SERIAL PRIMARY KEY,
    pipeline_id VARCHAR(255) REFERENCES pipeline_runs(pipeline_id),
    lab_scale_titer DECIMAL(10, 4),
    pilot_scale_titer DECIMAL(10, 4),
    production_scale_titer DECIMAL(10, 4),
    yield_loss_percent DECIMAL(5, 2),
    risk_level VARCHAR(20)
);

CREATE TABLE IF NOT EXISTS regulatory_assessments (
    id SERIAL PRIMARY KEY,
    pipeline_id VARCHAR(255) REFERENCES pipeline_runs(pipeline_id),
    biosafety_level VARCHAR(10),
    gras_status BOOLEAN,
    compliance_score INTEGER,
    allergenicity_risk VARCHAR(20),
    toxicity_risk VARCHAR(20)
);

CREATE TABLE IF NOT EXISTS dbtl_cycles (
    id SERIAL PRIMARY KEY,
    pipeline_id VARCHAR(255) REFERENCES pipeline_runs(pipeline_id),
    cycle_number INTEGER,
    constructs_tested INTEGER,
    best_titer_g_per_l DECIMAL(10, 4),
    best_construct_id VARCHAR(255),
    improvement_fold DECIMAL(6, 2)
);

-- Create indexes for faster queries
CREATE INDEX idx_pipeline_runs_organism ON pipeline_runs(organism_name);
CREATE INDEX idx_pipeline_runs_molecule ON pipeline_runs(target_molecule);
CREATE INDEX idx_pipeline_runs_timestamp ON pipeline_runs(timestamp);
CREATE INDEX idx_fermentation_titer ON fermentation_results(final_titer_g_per_l);

-- Grant permissions
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA public TO synbio;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA public TO synbio;

COMMENT ON TABLE pipeline_runs IS 'Stores metadata for each pipeline execution';
COMMENT ON TABLE fermentation_results IS 'Stores fermentation simulation results';
COMMENT ON TABLE scaleup_predictions IS 'Stores scale-up predictions from lab to production';
COMMENT ON TABLE regulatory_assessments IS 'Stores regulatory compliance assessments';
COMMENT ON TABLE dbtl_cycles IS 'Stores DBTL cycle iteration data';
