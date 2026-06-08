"""
Pipeline Report Generator for Industrial Metabolic Pathway Pipeline.
Compiles all stage outputs into comprehensive reports (JSON, HTML, text).
"""

import logging
import json
import os
from datetime import datetime
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Any, Optional
from pathlib import Path

# Internal imports
from pipeline_config import PipelineConfig
from exceptions import PipelineError
import logging

def get_logger(name: str):
    """Get a standard Python logger for the given module name."""
    return logging.getLogger(name)

logger = get_logger(__name__)

@dataclass
class PipelineReport:
    """Complete pipeline execution report."""
    pipeline_id: str
    timestamp: str
    organism: str
    target_molecule: str
    status: str  # SUCCESS, PARTIAL, FAILED
    
    # Stage summaries
    stage_1_summary: Optional[Dict[str, Any]] = None
    stage_2_summary: Optional[Dict[str, Any]] = None
    stage_3_summary: Optional[Dict[str, Any]] = None
    stage_4_summary: Optional[Dict[str, Any]] = None
    stage_5_summary: Optional[Dict[str, Any]] = None
    
    # Aggregated metrics
    final_predicted_titer_g_per_l: float = 0.0
    final_production_titer_g_per_l: float = 0.0
    downstream_cost_per_kg_usd: float = 0.0
    regulatory_feasibility: str = "UNKNOWN"
    overall_feasibility_score: float = 0.0
    
    # Errors and warnings
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    
    # Recommendations
    critical_recommendations: List[str] = field(default_factory=list)

class ReportGenerator:
    """
    Generates comprehensive reports from pipeline execution results.
    Supports JSON, HTML, and plain text formats.
    """
    
    def __init__(self, config: PipelineConfig, output_dir: str = "./pipeline_output"):
        self.config = config
        self.output_dir = Path(output_dir)
        self.logger = get_logger(__name__)
        
        # Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Template paths
        self.log_dir = Path("./logs")

    def compile_all_stages(self, 
                          stage_1_json: Dict[str, Any],
                          stage_2_json: Optional[Dict[str, Any]] = None,
                          stage_3_json: Optional[Dict[str, Any]] = None,
                          stage_4_json: Optional[Dict[str, Any]] = None,
                          stage_5_scaleup: Optional[Dict[str, Any]] = None,
                          stage_5_downstream: Optional[Dict[str, Any]] = None,
                          stage_5_regulatory: Optional[Dict[str, Any]] = None,
                          errors: Optional[List[str]] = None,
                          warnings: Optional[List[str]] = None) -> PipelineReport:
        """
        Compile all stage outputs into a unified report.
        """
        self.logger.info("=== COMPILING ALL STAGE RESULTS ===")
        
        # Extract basic info
        pipeline_id = stage_1_json.get('pipeline_id', 'unknown')
        organism = stage_1_json.get('organism', {}).get('name', 'unknown')
        molecule = stage_1_json.get('target_molecule', {}).get('name', 'unknown')
        
        # Determine overall status
        statuses = []
        if stage_1_json.get('stage_1_status') == 'PASS':
            statuses.append('S1_OK')
        elif stage_1_json.get('stage_1_status') == 'WARN':
            statuses.append('S1_WARN')
        else:
            statuses.append('S1_FAIL')
            
        if stage_2_json:
            if stage_2_json.get('stage_2_status') == 'PASS':
                statuses.append('S2_OK')
            elif stage_2_json.get('stage_2_status') == 'WARN':
                statuses.append('S2_WARN')
            else:
                statuses.append('S2_FAIL')
                
        if stage_3_json:
            if stage_3_json.get('stage_3_status') == 'PASS':
                statuses.append('S3_OK')
            elif stage_3_json.get('stage_3_status') == 'WARN':
                statuses.append('S3_WARN')
            else:
                statuses.append('S3_FAIL')
                
        if stage_4_json:
            if stage_4_json.get('stage_4_status') == 'PASS':
                statuses.append('S4_OK')
            elif stage_4_json.get('stage_4_status') == 'WARN':
                statuses.append('S4_WARN')
            else:
                statuses.append('S4_FAIL')
        
        # Determine overall status
        if 'FAIL' in str(statuses):
            overall_status = 'PARTIAL'
        elif 'WARN' in str(statuses):
            overall_status = 'SUCCESS_WITH_WARNINGS'
        else:
            overall_status = 'SUCCESS'
        
        # Extract key metrics
        final_titer = 0.0
        production_titer = 0.0
        ds_cost = 0.0
        reg_feas = "UNKNOWN"
        
        if stage_4_json and 'fermentation_simulation' in stage_4_json:
            final_titer = stage_4_json['fermentation_simulation'].get('final_titer_g_per_l', 0.0)
            
        if stage_5_scaleup:
            production_titer = stage_5_scaleup.get('final_production_titer_g_per_l', final_titer)
            
        if stage_5_downstream:
            ds_cost = stage_5_downstream.get('estimated_cost_per_kg_usd', 0.0)
            
        if stage_5_regulatory:
            reg_feas = stage_5_regulatory.get('compliance_summary', {}).get('regulatory_feasibility', 'UNKNOWN')
        
        # Calculate feasibility score (0-100)
        feasibility_score = self._calculate_feasibility_score(
            titer=production_titer,
            cost=ds_cost,
            regulatory=reg_feas,
            pathway_success=len(stage_2_json.get('pathway_candidates', [])) if stage_2_json else 0
        )
        
        # Collect recommendations
        recommendations = []
        if stage_5_regulatory:
            recommendations.extend(stage_5_regulatory.get('compliance_summary', {}).get('critical_recommendations', []))
        if stage_5_scaleup:
            for cascade in stage_5_scaleup.get('scale_up_cascade', []):
                recommendations.extend(cascade.get('recommendations', []))
        
        report = PipelineReport(
            pipeline_id=pipeline_id,
            timestamp=datetime.now().isoformat(),
            organism=organism,
            target_molecule=molecule,
            status=overall_status,
            stage_1_summary=stage_1_json,
            stage_2_summary=stage_2_json,
            stage_3_summary=stage_3_json,
            stage_4_summary=stage_4_json,
            stage_5_summary={
                'scaleup': stage_5_scaleup,
                'downstream': stage_5_downstream,
                'regulatory': stage_5_regulatory
            },
            final_predicted_titer_g_per_l=final_titer,
            final_production_titer_g_per_l=production_titer,
            downstream_cost_per_kg_usd=ds_cost,
            regulatory_feasibility=reg_feas,
            overall_feasibility_score=feasibility_score,
            errors=errors or [],
            warnings=warnings or [],
            critical_recommendations=recommendations
        )
        
        self.logger.info(f"Report compiled. Status: {overall_status}, Feasibility: {feasibility_score:.1f}/100")
        return report

    def _calculate_feasibility_score(self, titer: float, cost: float, 
                                     regulatory: str, pathway_success: int) -> float:
        """Calculate overall project feasibility score."""
        score = 0.0
        
        # Titer contribution (max 40 points)
        if titer > 50:
            score += 40
        elif titer > 20:
            score += 30
        elif titer > 10:
            score += 20
        elif titer > 5:
            score += 10
        elif titer > 1:
            score += 5
            
        # Cost contribution (max 30 points)
        if cost < 50:
            score += 30
        elif cost < 100:
            score += 25
        elif cost < 200:
            score += 20
        elif cost < 500:
            score += 10
            
        # Regulatory contribution (max 20 points)
        if regulatory == 'HIGH':
            score += 20
        elif regulatory == 'MEDIUM':
            score += 10
            
        # Pathway success contribution (max 10 points)
        if pathway_success > 3:
            score += 10
        elif pathway_success > 1:
            score += 5
            
        return min(100.0, score)

    def generate_summary_report(self, report: PipelineReport) -> str:
        """Generate human-readable text summary."""
        lines = [
            "=" * 80,
            "METABOLIC PATHWAY PIPELINE - FINAL REPORT",
            "=" * 80,
            "",
            f"Pipeline ID: {report.pipeline_id}",
            f"Timestamp: {report.timestamp}",
            f"Organism: {report.organism}",
            f"Target Molecule: {report.target_molecule}",
            f"Overall Status: {report.status}",
            "",
            "-" * 80,
            "KEY METRICS",
            "-" * 80,
            f"Lab-scale Predicted Titer: {report.final_predicted_titer_g_per_l:.2f} g/L",
            f"Production-scale Predicted Titer: {report.final_production_titer_g_per_l:.2f} g/L",
            f"Downstream Cost Estimate: ${report.downstream_cost_per_kg_usd:.2f}/kg",
            f"Regulatory Feasibility: {report.regulatory_feasibility}",
            f"Overall Feasibility Score: {report.overall_feasibility_score:.1f}/100",
            ""
        ]
        
        # Stage summaries
        lines.append("-" * 80)
        lines.append("STAGE SUMMARIES")
        lines.append("-" * 80)
        
        if report.stage_1_summary:
            s1 = report.stage_1_summary
            lines.append(f"Stage 1 (Data Layer): {s1.get('stage_1_status', 'UNKNOWN')}")
            lines.append(f"  - Completeness Score: {s1.get('data_quality_report', {}).get('completeness_score', 'N/A')}")
            
        if report.stage_2_summary:
            s2 = report.stage_2_summary
            lines.append(f"Stage 2 (Pathway Prediction): {s2.get('stage_2_status', 'UNKNOWN')}")
            pathways = s2.get('pathway_candidates', [])
            lines.append(f"  - Pathways Found: {len(pathways)}")
            if pathways:
                lines.append(f"  - Best Pathway Yield: {pathways[0].get('predicted_yield_mol_per_mol', 0):.3f} mol/mol")
                
        if report.stage_3_summary:
            s3 = report.stage_3_summary
            lines.append(f"Stage 3 (Flux Analysis): {s3.get('stage_3_status', 'UNKNOWN')}")
            strain = s3.get('strain_design', {})
            lines.append(f"  - Predicted Productivity: {strain.get('predicted_productivity_g_per_l_per_h', 0):.3f} g/L/h")
            lines.append(f"  - Knockouts Recommended: {len(strain.get('final_knockouts', []))}")
            
        if report.stage_4_summary:
            s4 = report.stage_4_summary
            lines.append(f"Stage 4 (DBTL + Fermentation): {s4.get('stage_4_status', 'UNKNOWN')}")
            ferm = s4.get('fermentation_simulation', {})
            lines.append(f"  - Final Titer: {ferm.get('final_titer_g_per_l', 0):.2f} g/L")
            lines.append(f"  - DBTL Cycles Simulated: {len(s4.get('dbtl_cycles', []))}")
            
        if report.stage_5_summary:
            lines.append("Stage 5 (Scale-up + Downstream + Regulatory): COMPLETE")
            if report.stage_5_summary.get('scaleup'):
                su = report.stage_5_summary['scaleup']
                lines.append(f"  - Scale-up Risk: {su.get('scale_up_risk', 'N/A')}")
            if report.stage_5_summary.get('downstream'):
                ds = report.stage_5_summary['downstream']
                lines.append(f"  - Overall Recovery: {ds.get('overall_recovery_percent', 0):.1f}%")
                lines.append(f"  - Economic Feasibility: {ds.get('economic_feasibility', 'N/A')}")
            if report.stage_5_summary.get('regulatory'):
                reg = report.stage_5_summary['regulatory']
                lines.append(f"  - Compliance Score: {reg.get('compliance_summary', {}).get('overall_score', 0)}/100")
                
        lines.append("")
        
        # Warnings and Errors
        if report.warnings:
            lines.append("-" * 80)
            lines.append("WARNINGS")
            lines.append("-" * 80)
            for i, warn in enumerate(report.warnings[:10], 1):  # Limit to 10
                lines.append(f"{i}. {warn}")
            if len(report.warnings) > 10:
                lines.append(f"... and {len(report.warnings) - 10} more warnings")
            lines.append("")
            
        if report.errors:
            lines.append("-" * 80)
            lines.append("ERRORS")
            lines.append("-" * 80)
            for i, err in enumerate(report.errors[:10], 1):
                lines.append(f"{i}. {err}")
            if len(report.errors) > 10:
                lines.append(f"... and {len(report.errors) - 10} more errors")
            lines.append("")
            
        # Recommendations
        if report.critical_recommendations:
            lines.append("-" * 80)
            lines.append("CRITICAL RECOMMENDATIONS")
            lines.append("-" * 80)
            for i, rec in enumerate(report.critical_recommendations, 1):
                lines.append(f"{i}. {rec}")
            lines.append("")
            
        lines.append("=" * 80)
        lines.append("END OF REPORT")
        lines.append("=" * 80)
        
        return "\n".join(lines)

    def generate_json_report(self, report: PipelineReport) -> Dict[str, Any]:
        """Generate machine-readable JSON report."""
        # Convert dataclass to dict
        report_dict = {
            'pipeline_id': report.pipeline_id,
            'timestamp': report.timestamp,
            'organism': report.organism,
            'target_molecule': report.target_molecule,
            'status': report.status,
            'key_metrics': {
                'lab_scale_titer_g_per_l': report.final_predicted_titer_g_per_l,
                'production_scale_titer_g_per_l': report.final_production_titer_g_per_l,
                'downstream_cost_per_kg_usd': report.downstream_cost_per_kg_usd,
                'regulatory_feasibility': report.regulatory_feasibility,
                'overall_feasibility_score': report.overall_feasibility_score
            },
            'stage_outputs': {
                'stage_1': report.stage_1_summary,
                'stage_2': report.stage_2_summary,
                'stage_3': report.stage_3_summary,
                'stage_4': report.stage_4_summary,
                'stage_5': report.stage_5_summary
            },
            'quality_indicators': {
                'errors': report.errors,
                'warnings': report.warnings,
                'critical_recommendations': report.critical_recommendations
            }
        }
        
        return report_dict

    def generate_error_report(self, errors: List[str], warnings: List[str]) -> str:
        """Generate dedicated error/warning report."""
        lines = [
            "PIPELINE ERROR AND WARNING REPORT",
            "=" * 60,
            f"Generated: {datetime.now().isoformat()}",
            "",
            f"Total Errors: {len(errors)}",
            f"Total Warnings: {len(warnings)}",
            ""
        ]
        
        if errors:
            lines.append("ERRORS:")
            lines.append("-" * 60)
            for i, err in enumerate(errors, 1):
                lines.append(f"[{i}] {err}")
            lines.append("")
            
        if warnings:
            lines.append("WARNINGS:")
            lines.append("-" * 60)
            for i, warn in enumerate(warnings, 1):
                lines.append(f"[{i}] {warn}")
                
        return "\n".join(lines)

    def export_to_files(self, report: PipelineReport, 
                       errors: Optional[List[str]] = None,
                       warnings: Optional[List[str]] = None) -> Dict[str, str]:
        """Export report to multiple file formats."""
        self.logger.info(f"Exporting reports to {self.output_dir}")
        
        saved_files = {}
        
        # 1. JSON Report
        json_path = self.output_dir / "final_report.json"
        json_data = self.generate_json_report(report)
        with open(json_path, 'w') as f:
            json.dump(json_data, f, indent=2, default=str)
        saved_files['json'] = str(json_path)
        self.logger.info(f"Saved JSON report: {json_path}")
        
        # 2. Text Summary
        txt_path = self.output_dir / "summary_report.txt"
        txt_content = self.generate_summary_report(report)
        with open(txt_path, 'w') as f:
            f.write(txt_content)
        saved_files['text'] = str(txt_path)
        self.logger.info(f"Saved text report: {txt_path}")
        
        # 3. Error Report (if any)
        if errors or warnings:
            err_path = self.output_dir / "error_report.txt"
            err_content = self.generate_error_report(errors or [], warnings or [])
            with open(err_path, 'w') as f:
                f.write(err_content)
            saved_files['errors'] = str(err_path)
            self.logger.info(f"Saved error report: {err_path}")
            
        # 4. HTML Report (simple version)
        html_path = self.output_dir / "report.html"
        html_content = self._generate_html_report(report)
        with open(html_path, 'w') as f:
            f.write(html_content)
        saved_files['html'] = str(html_path)
        self.logger.info(f"Saved HTML report: {html_path}")
        
        return saved_files

    def _generate_html_report(self, report: PipelineReport) -> str:
        """Generate simple HTML report."""
        html = f"""<!DOCTYPE html>
<html>
<head>
    <title>Pipeline Report - {report.pipeline_id}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        h1 {{ color: #2c3e50; }}
        h2 {{ color: #34495e; border-bottom: 2px solid #ecf0f1; }}
        .metric {{ background: #ecf0f1; padding: 10px; margin: 5px 0; border-radius: 5px; }}
        .success {{ color: #27ae60; }}
        .warning {{ color: #f39c12; }}
        .error {{ color: #e74c3c; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #3498db; color: white; }}
    </style>
</head>
<body>
    <h1>Metabolic Pathway Pipeline Report</h1>
    <p><strong>Pipeline ID:</strong> {report.pipeline_id}</p>
    <p><strong>Timestamp:</strong> {report.timestamp}</p>
    <p><strong>Status:</strong> <span class="{'success' if 'SUCCESS' in report.status else 'warning'}">{report.status}</span></p>
    
    <h2>Project Overview</h2>
    <table>
        <tr><th>Organism</th><td>{report.organism}</td></tr>
        <tr><th>Target Molecule</th><td>{report.target_molecule}</td></tr>
        <tr><th>Feasibility Score</th><td>{report.overall_feasibility_score:.1f}/100</td></tr>
    </table>
    
    <h2>Key Metrics</h2>
    <div class="metric"><strong>Lab-scale Titer:</strong> {report.final_predicted_titer_g_per_l:.2f} g/L</div>
    <div class="metric"><strong>Production Titer:</strong> {report.final_production_titer_g_per_l:.2f} g/L</div>
    <div class="metric"><strong>Downstream Cost:</strong> ${report.downstream_cost_per_kg_usd:.2f}/kg</div>
    <div class="metric"><strong>Regulatory Feasibility:</strong> {report.regulatory_feasibility}</div>
    
    <h2>Recommendations</h2>
    <ul>
"""
        for rec in report.critical_recommendations[:10]:
            html += f"        <li>{rec}</li>\n"
            
        html += """    </ul>
    
    <p><em>Full details available in JSON report.</em></p>
</body>
</html>"""
        
        return html


def run_stage_5_reporting(stage_1_json: Dict[str, Any],
                         stage_2_json: Optional[Dict[str, Any]],
                         stage_3_json: Optional[Dict[str, Any]],
                         stage_4_json: Optional[Dict[str, Any]],
                         stage_5_scaleup: Optional[Dict[str, Any]],
                         stage_5_downstream: Optional[Dict[str, Any]],
                         stage_5_regulatory: Optional[Dict[str, Any]],
                         errors: Optional[List[str]],
                         warnings: Optional[List[str]],
                         output_dir: str = "./pipeline_output") -> Dict[str, str]:
    """
    Master function to generate all Stage 5 reports.
    Returns paths to generated files.
    """
    logger = get_logger(__name__)
    logger.info("=== STAGE 5 REPORT GENERATION START ===")
    
    try:
        config = PipelineConfig()
        generator = ReportGenerator(config, output_dir)
        
        # Compile report
        report = generator.compile_all_stages(
            stage_1_json=stage_1_json,
            stage_2_json=stage_2_json,
            stage_3_json=stage_3_json,
            stage_4_json=stage_4_json,
            stage_5_scaleup=stage_5_scaleup,
            stage_5_downstream=stage_5_downstream,
            stage_5_regulatory=stage_5_regulatory,
            errors=errors,
            warnings=warnings
        )
        
        # Export files
        saved_files = generator.export_to_files(report, errors, warnings)
        
        logger.info(f"Reports exported successfully: {list(saved_files.keys())}")
        return saved_files
        
    except Exception as e:
        logger.error(f"Report generation failed: {e}", exc_info=True)
        raise PipelineError(f"Stage 5 Reporting failed: {str(e)}")
