//! # Telemetry Export & Persistence (`export`)
//!
//! Exports full `HumanZ` simulation bundles to structured JSON and CSV files for downstream
//! plotting, clinical auditing, and reproducible archival.

use crate::failure::FailureCampaignReport;
use crate::protocol::{NegligibleSenescenceRunSummary, PersonalizedComorbidityRescueReport};
use crate::stress_test::ClinicalStressChallengeReport;
use crate::validation::{ClinicalReferenceValidationReport, MonteCarloCohortReport};
use scies_bio_sim::prelude::*;
use serde::{Deserialize, Serialize};
use std::fs;
use std::io::Write;
use std::path::Path;

/// Complete serializable `HumanZ` simulation bundle combining all subsystems.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HumanzCompleteSimulationBundle {
    pub summary: NegligibleSenescenceRunSummary,
    pub stress_challenge: ClinicalStressChallengeReport,
    pub reference_validation: ClinicalReferenceValidationReport,
    pub monte_carlo_cohort: MonteCarloCohortReport,
    pub failure_campaign: FailureCampaignReport,
    pub personalized_rescue: PersonalizedComorbidityRescueReport,
}

/// Writes the complete `HumanzCompleteSimulationBundle` as pretty-printed JSON to `path`.
pub fn export_simulation_bundle_json(
    bundle: &HumanzCompleteSimulationBundle,
    path: &Path,
) -> Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .map_err(|e| BioSimError::Other(format!("mkdir error: {}", e)))?;
    }
    let json = serde_json::to_string_pretty(bundle)
        .map_err(|e| BioSimError::Other(format!("JSON serialization error: {}", e)))?;
    fs::write(path, json).map_err(|e| BioSimError::Other(format!("JSON write error: {}", e)))
}

/// Exports the year-by-year longitudinal biomarker trajectory (`control` vs `treated`) to CSV.
pub fn export_longitudinal_trajectory_csv(
    summary: &NegligibleSenescenceRunSummary,
    path: &Path,
) -> Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .map_err(|e| BioSimError::Other(format!("mkdir error: {}", e)))?;
    }
    let mut file = fs::File::create(path)
        .map_err(|e| BioSimError::Other(format!("CSV create error: {}", e)))?;
    writeln!(
        file,
        "year_index,chronological_age,control_bio_age,treated_bio_age,treated_telomere_bp,treated_viability,treated_egfr,treated_sbp_mmhg,treated_glymphatic_ul_min"
    )
    .map_err(|e| BioSimError::Other(format!("CSV write error: {}", e)))?;

    let n = summary
        .longitudinal_treated_history
        .len()
        .min(summary.longitudinal_control_history.len());
    for i in 0..n {
        let t = &summary.longitudinal_treated_history[i];
        let c = &summary.longitudinal_control_history[i];
        writeln!(
            file,
            "{},{:.1},{:.3},{:.3},{:.1},{:.4},{:.2},{:.2},{:.2}",
            i,
            t.chronological_age_years,
            c.mean_biological_age_years,
            t.mean_biological_age_years,
            t.mean_telomere_bp,
            t.mean_viability,
            t.circuits.renal_raas.egfr_ml_min,
            t.circuits.cardiovascular.systolic_bp_mmhg,
            t.circuits.neuroendocrine_glymphatic.cortical_glymphatic_clearance_rate,
        )
        .map_err(|e| BioSimError::Other(format!("CSV row error: {}", e)))?;
    }
    Ok(())
}

/// Exports the 8-scenario Failure Modes & Adverse Events campaign summary to CSV.
pub fn export_failure_campaign_csv(campaign: &FailureCampaignReport, path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .map_err(|e| BioSimError::Other(format!("mkdir error: {}", e)))?;
    }
    let mut file = fs::File::create(path)
        .map_err(|e| BioSimError::Other(format!("CSV create error: {}", e)))?;
    writeln!(
        file,
        "scenario_name,outcome,guard_action,efficacy_fraction,safety_activated,hazard_contained,events_count"
    )
    .map_err(|e| BioSimError::Other(format!("CSV header error: {}", e)))?;

    for r in &campaign.results {
        writeln!(
            file,
            "\"{}\",{:?},\"{:?}\",{:.4},{},{},{}",
            r.name,
            r.outcome,
            r.decision.action,
            r.efficacy_fraction,
            r.safety_activated,
            r.decision.hazard_contained(),
            r.events.len()
        )
        .map_err(|e| BioSimError::Other(format!("CSV row error: {}", e)))?;
    }
    Ok(())
}
