//! Whole-Body Safety Guardrails enforcing Zero-Teratoma, Zero-Oncogene, and Zero-Hyperplasia
//! invariants across all 6 organ compartments in `HumanZ`.

use crate::body::{HumanDigitalTwin, OrganSystem};
use serde::{Deserialize, Serialize};

/// Comprehensive whole-body safety audit across all 6 organs.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WholeBodySafetyAudit {
    pub all_organs_approved: bool,
    pub max_oncogenic_risk: f64,
    pub min_somatic_identity_retention: f64,
    pub max_stem_cell_fraction: f64,
    pub max_pluripotency: f64,
    pub organ_verdicts: Vec<(OrganSystem, bool)>,
}

/// Evaluates strict multi-organ safety guardrails across the `HumanDigitalTwin`.
pub fn audit_whole_body_safety(
    twin: &HumanDigitalTwin,
    max_oncogenic_risk_observed: f64,
    min_identity_observed: f64,
) -> WholeBodySafetyAudit {
    let snap = twin.snapshot_whole_body(true);
    let mut organ_verdicts = Vec::new();
    let mut all_ok = true;

    for (&organ, &col_id) in &twin.treated_organs {
        let rep = twin.sim.audit_colony_safety(col_id).unwrap();
        let organ_stat = &snap.organs[&organ].stats;
        let stem_safe =
            organ_stat.stem_cell_fraction <= scies_bio_th::aging::MAX_SAFE_STEM_CELL_FRACTION + 1e-9;
        let pluri_safe = organ_stat.mean_pluripotency < 0.30;
        let passed = rep.approved && stem_safe && pluri_safe;
        if !passed {
            all_ok = false;
        }
        organ_verdicts.push((organ, passed));
    }

    if max_oncogenic_risk_observed >= scies_bio_th::aging::MAX_SAFE_ONCOGENIC_RISK_THRESHOLD
        || min_identity_observed < 0.90
    {
        all_ok = false;
    }

    WholeBodySafetyAudit {
        all_organs_approved: all_ok,
        max_oncogenic_risk: max_oncogenic_risk_observed,
        min_somatic_identity_retention: min_identity_observed,
        max_stem_cell_fraction: snap.max_stem_cell_fraction,
        max_pluripotency: snap.max_pluripotency,
        organ_verdicts,
    }
}
