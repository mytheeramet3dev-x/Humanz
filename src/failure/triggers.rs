//! Trigger conditions and calibrated clinical/biological safety thresholds.
//!
//! All grading bands `[Mild, Moderate, Severe, Critical]` are anchored to established clinical
//! or experimental reference standards (ADA, KDIGO, CTCAE, Gill et al. 2022, Ocampo et al. 2016).

use super::adverse::{grade, AdverseEvent, AdverseEventKind, Severity};
use scies_bio_th::aging::{InVivoReprogrammingReport, OskPulseRegimen, SafetyVerdict};
use scies_bio_th::anatomy::{AnatomicalSpatialAtlas, OrganArchitecture};

/// Somatic dedifferentiation / teratoma hazard bands (`[0.10, 0.25, 0.50, 0.80]`), anchored to
/// enhancer bookmark retention (`<80%` Mild/Moderate caution, `<70%` Severe point-of-no-return;
/// Ocampo et al. 2016, Gill et al. 2022).
pub const TERATOMA_HAZARD_BANDS: [f64; 4] = [0.10, 0.25, 0.50, 0.80];

/// Composite oncogenic risk bands (`[0.10, 0.15, 0.30, 0.50]`), anchored to
/// `MAX_SAFE_ONCOGENIC_RISK_THRESHOLD = 0.15` as the Moderate lockout boundary.
pub const ONCOGENIC_RISK_BANDS: [f64; 4] = [0.10, 0.15, 0.30, 0.50];

/// Cumulative vector exposure ratio bands (`[1.8, 2.6, 3.6, 5.0]` x baseline 6.3 dose-days)
/// for anti-vector innate/adaptive immune activation.
pub const VECTOR_IMMUNE_EXPOSURE_BANDS: [f64; 4] = [1.8, 2.6, 3.6, 5.0];

/// Organ PBPK parenchymal exposure ratio bands (`[1.5, 2.0, 2.6, 3.4]` x reference bioavailability
/// 0.56), corresponding to CTCAE Grade 1–4 organ-specific toxicity.
pub const ORGAN_PBPK_OVERLOAD_BANDS: [f64; 4] = [1.5, 2.0, 2.6, 3.4];

/// Inter-organ functional spread bands (`[0.15, 0.30, 0.50, 0.70]`) across heart/kidney/brain.
pub const SYSTEMIC_IMBALANCE_SPREAD_BANDS: [f64; 4] = [0.15, 0.30, 0.50, 0.70];

/// ADA 2-hour 75g OGTT post-load glucose bands (`mg/dL`):
/// `<140` normal, `>=140` Mild IGT, `>=165` Moderate IGT, `>=200` Severe diabetic stress hyperglycemia.
pub const OGTT_2H_GLUCOSE_BANDS_MG_DL: [f64; 4] = [140.0, 165.0, 200.0, 260.0];

/// KDIGO eGFR deficit below `60.0 mL/min/1.73m^2` (`G3a < 60`, `G3b < 45`, `G4 < 30`, `G5 < 15`).
pub const KDIGO_EGFR_DEFICIT_BANDS_ML_MIN: [f64; 4] = [0.1, 15.0, 30.0, 45.0];

/// Core temperature deficit below `36.3 °C` (`>=0.3 °C` Mild `<36.0 °C`, `>=0.8 °C` Moderate `<35.5 °C`,
/// `>=1.5 °C` Severe `<34.8 °C` hypothermia).
pub const HYPOTHERMIA_DEFICIT_BANDS_C: [f64; 4] = [0.3, 0.8, 1.5, 2.5];

/// Blood-Brain Barrier TEER deficit below `1100.0 Ohm*cm^2` neurovascular integrity threshold.
pub const BBB_TEER_DEFICIT_BANDS_OHM_CM2: [f64; 4] = [30.0, 110.0, 220.0, 380.0];

/// Post-withdrawal aging rate ratio (`treated_rate / control_rate`) bands for senescence rebound.
pub const SENESCENCE_REBOUND_RATIO_BANDS: [f64; 4] = [1.15, 1.35, 1.65, 2.20];

/// Cumulative vector exposure of a regimen (baseline balanced regimen ≈ 6.3).
pub fn cumulative_vector_exposure(r: &OskPulseRegimen) -> f64 {
    r.total_cycles as f64 * r.factor_dosage * (r.pulse_duration_days / 2.0)
}

/// Relative per-cycle dose/expression load (balanced regimen = 1.0).
pub fn relative_dose_load(r: &OskPulseRegimen) -> f64 {
    (r.factor_dosage / 1.05) * (r.pulse_duration_days / 2.0)
}

/// Detects teratoma, clonal-expansion, immune and organ-toxicity events for a regimen.
///
/// Organ toxicity uses each organ's PDE parenchymal bioavailability relative to a 0.56
/// reference, so organs that take up the vector readily (liver, marrow) are exposed most.
pub fn detect_regimen_adverse_events(
    report: &InVivoReprogrammingReport,
    regimen: &OskPulseRegimen,
    atlas: &AnatomicalSpatialAtlas,
) -> Vec<AdverseEvent> {
    let mut ev = Vec::new();

    // 1. Teratoma / pluripotency risk
    let mut sev = grade(report.teratoma_hazard_index, TERATOMA_HAZARD_BANDS);
    match report.safety_verdict {
        SafetyVerdict::CarcinogenicMycRisk => sev = sev.max(Some(Severity::Critical)),
        SafetyVerdict::SevereTeratomaRiskExceededDedifferentiationThreshold => {
            sev = sev.max(Some(Severity::Severe))
        }
        _ => {}
    }
    if let Some(s) = sev {
        ev.push(AdverseEvent {
            kind: AdverseEventKind::TeratomaRisk,
            organ: None,
            severity: s,
            reversible: s < Severity::Severe,
            magnitude: report.teratoma_hazard_index,
            description: format!(
                "Teratoma hazard index {:.2}, somatic identity {:.0}%",
                report.teratoma_hazard_index,
                report.final_somatic_identity_retention * 100.0
            ),
        });
    }

    // 2. Clonal expansion / pre-cancerous growth
    let mut sev = grade(report.composite_oncogenic_risk, ONCOGENIC_RISK_BANDS);
    if report.safety_verdict == SafetyVerdict::OncogenicTelomeraseThresholdExceeded {
        sev = sev.max(Some(Severity::Moderate));
    }
    if let Some(s) = sev {
        ev.push(AdverseEvent {
            kind: AdverseEventKind::ClonalExpansion,
            organ: None,
            severity: s,
            reversible: false,
            magnitude: report.composite_oncogenic_risk,
            description: format!(
                "Composite oncogenic risk {:.3} (TERT x{:.1}, checkpoint x{:.2})",
                report.composite_oncogenic_risk,
                regimen.transient_tert_mrna_fold,
                regimen.tp53_cdkn2a_safety_fold
            ),
        });
    }

    // 3. Immune reaction against vector / factors (cumulative exposure)
    let exposure_ratio = cumulative_vector_exposure(regimen) / 6.3;
    if let Some(s) = grade(exposure_ratio, VECTOR_IMMUNE_EXPOSURE_BANDS) {
        ev.push(AdverseEvent {
            kind: AdverseEventKind::ImmuneReactionToVector,
            organ: None,
            severity: s,
            reversible: true,
            magnitude: exposure_ratio,
            description: format!(
                "Cumulative vector exposure {:.1}x baseline; anti-vector immunity lowers bioavailability",
                exposure_ratio
            ),
        });
    }

    // 4. Organ-specific toxicity
    let load = relative_dose_load(regimen);
    let organs = [
        (2usize, AdverseEventKind::HepaticToxicity, "Liver"),
        (6, AdverseEventKind::CardiacFibrosis, "Heart"),
        (1, AdverseEventKind::Neuroinflammation, "Brain"),
        (7, AdverseEventKind::RenalToxicity, "Kidney"),
    ];
    for (idx, kind, name) in organs {
        let bioavail = atlas.organs[idx].mean_parenchymal_bioavailability(
            scies_bio_th::aging::MW_MRNA_LNP_VECTOR_DA,
            0.92,
        );
        let ratio = load * (bioavail / 0.56);
        if let Some(s) = grade(ratio, ORGAN_PBPK_OVERLOAD_BANDS) {
            ev.push(AdverseEvent {
                kind,
                organ: Some(name.to_string()),
                severity: s,
                reversible: s < Severity::Critical,
                magnitude: ratio,
                description: format!(
                    "{} exposure {:.2}x reference (bioavailability {:.0}%)",
                    name,
                    ratio,
                    bioavail * 100.0
                ),
            });
        }
    }
    ev
}
