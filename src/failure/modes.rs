//! Failure-mode definitions (`FailureMode`) and outcome classification (`OutcomeClass`).

use scies_bio_th::aging::OskPulseRegimen;
use serde::{Deserialize, Serialize};

/// The 4 failure groups (A–D).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FailureGroup {
    /// A: over-intervention / overdose.
    OverIntervention,
    /// B: incomplete or partial response.
    IncompleteResponse,
    /// C: safety-system failure.
    SafetySystemFailure,
    /// D: external / physiological stress amplification and pre-existing damage.
    ExternalStress,
}

/// A concrete failure mode that can be injected into a subject / regimen.
///
/// `intensity` (see [`inject_failure_mode`]) scales the deviation from the balanced baseline:
/// `0.0` = baseline, `1.0` = the nominal severity encoded in the variant.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum FailureMode {
    /// Reference arm: the balanced protocol, no injected failure (must stay a Success).
    BalancedBaseline,
    /// Group A: factor dose and expression persistence multiplied by `dose_multiplier`.
    OskOverdose { dose_multiplier: f64 },
    /// Group A/C: strong TERT pulse with the tumor-suppressor checkpoint weakened (`tp53_fold < 1`).
    TelomeraseWithoutCheckpoint { tert_fold: f64, tp53_fold: f64 },
    /// Group C: same oncogenic insult, but the suicide switch (iCasp9) does not fire.
    SuicideSwitchFailure { tert_fold: f64, tp53_fold: f64 },
    /// Group B: only a subset of organs receives rejuvenation.
    PartialOrganResponse,
    /// Group D: hypotension + endotoxemia + cold + exercise applied while repair is incomplete.
    CombinedSevereStress { hypotension_map_mmhg: f64 },
    /// Group D: therapy starts at very advanced age / damage.
    LateStageHighDamage { initial_age_years: f64 },
    /// Group B/D: therapy stops; checks for senescence rebound.
    TreatmentWithdrawal { years_off_therapy: u32 },
}

impl FailureMode {
    /// Failure group of this mode.
    pub fn group(&self) -> FailureGroup {
        match self {
            FailureMode::BalancedBaseline => FailureGroup::ExternalStress,
            FailureMode::OskOverdose { .. } => FailureGroup::OverIntervention,
            FailureMode::TelomeraseWithoutCheckpoint { .. } => FailureGroup::OverIntervention,
            FailureMode::SuicideSwitchFailure { .. } => FailureGroup::SafetySystemFailure,
            FailureMode::PartialOrganResponse => FailureGroup::IncompleteResponse,
            FailureMode::CombinedSevereStress { .. } => FailureGroup::ExternalStress,
            FailureMode::LateStageHighDamage { .. } => FailureGroup::ExternalStress,
            FailureMode::TreatmentWithdrawal { .. } => FailureGroup::IncompleteResponse,
        }
    }
}

/// Standard balanced regimen (identical to the one used by the main protocol).
pub fn balanced_regimen(lineage: &str) -> OskPulseRegimen {
    OskPulseRegimen {
        somatic_lineage: lineage.to_string(),
        delivery_vector: scies_bio_th::aging::DeliveryVector::MrnaLipidNanoparticle,
        pulse_duration_days: 2.0,
        washout_interval_days: 5.0,
        total_cycles: 6,
        tet1_tet2_competence: 0.95,
        oncogene_cmyc_included: false,
        factor_dosage: 1.05,
        transient_tert_mrna_fold: 2.25,
        tp53_cdkn2a_safety_fold: 1.85,
    }
}

/// Injects `mode` into a regimen at the given `intensity` (regimen-level modes only;
/// organ-/subject-level modes are handled by the scenario runners).
pub fn inject_failure_mode(regimen: &mut OskPulseRegimen, mode: &FailureMode, intensity: f64) {
    let k = intensity.max(0.0);
    match mode {
        FailureMode::OskOverdose { dose_multiplier } => {
            let m = 1.0 + (dose_multiplier - 1.0) * k;
            // Higher mRNA dose also prolongs effective expression time.
            regimen.factor_dosage *= m;
            regimen.pulse_duration_days *= m;
        }
        FailureMode::TelomeraseWithoutCheckpoint { tert_fold, tp53_fold }
        | FailureMode::SuicideSwitchFailure { tert_fold, tp53_fold } => {
            regimen.transient_tert_mrna_fold = 1.0 + (tert_fold - 1.0) * k;
            regimen.tp53_cdkn2a_safety_fold = (1.0 + (tp53_fold - 1.0) * k).max(0.05);
        }
        _ => {}
    }
}

/// Final classification of a scenario.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum OutcomeClass {
    Success,
    Partial,
    Failure,
    Catastrophic,
}
