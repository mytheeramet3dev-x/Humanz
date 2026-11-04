//! # Failure Modes & Adverse Events (`failure`)
//!
//! Verifies that the Safety Guard actually reacts, that failure is possible and classified, and
//! that the balanced protocol is not harmed by the extra safeguards.
//!
//! - [`modes`]: `FailureMode`, `OutcomeClass`, `inject_failure_mode`
//! - [`adverse`]: `AdverseEvent`, `Severity`
//! - [`triggers`]: telemetry → adverse events
//! - [`safety_response`]: `SafetyAction` policy + decision log
//! - [`scenarios`]: runnable scenarios

pub mod adverse;
pub mod modes;
pub mod safety_response;
pub mod scenarios;
pub mod triggers;

use adverse::Severity;
use modes::{FailureMode, OutcomeClass};
use safety_response::SafetyAction;
use scenarios::FailureScenarioResult;
use scies_bio_sim::prelude::*;
use serde::{Deserialize, Serialize};

/// Aggregated campaign metrics.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FailureCampaignReport {
    pub results: Vec<FailureScenarioResult>,
    pub success: usize,
    pub partial: usize,
    pub failure: usize,
    pub catastrophic: usize,
    /// Fraction of scenarios not classified Success.
    pub failure_rate: f64,
    /// Among hazardous scenarios, fraction where the guard did not just `Continue`.
    pub safety_activation_rate: f64,
    /// Among hazardous scenarios, fraction where the hazard was actually contained.
    pub containment_rate: f64,
    /// Mean cycles until the guard saw a leave-safe-window (regimen scenarios only).
    pub mean_time_to_safety_activation: Option<f64>,
    /// Event counts `[Mild, Moderate, Severe, Critical]`.
    pub severity_histogram: [usize; 4],
}

/// Builds the aggregate report from scenario results.
pub fn summarize_campaign(results: Vec<FailureScenarioResult>) -> FailureCampaignReport {
    let count = |c: OutcomeClass| results.iter().filter(|r| r.outcome == c).count();
    let hazardous: Vec<_> = results.iter().filter(|r| r.hazardous).collect();
    let rate = |n: usize| if hazardous.is_empty() { 1.0 } else { n as f64 / hazardous.len() as f64 };
    let activated = hazardous.iter().filter(|r| r.decision.action != SafetyAction::Continue).count();
    let contained = hazardous.iter().filter(|r| r.decision.hazard_contained()).count();
    let times: Vec<f64> = results
        .iter()
        .filter_map(|r| r.time_to_safety_activation.map(|t| t as f64))
        .collect();
    let mut hist = [0usize; 4];
    for e in results.iter().flat_map(|r| r.events.iter()) {
        hist[match e.severity {
            Severity::Mild => 0,
            Severity::Moderate => 1,
            Severity::Severe => 2,
            Severity::Critical => 3,
        }] += 1;
    }
    let n = results.len().max(1) as f64;
    FailureCampaignReport {
        success: count(OutcomeClass::Success),
        partial: count(OutcomeClass::Partial),
        failure: count(OutcomeClass::Failure),
        catastrophic: count(OutcomeClass::Catastrophic),
        failure_rate: 1.0 - count(OutcomeClass::Success) as f64 / n,
        safety_activation_rate: rate(activated),
        containment_rate: rate(contained),
        mean_time_to_safety_activation: if times.is_empty() {
            None
        } else {
            Some(times.iter().sum::<f64>() / times.len() as f64)
        },
        severity_histogram: hist,
        results,
    }
}

/// Runs the full 8-scenario campaign (baseline + 7 failure scenarios).
pub fn run_failure_campaign() -> Result<FailureCampaignReport> {
    let results = vec![
        scenarios::run_regimen_scenario("Scenario 0: Balanced baseline", FailureMode::BalancedBaseline, 1.0, true),
        scenarios::run_regimen_scenario(
            "Scenario 1: OSK overdose (2.4x)",
            FailureMode::OskOverdose { dose_multiplier: 2.4 },
            1.0,
            true,
        ),
        scenarios::run_regimen_scenario(
            "Scenario 2: Telomerase without checkpoint",
            FailureMode::TelomeraseWithoutCheckpoint { tert_fold: 6.0, tp53_fold: 0.3 },
            1.0,
            true,
        ),
        scenarios::run_regimen_scenario(
            "Scenario 2b: Telomerase + suicide switch failure",
            FailureMode::SuicideSwitchFailure { tert_fold: 6.0, tp53_fold: 0.3 },
            1.0,
            false,
        ),
        scenarios::run_partial_organ_response_scenario(),
        scenarios::run_combined_severe_stress_scenario(66.0)?,
        scenarios::run_late_stage_scenario(100.0)?,
        scenarios::run_treatment_withdrawal_scenario(8)?,
    ];
    Ok(summarize_campaign(results))
}

#[cfg(test)]
mod tests {
    use super::*;
    use adverse::{max_severity, AdverseEventKind};
    use modes::balanced_regimen;

    #[test]
    fn test_balanced_baseline_is_not_harmed_by_safeguards() {
        let r = scenarios::run_regimen_scenario("base", FailureMode::BalancedBaseline, 1.0, true);
        assert_eq!(r.outcome, OutcomeClass::Success);
        assert!(r.events.is_empty(), "baseline raised events: {:?}", r.events);
        assert_eq!(r.decision.action, SafetyAction::Continue);
        let _ = balanced_regimen("x");
    }

    #[test]
    fn test_osk_overdose_triggers_guard_and_is_dose_dependent() {
        let r = scenarios::run_regimen_scenario(
            "od",
            FailureMode::OskOverdose { dose_multiplier: 2.4 },
            1.0,
            true,
        );
        assert!(r.events.iter().any(|e| e.kind == AdverseEventKind::TeratomaRisk));
        assert!(r.safety_activated && r.decision.hazard_contained());
        assert!(r.outcome >= OutcomeClass::Failure);
        assert!(r.time_to_safety_activation.is_some());
        // Identity collapsed to 0%, so functional rejuvenation & efficacy must be 0.0
        assert_eq!(r.efficacy_fraction, 0.0);

        let mut prev = None;
        let mut prev_eff = 1.0f64;
        for k in [0.0, 0.4, 0.7, 1.0] {
            let rr = scenarios::run_regimen_scenario(
                "dr",
                FailureMode::OskOverdose { dose_multiplier: 2.4 },
                k,
                true,
            );
            let s = max_severity(&rr.events);
            assert!(s >= prev, "severity must not decrease with dose");
            assert!(
                rr.efficacy_fraction <= prev_eff + 1e-6,
                "efficacy must decrease as identity erodes with overdose"
            );
            prev = s;
            prev_eff = rr.efficacy_fraction;
        }
    }

    #[test]
    fn test_stress_sweep_identifies_tipping_point_and_full_repair_shifts_threshold() {
        let partial_sweep = scenarios::run_combined_stress_sweep(0.35, 54.0).unwrap();
        let full_sweep = scenarios::run_combined_stress_sweep(0.95, 54.0).unwrap();

        assert_eq!(partial_sweep.first().unwrap().outcome, OutcomeClass::Success);
        assert!(partial_sweep.last().unwrap().outcome >= OutcomeClass::Failure);

        let partial_total_passed: usize = partial_sweep.iter().map(|p| p.criteria_passed).sum();
        let full_total_passed: usize = full_sweep.iter().map(|p| p.criteria_passed).sum();
        assert!(
            full_total_passed > partial_total_passed,
            "fully repaired twin must pass more stress sweep points than partially repaired twin"
        );

        let scen4 = scenarios::run_combined_severe_stress_scenario(54.0).unwrap();
        assert!(scen4.safety_activated);
        assert!(scen4.outcome >= OutcomeClass::Failure);
        assert!(scen4.time_to_safety_activation.is_some());
    }

    #[test]
    fn test_unconsolidated_treatment_withdrawal_triggers_senescence_rebound() {
        let r = scenarios::run_treatment_withdrawal_scenario(8).unwrap();
        assert!(r.events.iter().any(|e| e.kind == AdverseEventKind::SenescenceRebound));
        assert!(r.safety_activated);
        assert!(r.outcome >= OutcomeClass::Partial);
    }

    #[test]
    fn test_telomerase_without_checkpoint_raises_oncogenic_risk_and_guard_acts() {
        let r = scenarios::run_regimen_scenario(
            "tert",
            FailureMode::TelomeraseWithoutCheckpoint { tert_fold: 6.0, tp53_fold: 0.3 },
            1.0,
            true,
        );
        let clonal = r.events.iter().find(|e| e.kind == AdverseEventKind::ClonalExpansion).unwrap();
        assert!(clonal.magnitude >= scies_bio_th::aging::MAX_SAFE_ONCOGENIC_RISK_THRESHOLD);
        // Irreversible oncogenic hazard cannot be paused away: the switch must fire.
        assert_eq!(r.decision.action, SafetyAction::ActivateSuicideSwitch);
        assert!(r.decision.hazard_contained());
        assert_eq!(r.outcome, OutcomeClass::Failure);
    }

    #[test]
    fn test_suicide_switch_failure_is_catastrophic_and_flagged() {
        let r = scenarios::run_regimen_scenario(
            "sw",
            FailureMode::SuicideSwitchFailure { tert_fold: 6.0, tp53_fold: 0.3 },
            1.0,
            false,
        );
        assert!(r.events.iter().any(|e| e.kind == AdverseEventKind::ClonalExpansion));
        assert!(r.decision.safeguard_failure);
        assert_eq!(r.decision.action, SafetyAction::Abort);
        assert!(!r.decision.hazard_contained());
        assert_eq!(r.outcome, OutcomeClass::Catastrophic);
    }

    #[test]
    fn test_partial_response_targets_lagging_organs_not_dose() {
        let r = scenarios::run_partial_organ_response_scenario();
        assert!(r.events.iter().any(|e| e.kind == AdverseEventKind::SystemicImbalance));
        assert_eq!(r.decision.action, SafetyAction::TargetLaggingOrgans);
        assert_eq!(r.outcome, OutcomeClass::Partial);
    }

    #[test]
    fn test_late_stage_is_not_a_full_success() {
        let r = scenarios::run_late_stage_scenario(100.0).unwrap();
        assert!(r.outcome >= OutcomeClass::Partial);
    }

    #[test]
    fn test_failure_campaign_reports_non_trivial_outcomes() {
        let report = run_failure_campaign().unwrap();
        assert!(report.results.len() >= 6);
        // Not "everything passes": at least one failure and the baseline stays a success.
        assert_eq!(report.results[0].outcome, OutcomeClass::Success);
        assert!(report.failure + report.catastrophic >= 1);
        assert!(report.failure_rate > 0.0 && report.failure_rate < 1.0);
        // Every hazardous scenario must have triggered the guard.
        assert!((report.safety_activation_rate - 1.0).abs() < 1e-9);
        for r in &report.results {
            assert!(!r.decision.log.is_empty());
        }
    }
}
