//! Safety Guard response policy: maps adverse events to a `SafetyAction` and keeps a decision log.

use super::adverse::{max_severity, AdverseEvent, AdverseEventKind, Severity};
use scies_bio_th::aging::OskPulseRegimen;
use serde::{Deserialize, Serialize};

/// Action taken by the Safety Guard.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum SafetyAction {
    Continue,
    ReduceDose { factor: f64 },
    TargetLaggingOrgans,
    Pause,
    ActivateSuicideSwitch,
    Abort,
}

/// A logged Safety Guard decision.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SafetyDecision {
    pub action: SafetyAction,
    pub rationale: String,
    pub triggered_by: Vec<AdverseEventKind>,
    /// True when the guard needed a safeguard (e.g. suicide switch) that did not work.
    pub safeguard_failure: bool,
    pub log: Vec<String>,
}

impl SafetyDecision {
    /// True when the hazard was actually neutralised (action taken and safeguards worked).
    pub fn hazard_contained(&self) -> bool {
        self.action != SafetyAction::Continue && !self.safeguard_failure
    }
}

fn is_oncogenic(k: AdverseEventKind) -> bool {
    matches!(k, AdverseEventKind::TeratomaRisk | AdverseEventKind::ClonalExpansion)
}

fn is_coverage_issue(k: AdverseEventKind) -> bool {
    matches!(
        k,
        AdverseEventKind::SystemicImbalance | AdverseEventKind::IncompleteRejuvenation
    )
}

/// Policy:
/// - none / Mild → `Continue` (monitor)
/// - Moderate → `TargetLaggingOrgans` for coverage problems, otherwise `ReduceDose(0.6)`
/// - Severe → `Pause` (pausing cannot undo an event, so see the next rule)
/// - Critical **or Severe-and-irreversible** oncogenic event → `ActivateSuicideSwitch`; if the
///   switch is not functional → `Abort` with `safeguard_failure = true`
/// - Critical other → `Abort`
pub fn evaluate_safety_guard_response(
    events: &[AdverseEvent],
    suicide_switch_functional: bool,
) -> SafetyDecision {
    let mut log: Vec<String> = events
        .iter()
        .map(|e| {
            format!(
                "[{:?}] {:?}{} reversible={} :: {}",
                e.severity,
                e.kind,
                e.organ.as_ref().map(|o| format!("@{}", o)).unwrap_or_default(),
                e.reversible,
                e.description
            )
        })
        .collect();

    let triggered_by: Vec<AdverseEventKind> = events
        .iter()
        .filter(|e| e.severity >= Severity::Moderate)
        .map(|e| e.kind)
        .collect();

    // An irreversible oncogenic event at Severe or above cannot be "paused away": a transformed
    // clone persists, so the only effective containment is elimination of the cells.
    let needs_elimination = events.iter().any(|e| {
        is_oncogenic(e.kind)
            && (e.severity == Severity::Critical
                || (e.severity == Severity::Severe && !e.reversible))
    });

    let mut safeguard_failure = false;
    let (action, rationale) = if needs_elimination {
        if suicide_switch_functional {
            (
                SafetyAction::ActivateSuicideSwitch,
                "Irreversible/critical oncogenic hazard: iCasp9 suicide switch activated".to_string(),
            )
        } else {
            safeguard_failure = true;
            (
                SafetyAction::Abort,
                "Irreversible/critical oncogenic hazard and suicide switch NOT functional: abort, hazard uncontained"
                    .to_string(),
            )
        }
    } else {
        match max_severity(events) {
            None => (SafetyAction::Continue, "No adverse events".to_string()),
            Some(Severity::Mild) => (
                SafetyAction::Continue,
                "Mild events only: continue with enhanced monitoring".to_string(),
            ),
            Some(Severity::Moderate) => {
                let only_coverage = events
                    .iter()
                    .filter(|e| e.severity >= Severity::Moderate)
                    .all(|e| is_coverage_issue(e.kind));
                if only_coverage {
                    (
                        SafetyAction::TargetLaggingOrgans,
                        "Moderate coverage deficit: extend repair to lagging organs".to_string(),
                    )
                } else {
                    (
                        SafetyAction::ReduceDose { factor: 0.6 },
                        "Moderate toxicity: reduce dose to 60%".to_string(),
                    )
                }
            }
            Some(Severity::Severe) => (
                SafetyAction::Pause,
                "Severe event: pause therapy pending review".to_string(),
            ),
            Some(Severity::Critical) => {
                (SafetyAction::Abort, "Critical event: abort therapy".to_string())
            }
        }
    };
    log.push(format!("DECISION {:?} :: {}", action, rationale));

    SafetyDecision {
        action,
        rationale,
        triggered_by,
        safeguard_failure,
        log,
    }
}

/// Builds the mitigated regimen implied by the decision, if the decision allows continuing.
pub fn apply_safety_action(
    regimen: &OskPulseRegimen,
    decision: &SafetyDecision,
) -> Option<OskPulseRegimen> {
    match decision.action {
        SafetyAction::ReduceDose { factor } => {
            let mut r = regimen.clone();
            r.factor_dosage *= factor;
            r.pulse_duration_days = (r.pulse_duration_days * factor).clamp(1.5, 2.5);
            r.transient_tert_mrna_fold = (r.transient_tert_mrna_fold * factor).max(1.0);
            r.tp53_cdkn2a_safety_fold = r.tp53_cdkn2a_safety_fold.max(1.8);
            Some(r)
        }
        _ => None,
    }
}
