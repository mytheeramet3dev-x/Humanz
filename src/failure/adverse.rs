//! Adverse-event definitions with graded severity.

use serde::{Deserialize, Serialize};

/// Severity grade (ordered: `Mild < Moderate < Severe < Critical`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Severity {
    Mild,
    Moderate,
    Severe,
    Critical,
}

/// Kinds of adverse events the failure suite can raise.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AdverseEventKind {
    TeratomaRisk,
    ClonalExpansion,
    ImmuneReactionToVector,
    HepaticToxicity,
    RenalToxicity,
    CardiacFibrosis,
    Neuroinflammation,
    MetabolicInstability,
    ThermoregulatoryFailure,
    SenescenceRebound,
    SystemicImbalance,
    IncompleteRejuvenation,
}

/// One adverse event.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AdverseEvent {
    pub kind: AdverseEventKind,
    pub organ: Option<String>,
    pub severity: Severity,
    pub reversible: bool,
    /// Raw measured value that produced the grade (units depend on `kind`).
    pub magnitude: f64,
    pub description: String,
}

/// Grades `value` against ascending thresholds `[mild, moderate, severe, critical]`.
pub fn grade(value: f64, t: [f64; 4]) -> Option<Severity> {
    if value >= t[3] {
        Some(Severity::Critical)
    } else if value >= t[2] {
        Some(Severity::Severe)
    } else if value >= t[1] {
        Some(Severity::Moderate)
    } else if value >= t[0] {
        Some(Severity::Mild)
    } else {
        None
    }
}

/// Highest severity in `events`, if any.
pub fn max_severity(events: &[AdverseEvent]) -> Option<Severity> {
    events.iter().map(|e| e.severity).max()
}
