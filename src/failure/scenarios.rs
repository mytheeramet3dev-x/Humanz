//! Ready-made failure scenarios, parametric stress sweeps, and campaign evaluation.

use super::adverse::{grade, max_severity, AdverseEvent, AdverseEventKind, Severity};
use super::modes::{balanced_regimen, inject_failure_mode, FailureMode, OutcomeClass};
use super::safety_response::{
    apply_safety_action, evaluate_safety_guard_response, SafetyAction, SafetyDecision,
};
use super::triggers::{
    detect_regimen_adverse_events, BBB_TEER_DEFICIT_BANDS_OHM_CM2, HYPOTHERMIA_DEFICIT_BANDS_C,
    KDIGO_EGFR_DEFICIT_BANDS_ML_MIN, OGTT_2H_GLUCOSE_BANDS_MG_DL,
    SENESCENCE_REBOUND_RATIO_BANDS, SYSTEMIC_IMBALANCE_SPREAD_BANDS,
};
use crate::body::{HumanDigitalTwin, OrganSystem};
use crate::protocol::{HomeostasisSetpoints, NegligibleSenescenceController};
use scies_bio_sim::prelude::*;
use scies_bio_th::aging::simulate_in_vivo_osk_pulse_regimen;
use scies_bio_th::anatomy::{
    AnatomicalOrgan, AnatomicalSpatialAtlas, OrganArchitecture, WholeBodySystemicCircuits,
};
use serde::{Deserialize, Serialize};

/// One point of a failure trajectory or parametric stress sweep.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TrajectoryPoint {
    pub step: usize,
    pub primary_value: f64,
    pub hazard: f64,
    pub safe: bool,
}

/// Telemetry of a single stress level $s \in [0.0, 1.0]$ during the Scenario 4 Stress Sweep.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StressSweepPoint {
    pub stress_intensity: f64,
    pub map_mmhg: f64,
    pub ambient_temp_c: f64,
    pub vo2_max_ml_kg_min: f64,
    pub ogtt_2h_glucose_mg_dl: f64,
    pub defended_core_temp_c: f64,
    pub stressed_egfr_ml_min: f64,
    pub stressed_bbb_teer_ohm_cm2: f64,
    pub criteria_passed: usize,
    pub efficacy: f64,
    pub events: Vec<AdverseEvent>,
    pub outcome: OutcomeClass,
}

/// Result of one failure scenario.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FailureScenarioResult {
    pub name: String,
    pub mode: FailureMode,
    pub events: Vec<AdverseEvent>,
    pub decision: SafetyDecision,
    pub outcome: OutcomeClass,
    /// Guard did something other than `Continue`.
    pub safety_activated: bool,
    /// At least a Moderate event occurred.
    pub hazardous: bool,
    /// Cycle/step at which the guard first saw the system leave its safe window.
    pub time_to_safety_activation: Option<usize>,
    pub efficacy_fraction: f64,
    pub trajectory: Vec<TrajectoryPoint>,
    /// Outcome after applying the guard's dose reduction (when the action allows continuing).
    pub mitigated_outcome: Option<OutcomeClass>,
    pub notes: Vec<String>,
}

/// Classifies an outcome from the worst severity, containment, and efficacy.
pub fn classify_outcome(
    max_sev: Option<Severity>,
    contained: bool,
    efficacy: f64,
) -> OutcomeClass {
    match max_sev {
        Some(Severity::Critical) | Some(Severity::Severe) => {
            if contained {
                OutcomeClass::Failure
            } else {
                OutcomeClass::Catastrophic
            }
        }
        Some(Severity::Moderate) => {
            if contained {
                OutcomeClass::Partial
            } else {
                OutcomeClass::Failure
            }
        }
        Some(Severity::Mild) | None => {
            if efficacy >= 0.7 {
                OutcomeClass::Success
            } else {
                OutcomeClass::Partial
            }
        }
    }
}

fn finish(
    name: &str,
    mode: FailureMode,
    events: Vec<AdverseEvent>,
    decision: SafetyDecision,
    efficacy: f64,
    trajectory: Vec<TrajectoryPoint>,
    time_to_safety: Option<usize>,
    mitigated: Option<OutcomeClass>,
    notes: Vec<String>,
) -> FailureScenarioResult {
    let max_sev = max_severity(&events);
    let outcome = classify_outcome(max_sev, decision.hazard_contained(), efficacy);
    FailureScenarioResult {
        name: name.to_string(),
        mode,
        safety_activated: decision.action != SafetyAction::Continue,
        hazardous: max_sev.map_or(false, |s| s >= Severity::Moderate),
        events,
        decision,
        outcome,
        time_to_safety_activation: time_to_safety,
        efficacy_fraction: efficacy,
        trajectory,
        mitigated_outcome: mitigated,
        notes,
    }
}

/// Runs a regimen-level scenario (overdose, telomerase, switch failure, baseline).
pub fn run_regimen_scenario(
    name: &str,
    mode: FailureMode,
    intensity: f64,
    suicide_switch_functional: bool,
) -> FailureScenarioResult {
    let initial_age = 72.0;
    let atlas = AnatomicalSpatialAtlas::from_somatic_age(26.0);
    let baseline = balanced_regimen("Fibroblast");
    let baseline_report = simulate_in_vivo_osk_pulse_regimen(initial_age, &baseline);

    let mut regimen = baseline.clone();
    inject_failure_mode(&mut regimen, &mode, intensity);
    let report = simulate_in_vivo_osk_pulse_regimen(initial_age, &regimen);

    let events = detect_regimen_adverse_events(&report, &regimen, &atlas);
    let decision = evaluate_safety_guard_response(&events, suicide_switch_functional);

    // Therapeutic efficacy scales with identity-coupled biological age reversal and is
    // penalized if the guard had to ablate the transduced clone (iCasp9) or abort.
    let ablation_factor = match decision.action {
        SafetyAction::ActivateSuicideSwitch => 0.15,
        SafetyAction::Abort => 0.0,
        _ => 1.0,
    };
    let efficacy = ((report.total_rejuvenation_years
        / baseline_report.total_rejuvenation_years.max(1e-6))
        * report.final_somatic_identity_retention.clamp(0.0, 1.0)
        * ablation_factor)
    .clamp(0.0, 1.0);

    let trajectory: Vec<TrajectoryPoint> = report
        .per_cycle_metrics
        .iter()
        .map(|c| TrajectoryPoint {
            step: c.cycle,
            primary_value: c.somatic_identity_retention,
            hazard: c.teratoma_hazard_score,
            safe: c.in_safe_operating_window,
        })
        .collect();
    let time_to_safety = report
        .per_cycle_metrics
        .iter()
        .find(|c| !c.in_safe_operating_window)
        .map(|c| c.cycle);

    // Closed-loop mitigation: re-run with the guard's reduced regimen, if it allows continuing.
    let mitigated = apply_safety_action(&regimen, &decision).map(|r| {
        let rep = simulate_in_vivo_osk_pulse_regimen(initial_age, &r);
        let ev = detect_regimen_adverse_events(&rep, &r, &atlas);
        let dec = evaluate_safety_guard_response(&ev, suicide_switch_functional);
        let eff = ((rep.total_rejuvenation_years
            / baseline_report.total_rejuvenation_years.max(1e-6))
            * rep.final_somatic_identity_retention)
            .clamp(0.0, 1.0);
        classify_outcome(max_severity(&ev), dec.hazard_contained() || ev.is_empty(), eff)
    });

    let notes = vec![format!(
        "verdict={:?}, rejuvenation={:.1}y (baseline {:.1}y), identity={:.0}%",
        report.safety_verdict,
        report.total_rejuvenation_years,
        baseline_report.total_rejuvenation_years,
        report.final_somatic_identity_retention * 100.0
    )];
    finish(name, mode, events, decision, efficacy, trajectory, time_to_safety, mitigated, notes)
}

fn steady_state(atlas: &AnatomicalSpatialAtlas) -> WholeBodySystemicCircuits {
    let mut c = WholeBodySystemicCircuits::new();
    for _ in 0..30 {
        c.step_synchronized(atlas, 1.0, 22.0, 75.0, 0.10);
    }
    c
}

/// Scenario 3: only heart/artery/liver/marrow/thymus/skin are rejuvenated; kidney, lung, brain,
/// choroid, hippocampus and adrenal remain aged. Measures cross-system discordance.
pub fn run_partial_organ_response_scenario() -> FailureScenarioResult {
    let aged_age = 95.0;
    let treat = |atlas: &mut AnatomicalSpatialAtlas, idx: &[usize]| {
        for &i in idx {
            for _ in 0..2 {
                atlas.organs[i].rejuvenate_architecture(0.95, 0.92, 0.90);
            }
        }
    };
    let mut full = AnatomicalSpatialAtlas::from_somatic_age(aged_age);
    let all: Vec<usize> = (0..20).collect();
    treat(&mut full, &all);
    let mut partial = AnatomicalSpatialAtlas::from_somatic_age(aged_age);
    treat(&mut partial, &[0, 6, 2, 3, 4, 5, 9, 10]);

    let c_full = steady_state(&full);
    let c_part = steady_state(&partial);
    let heart_norm = (c_part.cardiovascular.ejection_fraction_pct
        / c_full.cardiovascular.ejection_fraction_pct)
        .clamp(0.0, 1.2);
    let kidney_norm =
        (c_part.renal_raas.egfr_ml_min / c_full.renal_raas.egfr_ml_min).clamp(0.0, 1.2);
    let brain_norm = (c_part.neuroendocrine_glymphatic.cortical_glymphatic_clearance_rate
        / c_full.neuroendocrine_glymphatic.cortical_glymphatic_clearance_rate)
        .clamp(0.0, 1.2);
    let norms = [heart_norm, kidney_norm, brain_norm];
    let spread = norms.iter().cloned().fold(f64::MIN, f64::max)
        - norms.iter().cloned().fold(f64::MAX, f64::min);
    let efficacy = (norms.iter().sum::<f64>() / 3.0).clamp(0.0, 1.0);

    let mut events = Vec::new();
    if let Some(s) = grade(spread, SYSTEMIC_IMBALANCE_SPREAD_BANDS) {
        events.push(AdverseEvent {
            kind: AdverseEventKind::SystemicImbalance,
            organ: None,
            severity: s,
            reversible: true,
            magnitude: spread,
            description: format!(
                "Heart {:.0}% / Kidney {:.0}% / Glymphatic {:.0}% of fully-rejuvenated function",
                heart_norm * 100.0,
                kidney_norm * 100.0,
                brain_norm * 100.0
            ),
        });
    }
    let decision = evaluate_safety_guard_response(&events, true);
    let traj = norms
        .iter()
        .enumerate()
        .map(|(i, v)| TrajectoryPoint {
            step: i + 1,
            primary_value: *v,
            hazard: spread,
            safe: *v >= 0.7,
        })
        .collect();
    finish(
        "Scenario 3: Partial organ response",
        FailureMode::PartialOrganResponse,
        events,
        decision,
        efficacy,
        traj,
        None,
        None,
        vec![format!("function spread across heart/kidney/brain = {:.2}", spread)],
    )
}

/// Evaluates all 5 acute physiological resilience axes on `twin` at a normalized combined stress
/// intensity $s \in [0.0, 1.0]$ and target hypotensive floor `min_map_mmhg`.
pub fn evaluate_stress_at_intensity(
    twin: &HumanDigitalTwin,
    stress_intensity: f64,
    min_map_mmhg: f64,
) -> StressSweepPoint {
    let s = stress_intensity.clamp(0.0, 1.0);
    let atlas = &twin.treated_atlas;
    let circuits = &twin.treated_circuits;
    let snap = twin.snapshot_whole_body(true);

    let heart = match &atlas.organs[6] {
        AnatomicalOrgan::MyocardialVentricular(o) => o,
        _ => unreachable!(),
    };
    let lung = match &atlas.organs[8] {
        AnatomicalOrgan::PulmonaryAlveolar(o) => o,
        _ => unreachable!(),
    };
    let muscle = match &atlas.organs[10] {
        AnatomicalOrgan::SkeletalMuscleFascicle(o) => o,
        _ => unreachable!(),
    };
    let gut = match &atlas.organs[11] {
        AnatomicalOrgan::GastrointestinalCrypt(o) => o,
        _ => unreachable!(),
    };
    let pancreas = match &atlas.organs[9] {
        AnatomicalOrgan::PancreaticIslet(o) => o,
        _ => unreachable!(),
    };
    let liver = match &atlas.organs[2] {
        AnatomicalOrgan::HepaticLobule(o) => o,
        _ => unreachable!(),
    };
    let bat = match &atlas.organs[19] {
        AnatomicalOrgan::ThermogenicBrownAdipose(o) => o,
        _ => unreachable!(),
    };
    let skin = match &atlas.organs[5] {
        AnatomicalOrgan::DermalCutaneous(o) => o,
        _ => unreachable!(),
    };
    let kidney = match &atlas.organs[7] {
        AnatomicalOrgan::RenalNephron(o) => o,
        _ => unreachable!(),
    };
    let cortex = match &atlas.organs[1] {
        AnatomicalOrgan::CerebralCortex(o) => o,
        _ => unreachable!(),
    };

    // 1. Cardiopulmonary + Hypoxic Exercise Challenge:
    let fio2 = 0.21 - 0.045 * s;
    let peak_hr = (208.0 - 0.7 * snap.mean_biological_age_years).clamp(115.0, 192.0);
    let beat = heart.simulate_cardiac_cycle(peak_hr, 16.5 + 3.5 * s, 1.15 * (1.0 - 0.15 * s));
    let filling = (0.55 + 0.55 * beat.diastolic_compliance_index).clamp(0.62, 1.15);
    let septic_fatigue = 1.0 - 0.28 * s;
    let peak_co = ((beat.stroke_volume_ml * filling * peak_hr / 1000.0) * 1.48 * septic_fatigue)
        .clamp(4.0, 22.5);
    let gas = lung.simulate_respiratory_gas_exchange(1050.0, 36.0, fio2, peak_co);
    let o2_ext = (0.55 + 0.35 * (muscle.vascular.capillary_density / 0.90)).clamp(0.45, 0.92);
    let vo2_max_ml_kg_min =
        (((peak_co * 0.20 * (gas.arterial_o2_saturation_pct / 100.0) * o2_ext) * 1000.0) / 70.0)
            .clamp(10.0, 65.0);

    // 2. Metabolic OGTT + Stress Cortisol/Catecholamine Insulin Resistance:
    let gut_res =
        gut.simulate_digestion_and_barrier(750.0 * (1.0 + 0.45 * s), 28.0 * (1.0 - 0.65 * s));
    let panc_res =
        pancreas.simulate_islet_secretion(180.0 + 45.0 * s, gut_res.glp1_incretin_secretion_pmol_l);
    let mus_res =
        muscle.simulate_locomotor_contraction(0.35, panc_res.insulin_secretion_rate_pmol_l);
    let hep_res = liver.simulate_zonated_metabolism(
        1.0 + 0.4 * s,
        panc_res.insulin_secretion_rate_pmol_l,
        1.0 + 0.5 * s,
        circuits.enterohepatic.systemic_ammonia_umol_l,
    );
    let stress_insulin_resistance = 1.0 - 0.52 * s;
    let ins_sens = ((mus_res.glucose_disposal_rate_mg_min / 210.0) * stress_insulin_resistance)
        .clamp(0.14, 1.25);
    let disp_idx = (panc_res.beta_cell_functional_mass * ins_sens).clamp(0.10, 1.25);
    let peak_glu = 92.0 + (48.0 + 26.0 * s) / disp_idx.sqrt();
    let clearance = 24.0
        + 34.0 * disp_idx
        + (160.0 - hep_res.hepatic_glucose_output_mg_min).max(0.0) * 0.12;
    let ogtt_2h = (peak_glu - clearance).max(78.0);

    // 3. Cold Shock + Exhaustive Convective Drain (+4 °C down to -10 °C as s -> 1.0):
    let ambient_temp_c = 4.0 - 14.0 * s;
    let bat_res = bat.simulate_thermogenesis(2.2 * (1.0 - 0.15 * s), ambient_temp_c.max(0.0));
    let mus_therm = muscle.simulate_locomotor_contraction(0.75 * (1.0 - 0.25 * s), 45.0);
    let skin_res = skin.simulate_barrier_and_thermoregulation(ambient_temp_c.max(0.0), 36.8);
    let total_heat =
        bat_res.non_shivering_thermogenesis_watts + mus_therm.mechanical_heat_generation_watts;
    let cold_drain = 82.0 + 38.0 * s + 0.05 * skin_res.cutaneous_heat_dissipation_flow_ml_min;
    let core_temp = (36.8 + (total_heat - cold_drain) * 0.028).clamp(34.0, 37.4);

    // 4 & 5. Hypotensive Renal Perfusion + Endotoxic Neurovascular BBB Challenge:
    let map_mmhg = 75.0 - (75.0 - min_map_mmhg) * s;
    let stressed_co = circuits.cardiovascular.cardiac_output_l_min * (1.0 - 0.30 * s);
    let barotrauma = snap.hemodynamics.microvascular_barotrauma_index + 0.18 * s;
    let renal = kidney.simulate_renal_filtration(map_mmhg, stressed_co, barotrauma);

    let acute_sasp = (snap.systemic_inflammaging_sasp + 0.35 + 0.52 * s).clamp(0.10, 0.98);
    let cort = cortex.simulate_cortical_computation(map_mmhg, 0.90 - 0.25 * s, acute_sasp);

    let mut events = Vec::new();
    let mut passed = 0usize;

    if vo2_max_ml_kg_min >= 38.0 {
        passed += 1;
    } else if let Some(sev) = grade(38.0 - vo2_max_ml_kg_min, [1.0, 5.0, 10.0, 16.0]) {
        events.push(AdverseEvent {
            kind: AdverseEventKind::CardiacFibrosis,
            organ: Some("Heart/Lung".into()),
            severity: sev,
            reversible: true,
            magnitude: vo2_max_ml_kg_min,
            description: format!("VO2max {:.1} mL/kg/min under combined stress (s={:.1})", vo2_max_ml_kg_min, s),
        });
    }

    if let Some(sev) = grade(ogtt_2h, OGTT_2H_GLUCOSE_BANDS_MG_DL) {
        events.push(AdverseEvent {
            kind: AdverseEventKind::MetabolicInstability,
            organ: Some("Pancreas/Liver".into()),
            severity: sev,
            reversible: true,
            magnitude: ogtt_2h,
            description: format!("OGTT 2h glucose {:.0} mg/dL under stress (s={:.1})", ogtt_2h, s),
        });
    } else {
        passed += 1;
    }

    if core_temp >= 36.0 {
        passed += 1;
    } else if let Some(sev) = grade(36.3 - core_temp, HYPOTHERMIA_DEFICIT_BANDS_C) {
        events.push(AdverseEvent {
            kind: AdverseEventKind::ThermoregulatoryFailure,
            organ: Some("BAT/Muscle".into()),
            severity: sev,
            reversible: true,
            magnitude: core_temp,
            description: format!(
                "Core temp {:.2} °C at {:.1} °C ambient cold (s={:.1})",
                core_temp, ambient_temp_c, s
            ),
        });
    }

    if renal.egfr_ml_min >= 60.0 {
        passed += 1;
    } else if let Some(sev) = grade(60.0 - renal.egfr_ml_min, KDIGO_EGFR_DEFICIT_BANDS_ML_MIN) {
        events.push(AdverseEvent {
            kind: AdverseEventKind::RenalToxicity,
            organ: Some("Kidney".into()),
            severity: sev,
            reversible: true,
            magnitude: renal.egfr_ml_min,
            description: format!(
                "eGFR {:.1} mL/min at MAP {:.0} mmHg (s={:.1})",
                renal.egfr_ml_min, map_mmhg, s
            ),
        });
    }

    if cort.bbb_teer_ohm_cm2 >= 1100.0 {
        passed += 1;
    } else if let Some(sev) = grade(1100.0 - cort.bbb_teer_ohm_cm2, BBB_TEER_DEFICIT_BANDS_OHM_CM2)
    {
        events.push(AdverseEvent {
            kind: AdverseEventKind::Neuroinflammation,
            organ: Some("Brain".into()),
            severity: sev,
            reversible: true,
            magnitude: cort.bbb_teer_ohm_cm2,
            description: format!(
                "BBB TEER {:.0} Ohm*cm2 under endotoxemia SASP {:.2} (s={:.1})",
                cort.bbb_teer_ohm_cm2, acute_sasp, s
            ),
        });
    }

    let efficacy = passed as f64 / 5.0;
    let dec = evaluate_safety_guard_response(&events, true);
    let outcome = classify_outcome(max_severity(&events), dec.hazard_contained() || events.is_empty(), efficacy);

    StressSweepPoint {
        stress_intensity: s,
        map_mmhg,
        ambient_temp_c,
        vo2_max_ml_kg_min,
        ogtt_2h_glucose_mg_dl: ogtt_2h,
        defended_core_temp_c: core_temp,
        stressed_egfr_ml_min: renal.egfr_ml_min,
        stressed_bbb_teer_ohm_cm2: cort.bbb_teer_ohm_cm2,
        criteria_passed: passed,
        efficacy,
        events,
        outcome,
    }
}

/// Builds a twin at `repair_fraction` in `[0.0, 1.0]` (`0.35` = early 1-round partial repair,
/// `0.95` = full 4-round consolidation) and evaluates the 6-level Stress Sweep $s \in \{0.0, 0.2, 0.4, 0.6, 0.8, 1.0\}$.
pub fn run_combined_stress_sweep(
    repair_fraction: f64,
    min_map_mmhg: f64,
) -> Result<Vec<StressSweepPoint>> {
    let mut twin = HumanDigitalTwin::new("Failure-StressSweep", 72.0, 12, 4242);
    let controller = NegligibleSenescenceController::new(HomeostasisSetpoints::default());
    let f = repair_fraction.clamp(0.10, 1.0);

    if f >= 0.80 {
        for round in 0..4 {
            controller.execute_six_hallmark_rejuvenation_sweep(
                &mut twin,
                6,
                if round == 0 { 0.75 } else { 0.25 },
            )?;
        }
    } else {
        // Only ONE early repair sweep (2 OSK cycles): passes baseline stress (s=0.0),
        // but breaks down as combined stressor intensity s ramps toward 1.0.
        controller.execute_six_hallmark_rejuvenation_sweep(&mut twin, 2, 0.20)?;
    }

    let levels = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0];
    Ok(levels
        .iter()
        .map(|&s| evaluate_stress_at_intensity(&twin, s, min_map_mmhg))
        .collect())
}

/// Scenario 4: executes a 6-step parametric Stress Sweep ($s \in [0.0, 1.0]$) on an incompletely
/// repaired twin (`repair_fraction = 0.35`), identifying the exact tipping point $s^*$ where
/// multi-organ physiological reserve breaks down.
pub fn run_combined_severe_stress_scenario(hypotension_map_mmhg: f64) -> Result<FailureScenarioResult> {
    let min_map = hypotension_map_mmhg.min(58.0);
    let sweep = run_combined_stress_sweep(0.35, min_map)?;

    let tipping_step = sweep
        .iter()
        .position(|p| p.outcome != OutcomeClass::Success)
        .map(|idx| idx + 1);
    let tipping_s = sweep
        .iter()
        .find(|p| p.outcome != OutcomeClass::Success)
        .map(|p| p.stress_intensity);
    let failure_s = sweep
        .iter()
        .find(|p| p.outcome >= OutcomeClass::Failure)
        .map(|p| p.stress_intensity);

    let traj: Vec<TrajectoryPoint> = sweep
        .iter()
        .enumerate()
        .map(|(i, p)| TrajectoryPoint {
            step: i + 1,
            primary_value: p.efficacy,
            hazard: p.stress_intensity,
            safe: p.outcome == OutcomeClass::Success,
        })
        .collect();

    let peak = sweep.last().unwrap().clone();
    let decision = evaluate_safety_guard_response(&peak.events, true);
    let curve_summary: Vec<String> = sweep
        .iter()
        .map(|p| format!("s={:.1}:{}/5({:?})", p.stress_intensity, p.criteria_passed, p.outcome))
        .collect();

    let notes = vec![
        format!("stress sweep [{}]", curve_summary.join(", ")),
        format!(
            "tipping point s* = {:?} (breakdown to Failure at s = {:?}); peak stress eGFR={:.1}, Core={:.2}°C, 2hGlu={:.0}, BBB={:.0}",
            tipping_s,
            failure_s,
            peak.stressed_egfr_ml_min,
            peak.defended_core_temp_c,
            peak.ogtt_2h_glucose_mg_dl,
            peak.stressed_bbb_teer_ohm_cm2
        ),
    ];

    Ok(finish(
        "Scenario 4: Combined severe stress (sweep s=0->1)",
        FailureMode::CombinedSevereStress { hypotension_map_mmhg },
        peak.events,
        decision,
        peak.efficacy,
        traj,
        tipping_step,
        None,
        notes,
    ))
}

/// Scenario 5: therapy begins at very advanced age / damage.
pub fn run_late_stage_scenario(initial_age_years: f64) -> Result<FailureScenarioResult> {
    let mut twin = HumanDigitalTwin::new("Failure-Late", initial_age_years, 12, 777);
    let controller = NegligibleSenescenceController::new(HomeostasisSetpoints::default());
    let summary = controller.run_negligible_senescence_lifespan_simulation(&mut twin, 3)?;
    let start = summary.baseline_snapshot.mean_biological_age_years;
    let post = summary.post_reversal_snapshot.mean_biological_age_years;
    let efficacy = ((start - post) / (start - 27.0).max(1.0)).clamp(0.0, 1.0);

    let mut events = Vec::new();
    if let Some(s) = grade(post - 27.0, [3.0, 8.0, 15.0, 30.0]) {
        events.push(AdverseEvent {
            kind: AdverseEventKind::IncompleteRejuvenation,
            organ: None,
            severity: s,
            reversible: true,
            magnitude: post,
            description: format!("Bio age {:.1}y -> {:.1}y (target <= 27y)", start, post),
        });
    }
    if !summary.safety_audit.all_organs_approved {
        events.push(AdverseEvent {
            kind: AdverseEventKind::ClonalExpansion,
            organ: None,
            severity: Severity::Severe,
            reversible: false,
            magnitude: summary.safety_audit.max_oncogenic_risk,
            description: "Whole-body safety audit rejected".into(),
        });
    }
    let decision = evaluate_safety_guard_response(&events, true);
    let traj = summary
        .longitudinal_treated_history
        .iter()
        .enumerate()
        .map(|(i, s)| TrajectoryPoint {
            step: i,
            primary_value: s.mean_biological_age_years,
            hazard: 0.0,
            safe: s.mean_biological_age_years <= 30.0,
        })
        .collect();
    Ok(finish(
        "Scenario 5: Late-stage high damage",
        FailureMode::LateStageHighDamage { initial_age_years },
        events,
        decision,
        efficacy,
        traj,
        None,
        None,
        vec![format!("start age {:.0}y, bio age {:.1}y -> {:.1}y", initial_age_years, start, post)],
    ))
}

/// Scenario 6: therapy stops abruptly after an unconsolidated / OSK-heavy reprogramming induction
/// (epigenetic clock reset without full ECM matrix de-stiffening or senolytic clearance).
/// Residual `YAP/TAZ` mechanotransduction, `ROS`, and paracrine `SASP` drive accelerated
/// post-withdrawal senescence rebound relative to untreated control.
pub fn run_treatment_withdrawal_scenario(years_off_therapy: u32) -> Result<FailureScenarioResult> {
    let mut twin = HumanDigitalTwin::new("Failure-Withdrawal", 72.0, 12, 99);

    // Unconsolidated OSK-dominant induction: resets epigenetic age across 3 rounds of 6-cycle OSK
    // pulses (raising chromatin plasticity), but performs only minimal (0.15) ECM/mitochondrial
    // consolidation and omits senolytics so residual YAP/TAZ, ROS, and SASP remain in the niche.
    for _round in 0..3 {
        for organ in OrganSystem::all() {
            let col_id = twin.treated_organs[&organ];
            if let Some(anat) = twin.treated_anatomy.get_mut(&organ) {
                anat.rejuvenate_architecture(0.15, 0.15, 0.15);
            }
            let mut reg = balanced_regimen(organ.lineage_name());
            reg.factor_dosage = 1.15;
            reg.pulse_duration_days = 2.4;
            reg.total_cycles = 6;
            twin.sim.run_in_vivo_osk_partial_reprogramming(col_id, &reg)?;
        }
    }
    twin.propagate_inter_organ_systemic_coupling(true);

    let start_t = twin.snapshot_whole_body(true).mean_biological_age_years;
    let start_c = twin.snapshot_whole_body(false).mean_biological_age_years;
    let mut traj = vec![TrajectoryPoint {
        step: 0,
        primary_value: start_t,
        hazard: 1.0,
        safe: true,
    }];
    let mut first_unsafe_year = None;

    for y in 1..=years_off_therapy {
        let prev_t = twin.snapshot_whole_body(true).mean_biological_age_years;
        let prev_c = twin.snapshot_whole_body(false).mean_biological_age_years;
        twin.chronological_age_years += 1.0;
        twin.apply_chronological_entropy(1.0, true);
        twin.propagate_inter_organ_systemic_coupling(true);
        twin.apply_chronological_entropy(1.0, false);
        twin.propagate_inter_organ_systemic_coupling(false);
        let cur_t = twin.snapshot_whole_body(true).mean_biological_age_years;
        let cur_c = twin.snapshot_whole_body(false).mean_biological_age_years;
        let step_ratio = (cur_t - prev_t) / (cur_c - prev_c).max(1e-6);
        let safe = step_ratio < SENESCENCE_REBOUND_RATIO_BANDS[1];
        if !safe && first_unsafe_year.is_none() {
            first_unsafe_year = Some(y as usize);
        }
        traj.push(TrajectoryPoint {
            step: y as usize,
            primary_value: cur_t,
            hazard: step_ratio,
            safe,
        });
    }

    let end_t = twin.snapshot_whole_body(true).mean_biological_age_years;
    let end_c = twin.snapshot_whole_body(false).mean_biological_age_years;
    let yrs = years_off_therapy.max(1) as f64;
    let rate_t = (end_t - start_t) / yrs;
    let rate_c = ((end_c - start_c) / yrs).max(1e-6);
    let ratio = rate_t / rate_c;
    let retained = ((end_c - end_t) / (end_c - start_t).max(1.0)).clamp(0.0, 1.0);

    let mut events = Vec::new();
    if let Some(s) = grade(ratio, SENESCENCE_REBOUND_RATIO_BANDS) {
        events.push(AdverseEvent {
            kind: AdverseEventKind::SenescenceRebound,
            organ: None,
            severity: s,
            reversible: true,
            magnitude: ratio,
            description: format!(
                "Post-withdrawal aging {:.2} y/y vs control {:.2} y/y (rebound x{:.2})",
                rate_t, rate_c, ratio
            ),
        });
    }
    let decision = evaluate_safety_guard_response(&events, true);
    Ok(finish(
        "Scenario 6: Treatment withdrawal (rebound check)",
        FailureMode::TreatmentWithdrawal { years_off_therapy },
        events,
        decision,
        retained,
        traj,
        first_unsafe_year,
        None,
        vec![format!(
            "rebound rate ratio treated/control = {:.2}x ({:.2} vs {:.2} bio-y/y); bio age {:.1}y -> {:.1}y (control {:.1}y)",
            ratio, rate_t, rate_c, start_t, end_t, end_c
        )],
    ))
}
