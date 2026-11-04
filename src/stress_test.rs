//! # Clinical Stress & Acute Physiological Challenge Suite (`stress_test`)
//!
//! Subjects both the untreated `Control` twin and the rejuvenated `HumanZ` twin (at the end of the
//! multi-decade simulation, e.g. Chronological Age 102 years) to 4 standardized, non-destructive
//! acute physiological stress tests (executed on cloned twin states so baseline state is preserved):
//!
//! 1. **Cardiopulmonary & Locomotor $\text{VO}_{2\max}$ Exercise Stress Test** (Bruce Protocol Ramp)
//! 2. **Oral Glucose Tolerance & Enterohepatic Portal Load Test (`OGTT` + Endotoxin Challenge)**
//! 3. **Acute Cold Exposure ($4^\circ\text{C}$) & Thermoregulatory Defense Challenge**
//! 4. **Systemic Endotoxemia, Acute Renal Perfusion & Neuro-Glymphatic Recovery Challenge**

use crate::body::HumanDigitalTwin;
use scies_bio_th::anatomy::AnatomicalOrgan;
use serde::{Deserialize, Serialize};

/// Telemetry of a single arm (`Control` or `HumanZ`) during the Cardiopulmonary $\text{VO}_{2\max}$ Exercise Challenge.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ExerciseStressMetrics {
    pub peak_heart_rate_bpm: f64,
    pub peak_stroke_volume_ml: f64,
    pub peak_cardiac_output_l_min: f64,
    pub vo2_max_ml_kg_min: f64,
    pub exercise_arterial_sao2_pct: f64,
    pub peak_muscle_contractile_force_n_cm2: f64,
    pub peak_irisin_myokine_ng_ml: f64,
}

/// Telemetry of a single arm during the Oral Glucose Tolerance & Enterohepatic Portal Challenge (`OGTT`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OgttMetabolicChallengeMetrics {
    pub peak_postprandial_glucose_mg_dl: f64,
    pub two_hour_post_load_glucose_mg_dl: f64,
    pub peak_biphasic_insulin_pmol_l: f64,
    pub hepatic_urea_detox_efficiency_pct: f64,
    pub hepatic_cyp450_clearance_multiplier: f64,
    pub portal_lps_endotoxin_leak_eu_ml: f64,
}

/// Telemetry of a single arm during the $4^\circ\text{C}$ Acute Cold Exposure Challenge.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ColdThermoregulatoryChallengeMetrics {
    pub defended_core_temperature_c: f64,
    pub bat_ucp1_non_shivering_heat_watts: f64,
    pub muscle_thermogenic_heat_watts: f64,
    pub total_thermogenic_power_watts: f64,
    pub hypothermia_avoided: bool,
}

/// Telemetry of a single arm during the Systemic Endotoxemia & Neuro-Renal Recovery Challenge.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EndotoxemiaNeuroRenalChallengeMetrics {
    pub stressed_egfr_ml_min: f64,
    pub acute_kidney_injury_avoided: bool,
    pub stressed_bbb_teer_ohm_cm2: f64,
    pub sws_glymphatic_clearance_ul_min: f64,
    pub hippocampal_neurogenesis_cells_day: f64,
    pub immune_patrolling_competence: f64,
}

/// Complete comparative report across all 4 clinical stress challenges (`Control` vs `HumanZ`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ClinicalStressChallengeReport {
    pub chronological_age_years: f64,
    pub control_exercise: ExerciseStressMetrics,
    pub humanz_exercise: ExerciseStressMetrics,
    pub control_ogtt: OgttMetabolicChallengeMetrics,
    pub humanz_ogtt: OgttMetabolicChallengeMetrics,
    pub control_cold: ColdThermoregulatoryChallengeMetrics,
    pub humanz_cold: ColdThermoregulatoryChallengeMetrics,
    pub control_endotoxemia: EndotoxemiaNeuroRenalChallengeMetrics,
    pub humanz_endotoxemia: EndotoxemiaNeuroRenalChallengeMetrics,
    pub all_humanz_resilience_criteria_passed: bool,
}

/// Runs the complete 4-Stage Clinical Stress & Acute Physiological Challenge Suite on `twin`.
pub fn run_clinical_stress_challenge_suite(
    twin: &HumanDigitalTwin,
) -> ClinicalStressChallengeReport {
    let control_exercise = evaluate_exercise_stress_arm(twin, false);
    let humanz_exercise = evaluate_exercise_stress_arm(twin, true);

    let control_ogtt = evaluate_ogtt_metabolic_arm(twin, false);
    let humanz_ogtt = evaluate_ogtt_metabolic_arm(twin, true);

    let control_cold = evaluate_cold_exposure_arm(twin, false);
    let humanz_cold = evaluate_cold_exposure_arm(twin, true);

    let control_endotoxemia = evaluate_endotoxemia_neuro_renal_arm(twin, false);
    let humanz_endotoxemia = evaluate_endotoxemia_neuro_renal_arm(twin, true);

    let all_passed = humanz_exercise.vo2_max_ml_kg_min >= 38.0
        && humanz_exercise.peak_cardiac_output_l_min > control_exercise.peak_cardiac_output_l_min * 1.8
        && humanz_ogtt.two_hour_post_load_glucose_mg_dl < 140.0
        && humanz_cold.hypothermia_avoided
        && humanz_endotoxemia.acute_kidney_injury_avoided
        && humanz_endotoxemia.stressed_bbb_teer_ohm_cm2 > 1100.0;

    ClinicalStressChallengeReport {
        chronological_age_years: twin.chronological_age_years,
        control_exercise,
        humanz_exercise,
        control_ogtt,
        humanz_ogtt,
        control_cold,
        humanz_cold,
        control_endotoxemia,
        humanz_endotoxemia,
        all_humanz_resilience_criteria_passed: all_passed,
    }
}

fn evaluate_exercise_stress_arm(
    twin: &HumanDigitalTwin,
    is_treated: bool,
) -> ExerciseStressMetrics {
    let atlas = if is_treated {
        &twin.treated_atlas
    } else {
        &twin.control_atlas
    };
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

    // Maximal sympathetic exercise drive: target HR = 208 - 0.7 * effective_bio_age
    // In geriatric HFpEF (low diastolic compliance), rapid tachycardia shortens diastole and impairs ventricular filling
    let snap = twin.snapshot_whole_body(is_treated);
    let peak_hr = (208.0 - 0.7 * snap.mean_biological_age_years).clamp(115.0, 192.0);
    let beat = heart.simulate_cardiac_cycle(peak_hr, 16.5, 1.15);
    let tachycardia_filling_factor =
        (0.55 + 0.55 * beat.diastolic_compliance_index).clamp(0.62, 1.15);
    let peak_sv = beat.stroke_volume_ml * tachycardia_filling_factor;
    let peak_co = ((peak_sv * peak_hr / 1000.0) * 1.48).clamp(4.0, 22.5);

    // Exercise hyperpnea (RR = 36 bpm, Vt = 1050 mL)
    let gas = lung.simulate_respiratory_gas_exchange(1050.0, 36.0, 0.21, peak_co);

    // Maximal voluntary & sympathetic motor unit recruitment (0.95)
    let mus = muscle.simulate_locomotor_contraction(0.95, 65.0);

    // Fick Principle VO2max (mL/kg/min for 70 kg adult): CO * arteriovenous O2 difference * extraction efficiency
    let o2_extraction = (0.55 + 0.35 * (muscle.vascular.capillary_density / 0.90)).clamp(0.45, 0.92);
    let vo2_max_l_min =
        peak_co * 0.20 * (gas.arterial_o2_saturation_pct / 100.0) * o2_extraction;
    let vo2_max_ml_kg_min = ((vo2_max_l_min * 1000.0) / 70.0).clamp(10.0, 65.0);

    ExerciseStressMetrics {
        peak_heart_rate_bpm: peak_hr,
        peak_stroke_volume_ml: beat.stroke_volume_ml,
        peak_cardiac_output_l_min: peak_co,
        vo2_max_ml_kg_min,
        exercise_arterial_sao2_pct: gas.arterial_o2_saturation_pct,
        peak_muscle_contractile_force_n_cm2: mus.specific_contractile_force_n_cm2,
        peak_irisin_myokine_ng_ml: mus.irisin_secretion_ng_ml,
    }
}

fn evaluate_ogtt_metabolic_arm(
    twin: &HumanDigitalTwin,
    is_treated: bool,
) -> OgttMetabolicChallengeMetrics {
    let atlas = if is_treated {
        &twin.treated_atlas
    } else {
        &twin.control_atlas
    };
    let circuits = if is_treated {
        &twin.treated_circuits
    } else {
        &twin.control_circuits
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
    let muscle = match &atlas.organs[10] {
        AnatomicalOrgan::SkeletalMuscleFascicle(o) => o,
        _ => unreachable!(),
    };

    // 1. Acute 75g Oral Glucose + Lipid Meal Challenge:
    // Gut incretin & barrier response under standard 28g prebiotic fiber buffering
    let gut_res = gut.simulate_digestion_and_barrier(750.0, 28.0);
    let panc_res =
        pancreas.simulate_islet_secretion(180.0, gut_res.glp1_incretin_secretion_pmol_l);
    let mus_res =
        muscle.simulate_locomotor_contraction(0.35, panc_res.insulin_secretion_rate_pmol_l);
    let hep_res = liver.simulate_zonated_metabolism(
        1.0,
        panc_res.insulin_secretion_rate_pmol_l,
        1.0,
        circuits.enterohepatic.systemic_ammonia_umol_l,
    );

    // 1-Hour Peak Postprandial Glucose (ADA Reference: Youthful 125-155 mg/dL; Aged IGT/T2D 195-235 mg/dL)
    let insulin_sensitivity = (mus_res.glucose_disposal_rate_mg_min / 210.0).clamp(0.28, 1.25);
    let beta_disposition_index =
        (panc_res.beta_cell_functional_mass * insulin_sensitivity).clamp(0.14, 1.25);
    let peak_glucose = 92.0 + 48.0 / beta_disposition_index.sqrt();

    // 2-Hour Post-Load Glucose (ADA Reference: Normal < 140 mg/dL [youthful ~88-105]; IGT 140-199 mg/dL)
    let two_hr_clearance = 24.0 + 34.0 * beta_disposition_index
        + (160.0 - hep_res.hepatic_glucose_output_mg_min).max(0.0) * 0.12;
    let two_hr_glucose = (peak_glucose - two_hr_clearance).max(78.0);

    OgttMetabolicChallengeMetrics {
        peak_postprandial_glucose_mg_dl: peak_glucose,
        two_hour_post_load_glucose_mg_dl: two_hr_glucose,
        peak_biphasic_insulin_pmol_l: panc_res.insulin_secretion_rate_pmol_l,
        hepatic_urea_detox_efficiency_pct: hep_res.urea_cycle_clearance_efficiency_pct,
        hepatic_cyp450_clearance_multiplier: hep_res.cyp_detoxification_capacity_multiplier,
        portal_lps_endotoxin_leak_eu_ml: gut_res.portal_lps_endotoxin_translocation_eu_ml,
    }
}

fn evaluate_cold_exposure_arm(
    twin: &HumanDigitalTwin,
    is_treated: bool,
) -> ColdThermoregulatoryChallengeMetrics {
    let atlas = if is_treated {
        &twin.treated_atlas
    } else {
        &twin.control_atlas
    };
    let bat = match &atlas.organs[19] {
        AnatomicalOrgan::ThermogenicBrownAdipose(o) => o,
        _ => unreachable!(),
    };
    let muscle = match &atlas.organs[10] {
        AnatomicalOrgan::SkeletalMuscleFascicle(o) => o,
        _ => unreachable!(),
    };
    let skin = match &atlas.organs[5] {
        AnatomicalOrgan::DermalCutaneous(o) => o,
        _ => unreachable!(),
    };

    // Acute 4.0 °C cold shock with intense beta-3 sympathetic norepinephrine drive (2.2)
    let bat_res = bat.simulate_thermogenesis(2.2, 4.0);
    let mus_res = muscle.simulate_locomotor_contraction(0.75, 45.0);
    let skin_res = skin.simulate_barrier_and_thermoregulation(4.0, 36.8);

    let total_heat =
        bat_res.non_shivering_thermogenesis_watts + mus_res.mechanical_heat_generation_watts;
    let cold_drain_watts = 82.0 + 0.05 * skin_res.cutaneous_heat_dissipation_flow_ml_min;
    let core_temp = (36.8 + (total_heat - cold_drain_watts) * 0.028).clamp(34.2, 37.4);

    ColdThermoregulatoryChallengeMetrics {
        defended_core_temperature_c: core_temp,
        bat_ucp1_non_shivering_heat_watts: bat_res.non_shivering_thermogenesis_watts,
        muscle_thermogenic_heat_watts: mus_res.mechanical_heat_generation_watts,
        total_thermogenic_power_watts: total_heat,
        hypothermia_avoided: core_temp >= 36.0,
    }
}

fn evaluate_endotoxemia_neuro_renal_arm(
    twin: &HumanDigitalTwin,
    is_treated: bool,
) -> EndotoxemiaNeuroRenalChallengeMetrics {
    let atlas = if is_treated {
        &twin.treated_atlas
    } else {
        &twin.control_atlas
    };
    let mut circuits = if is_treated {
        twin.treated_circuits.clone()
    } else {
        twin.control_circuits.clone()
    };
    let snap = twin.snapshot_whole_body(is_treated);

    let cortex = match &atlas.organs[1] {
        AnatomicalOrgan::CerebralCortex(o) => o,
        _ => unreachable!(),
    };
    let kidney = match &atlas.organs[7] {
        AnatomicalOrgan::RenalNephron(o) => o,
        _ => unreachable!(),
    };
    let hypo = match &atlas.organs[12] {
        AnatomicalOrgan::HypothalamicNeuroendocrine(o) => o,
        _ => unreachable!(),
    };
    let hippo = match &atlas.organs[13] {
        AnatomicalOrgan::HippocampalDentateGyrus(o) => o,
        _ => unreachable!(),
    };
    let choroid = match &atlas.organs[14] {
        AnatomicalOrgan::ChoroidPlexusCsf(o) => o,
        _ => unreachable!(),
    };
    let pit = match &atlas.organs[15] {
        AnatomicalOrgan::PituitaryHypophyseal(o) => o,
        _ => unreachable!(),
    };
    let adrenal = match &atlas.organs[16] {
        AnatomicalOrgan::AdrenalSuprarenal(o) => o,
        _ => unreachable!(),
    };
    let marrow = match &atlas.organs[3] {
        AnatomicalOrgan::BoneMarrowNiche(o) => o,
        _ => unreachable!(),
    };
    let thymus = match &atlas.organs[4] {
        AnatomicalOrgan::ThymicLymphoid(o) => o,
        _ => unreachable!(),
    };
    let spleen = match &atlas.organs[17] {
        AnatomicalOrgan::SplenicLymphoidPulp(o) => o,
        _ => unreachable!(),
    };

    // Acute endotoxic / inflammatory surge (+0.38 SASP) + hypotensive MAP challenge (74 mmHg)
    let acute_sasp = (snap.systemic_inflammaging_sasp + 0.38).clamp(0.10, 0.98);
    let barotrauma = snap.hemodynamics.microvascular_barotrauma_index;
    let renal = kidney.simulate_renal_filtration(74.0, circuits.cardiovascular.cardiac_output_l_min, barotrauma);
    let cort = cortex.simulate_cortical_computation(74.0, 0.90, acute_sasp);

    // Advance into Slow-Wave Sleep (SWS, circadian hour 2.0 AM) to measure glymphatic clearance & neurogenesis
    circuits.neuroendocrine_glymphatic.circadian_hour = 1.0;
    circuits.neuroendocrine_glymphatic.step(
        hypo,
        pit,
        adrenal,
        choroid,
        cortex,
        hippo,
        1.0,
        acute_sasp,
    );
    circuits.immune_lymphatic.step(marrow, thymus, spleen);

    EndotoxemiaNeuroRenalChallengeMetrics {
        stressed_egfr_ml_min: renal.egfr_ml_min,
        acute_kidney_injury_avoided: renal.egfr_ml_min >= 60.0,
        stressed_bbb_teer_ohm_cm2: cort.bbb_teer_ohm_cm2,
        sws_glymphatic_clearance_ul_min: circuits
            .neuroendocrine_glymphatic
            .cortical_glymphatic_clearance_rate,
        hippocampal_neurogenesis_cells_day: circuits
            .neuroendocrine_glymphatic
            .adult_hippocampal_neurogenesis_cells_day,
        immune_patrolling_competence: circuits
            .immune_lymphatic
            .peripheral_immune_patrolling_competence,
    }
}
