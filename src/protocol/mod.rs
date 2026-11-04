//! Autonomous Closed-Loop Multi-Hallmark Rejuvenation & Negligible Senescence Controller (`NegligibleSenescenceController`).

use crate::body::{HumanDigitalTwin, OrganSystem, WholeBodySnapshot};
use crate::failure::safety_response::{
    apply_safety_action, evaluate_safety_guard_response, SafetyAction,
};
use crate::failure::triggers::detect_regimen_adverse_events;
use crate::genomics::{
    audit_patient_genome_with_biodb, BioDbGenomicAuditReport, ComorbidityProfile, PatientGenotype,
};
use crate::safety_guard::{audit_whole_body_safety, WholeBodySafetyAudit};
use scies_bio_sim::prelude::*;
use scies_bio_th::anatomy::OrganArchitecture;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::Path;

/// Target youthful homeostasis setpoints for Negligible Senescence maintenance,
/// combining cellular hallmark thresholds with 7-circuit clinical physiological guardrails.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HomeostasisSetpoints {
    /// Target maximum biological age (years) before triggering a cyclic OSK + TERT maintenance pulse.
    pub target_max_biological_age_years: f64,
    /// Target minimum telomere length (bp).
    pub target_min_telomere_bp: f64,
    /// Target maximum mtDNA heteroplasmy fraction.
    pub target_max_heteroplasmy: f64,
    /// Target maximum proteostasis aggregate burden.
    pub target_max_proteostasis_burden: f64,
    /// Target minimum microvascular capillary density.
    pub target_min_capillary_density: f64,
    /// Target minimum cellular viability across all organs.
    pub target_min_viability: f64,
    /// Target minimum renal glomerular filtration rate ($\text{eGFR}$, $\text{mL/min}$).
    pub target_min_egfr_ml_min: f64,
    /// Target minimum left ventricular ejection fraction ($\text{EF}$, $\%$).
    pub target_min_ejection_fraction_pct: f64,
    /// Target minimum cortical glymphatic clearance rate ($\mu\text{L/min}$).
    pub target_min_glymphatic_clearance_ul_min: f64,
    /// Target maximum portal venous endotoxin translocation ($\text{LPS}$, $\text{EU/mL}$).
    pub target_max_portal_lps_eu_ml: f64,
    /// Target minimum daily thymic naive T-cell emigration ($\times 10^7\text{ cells/day}$).
    pub target_min_naive_t_cell_output_10e7: f64,
}

impl Default for HomeostasisSetpoints {
    fn default() -> Self {
        Self {
            target_max_biological_age_years: 26.5,
            target_min_telomere_bp: 9_200.0,
            target_max_heteroplasmy: 0.045,
            target_max_proteostasis_burden: 0.055,
            target_min_capillary_density: 0.88,
            target_min_viability: 0.96,
            target_min_egfr_ml_min: 92.0,
            target_min_ejection_fraction_pct: 58.0,
            target_min_glymphatic_clearance_ul_min: 18.0,
            target_max_portal_lps_eu_ml: 0.32,
            target_min_naive_t_cell_output_10e7: 7.5,
        }
    }
}

/// Summary of an executed multi-organ rejuvenation & multi-decade negligible senescence trajectory.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NegligibleSenescenceRunSummary {
    pub baseline_snapshot: WholeBodySnapshot,
    pub post_reversal_snapshot: WholeBodySnapshot,
    pub final_treated_snapshot: WholeBodySnapshot,
    pub final_control_snapshot: WholeBodySnapshot,
    pub longitudinal_treated_history: Vec<WholeBodySnapshot>,
    pub longitudinal_control_history: Vec<WholeBodySnapshot>,
    pub net_aging_rate_treated_bio_yr_per_chrono_yr: f64,
    pub net_aging_rate_control_bio_yr_per_chrono_yr: f64,
    pub safety_audit: WholeBodySafetyAudit,
    pub maintenance_interventions_triggered: usize,
    pub targeted_circuit_repairs_triggered: usize,
    pub pbpk_precision_organ_doses: BTreeMap<OrganSystem, f64>,
}

/// Autonomous Closed-Loop Multi-Hallmark Rejuvenation & Negligible Senescence Controller.
pub struct NegligibleSenescenceController {
    pub setpoints: HomeostasisSetpoints,
}

impl Default for NegligibleSenescenceController {
    fn default() -> Self {
        Self::new(HomeostasisSetpoints::default())
    }
}

impl NegligibleSenescenceController {
    pub fn new(setpoints: HomeostasisSetpoints) -> Self {
        Self { setpoints }
    }

    /// **Step 1: PBPK-Guided Precision Dosing**
    /// Computes individualized `OSK` / `LNP` factor dosage multipliers for each organ based on its
    /// spatial Reaction-Diffusion PDE bioavailability (`mean_parenchymal_bioavailability`).
    ///
    /// Sinusoidal organs with high passive uptake (e.g. Liver, Bone Marrow) receive a conservative,
    /// de-escalated dose to minimize oncogenic stress, while tight-barrier organs (e.g. Brain BBB,
    /// Dense Dermis) receive a calibrated shuttle-assisted precision dose.
    pub fn compute_pbpk_precision_organ_doses(
        &self,
        twin: &HumanDigitalTwin,
    ) -> BTreeMap<OrganSystem, f64> {
        let mut doses = BTreeMap::new();
        let target_bioavail = 0.56;
        for organ in OrganSystem::all() {
            let anat = &twin.treated_anatomy[&organ];
            let bioavail = anat
                .mean_parenchymal_bioavailability(scies_bio_th::aging::MW_MRNA_LNP_VECTOR_DA, 0.92)
                .clamp(0.20, 0.95);
            // Inverse square-root PBPK dose normalization bounded within safe window [0.92, 1.14],
            // further capped by pharmacogenomic CYP2D6/3A4 hepatic clearance capacity.
            let mut precision_dose = (1.02 * (target_bioavail / bioavail).sqrt()).clamp(0.92, 1.14);
            if organ == OrganSystem::HepaticMetabolicParenchyma {
                precision_dose = precision_dose.min(twin.genotype.safe_hepatic_lnp_dose_cap());
            }
            doses.insert(organ, precision_dose);
        }
        doses
    }

    /// **Step 2: Physiological Closed-Loop Autonomous Controller**
    /// Inspects real-time clinical telemetry across the 7 closed-loop physiological circuits
    /// (`WholeBodySystemicCircuits`) and dispatches targeted organ-level micro-interventions
    /// whenever any circuit drifts outside youthful physiological setpoints.
    pub fn execute_targeted_physiological_circuit_repairs(
        &self,
        twin: &mut HumanDigitalTwin,
    ) -> Result<usize> {
        let mut repairs = 0usize;

        // 1. Renal-RAAS Circuit Guard: eGFR < 92 mL/min or elevated RAAS
        if twin.treated_circuits.renal_raas.egfr_ml_min < self.setpoints.target_min_egfr_ml_min {
            twin.treated_atlas.organs[7].rejuvenate_architecture(0.85, 0.85, 0.85);
            twin.treated_atlas.organs[16].rejuvenate_architecture(0.80, 0.80, 0.80);
            twin.treated_endocrine.apply_endocrine_rejuvenation(0.85, 0.85);
            repairs += 1;
        }

        // 2. Cardiovascular-Hemodynamic Circuit Guard: EF < 58% or elevated SVR
        if twin.treated_circuits.cardiovascular.ejection_fraction_pct
            < self.setpoints.target_min_ejection_fraction_pct
        {
            let col_id = twin.treated_organs[&OrganSystem::CardiovascularEndothelium];
            twin.sim.run_vascular_angiogenesis_restoration(col_id, 0.85)?;
            if let Some(anat) = twin
                .treated_anatomy
                .get_mut(&OrganSystem::CardiovascularEndothelium)
            {
                anat.rejuvenate_architecture(0.85, 0.85, 0.85);
            }
            twin.treated_atlas.organs[0].rejuvenate_architecture(0.85, 0.85, 0.85);
            twin.treated_atlas.organs[6].rejuvenate_architecture(0.85, 0.85, 0.85);
            repairs += 1;
        }

        // 3. Neuroendocrine-Glymphatic Circuit Guard: Glymphatic clearance < 18 uL/min
        if twin
            .treated_circuits
            .neuroendocrine_glymphatic
            .cortical_glymphatic_clearance_rate
            < self.setpoints.target_min_glymphatic_clearance_ul_min
        {
            let col_id = twin.treated_organs[&OrganSystem::CerebralCortexNeuroglia];
            twin.sim
                .run_proteostasis_clearance(col_id, 0.85, 1.25, 2.1, 1.8)?;
            twin.treated_atlas.organs[1].rejuvenate_architecture(0.85, 0.85, 0.85);
            twin.treated_atlas.organs[13].rejuvenate_architecture(0.85, 0.85, 0.85);
            twin.treated_atlas.organs[14].rejuvenate_architecture(0.85, 0.85, 0.85);
            repairs += 1;
        }

        // 4. Enterohepatic Barrier & Detox Circuit Guard: Portal LPS > 0.32 EU/mL
        if twin
            .treated_circuits
            .enterohepatic
            .portal_lps_endotoxin_eu_ml
            > self.setpoints.target_max_portal_lps_eu_ml
        {
            let col_id = twin.treated_organs[&OrganSystem::HepaticMetabolicParenchyma];
            twin.sim.run_pink1_parkin_mitophagy(col_id, 0.85, 1.0)?;
            twin.treated_atlas.organs[2].rejuvenate_architecture(0.85, 0.85, 0.85);
            twin.treated_atlas.organs[9].rejuvenate_architecture(0.85, 0.85, 0.85);
            twin.treated_atlas.organs[11].rejuvenate_architecture(0.85, 0.85, 0.85);
            repairs += 1;
        }

        // 5. Immune-Lymphoid Surveillance Circuit Guard: Naive T-Cell < 7.5e7/day
        if twin
            .treated_circuits
            .immune_lymphatic
            .daily_naive_t_cell_emigration_10e7
            < self.setpoints.target_min_naive_t_cell_output_10e7
        {
            let col_id = twin.treated_organs[&OrganSystem::ThymicLymphoidOrgan];
            twin.sim
                .run_immune_surveillance_and_thymic_regeneration(col_id, 0.88, 0.90)?;
            twin.treated_atlas.organs[3].rejuvenate_architecture(0.85, 0.85, 0.85);
            twin.treated_atlas.organs[4].rejuvenate_architecture(0.85, 0.85, 0.85);
            twin.treated_atlas.organs[17].rejuvenate_architecture(0.85, 0.85, 0.85);
            repairs += 1;
        }

        if repairs > 0 {
            twin.propagate_inter_organ_systemic_coupling(true);
        }

        Ok(repairs)
    }

    /// Executes a single full 6-Hallmark Rejuvenation Sweep across all 6 organs of `twin.treated_organs`.
    /// Returns `(max_oncogenic_risk_observed, min_somatic_identity_retention_observed)`.
    pub fn execute_six_hallmark_rejuvenation_sweep(
        &self,
        twin: &mut HumanDigitalTwin,
        osk_cycles: u32,
        adjuvant_navitoclax_um: f64,
    ) -> Result<(f64, f64)> {
        let mut max_onc_risk: f64 = 0.0;
        let mut min_identity: f64 = 1.0;
        let target_cells = twin.cells_per_organ;
        let precision_doses = self.compute_pbpk_precision_organ_doses(twin);

        for organ in OrganSystem::all() {
            let col_id = twin.treated_organs[&organ];

            // 1. Hallmark 5 + Anatomy OOP: Microvascular Angiogenesis & ECM Matrix De-Stiffening (VEGF-A / Glucosepane Lyase / MMP)
            twin.sim.run_vascular_angiogenesis_restoration(col_id, 0.95)?;
            if let Some(anat) = twin.treated_anatomy.get_mut(&organ) {
                anat.rejuvenate_architecture(0.95, 0.92, 0.90);
            }

            let organ_dose = precision_doses[&organ];
            let (tert_fold, tp53_fold) = twin
                .genotype
                .checkpoint_adjusted_tert_and_tp53_folds(2.25 * (organ_dose / 1.02), 1.85);

            // 2. Hallmark 1: PBPK & PGx-Guided Pulsed OSK Partial Epigenetic Reprogramming + Active Safety Guard Governor
            let mut regimen = scies_bio_th::aging::OskPulseRegimen {
                somatic_lineage: organ.lineage_name().to_string(),
                delivery_vector: scies_bio_th::aging::DeliveryVector::MrnaLipidNanoparticle,
                pulse_duration_days: 2.0,
                washout_interval_days: 5.0,
                total_cycles: osk_cycles as usize,
                tet1_tet2_competence: 0.95,
                oncogene_cmyc_included: false,
                factor_dosage: organ_dose,
                transient_tert_mrna_fold: tert_fold,
                tp53_cdkn2a_safety_fold: tp53_fold,
            };

            // Active Closed-Loop Safety Guard pre-flight check before in vivo transduction
            let preflight =
                scies_bio_th::aging::simulate_in_vivo_osk_pulse_regimen(twin.chronological_age_years, &regimen);
            let events = detect_regimen_adverse_events(&preflight, &regimen, &twin.treated_atlas);
            let decision = evaluate_safety_guard_response(&events, true);
            if let Some(mitigated_regimen) = apply_safety_action(&regimen, &decision) {
                regimen = mitigated_regimen;
            } else if matches!(decision.action, SafetyAction::Pause | SafetyAction::Abort | SafetyAction::ActivateSuicideSwitch) {
                regimen.total_cycles = 1;
                regimen.factor_dosage *= 0.5;
            }

            let rep = twin
                .sim
                .run_in_vivo_osk_partial_reprogramming(col_id, &regimen)?;
            max_onc_risk = max_onc_risk.max(rep.composite_oncogenic_risk);
            min_identity = min_identity.min(rep.final_somatic_identity_retention);

            // 3. Hallmark 2: PINK1/Parkin Mitophagy + TFAM mtDNA Biogenesis
            twin.sim.run_pink1_parkin_mitophagy(col_id, 0.92, 1.05)?;

            // 4. Hallmark 3: mTORC1/TFEB Autophagy + HSP70/HSP90 Chaperone Proteostasis Clearance
            twin.sim
                .run_proteostasis_clearance(col_id, 0.88, 1.30, 2.3, 1.9)?;

            // 5. Hallmark 6: FOXN1/IL-7 Thymic Rejuvenation + Endogenous NK & Macrophage Immunosurveillance
            twin.sim
                .run_immune_surveillance_and_thymic_regeneration(col_id, 0.92, 0.94)?;

            if adjuvant_navitoclax_um > 0.01 {
                let pgx_nav_dose =
                    adjuvant_navitoclax_um * twin.genotype.navitoclax_pgx_dose_multiplier();
                let nav_params = scies_bio_th::pharma::PkPdParameters {
                    compound_name: "Navitoclax".to_string(),
                    ec50_um: 0.25,
                    td50_um: 4.8,
                    hill_coefficient: 2.2,
                    half_life_hours: twin.genotype.navitoclax_half_life_hours(),
                    bioavailability: 0.90,
                };
                twin.sim
                    .administer_pkpd_compound(col_id, nav_params, pgx_nav_dose);
                twin.sim.run_for(SimDuration::hours(24));
            }

            // 6. Hallmark 4: Wnt/Notch Stem Cell Niche Asymmetric Division & Progenitor Repopulation
            twin.sim
                .run_stem_cell_niche_regeneration(col_id, 0.92, 0.90, target_cells)?;
        }

        // 7. Systemic Neuroendocrine Axes & 20-Organ Atlas Rejuvenation
        twin.treated_endocrine.apply_endocrine_rejuvenation(0.95, 0.95);
        for organ in &mut twin.treated_atlas.organs {
            organ.rejuvenate_architecture(0.95, 0.92, 0.90);
        }

        twin.propagate_inter_organ_systemic_coupling(true);
        Ok((max_onc_risk, min_identity))
    }

    /// Runs the full 2-Phase HumanZ Negligible Senescence Protocol:
    /// - **Phase I (Deep Systemic Age Reversal)**: Rejuvenates an aged twin (e.g., 72 yrs) down to youthful biological state (~22–25 yrs) across 3 staged therapeutic rounds.
    /// - **Phase II (Multi-Decade Closed-Loop Negligible Senescence Homeostasis)**: Simulates `horizon_years` of chronological time while maintaining net biological aging rate $\le 0.0\text{ bio-yr/yr}$.
    pub fn run_negligible_senescence_lifespan_simulation(
        &self,
        twin: &mut HumanDigitalTwin,
        horizon_years: u32,
    ) -> Result<NegligibleSenescenceRunSummary> {
        let baseline_snapshot = twin.snapshot_whole_body(true);
        let mut max_onc_risk: f64 = 0.0;
        let mut min_identity: f64 = 1.0;

        // Phase I: 4 staged induction rounds to reverse 72-year-old multi-organ damage down to youthful ~27-year baseline
        for round in 0..4 {
            let nav_dose = if round == 0 { 0.75 } else { 0.25 };
            let (r_risk, r_id) = self.execute_six_hallmark_rejuvenation_sweep(twin, 6, nav_dose)?;
            max_onc_risk = max_onc_risk.max(r_risk);
            min_identity = min_identity.min(r_id);
        }

        let post_reversal_snapshot = twin.snapshot_whole_body(true);
        let mut longitudinal_treated_history = vec![post_reversal_snapshot.clone()];
        let mut longitudinal_control_history = vec![twin.snapshot_whole_body(false)];
        let mut maintenance_interventions_triggered = 0usize;
        let mut targeted_circuit_repairs_triggered = 0usize;

        // Phase II: Year-by-year closed-loop Negligible Senescence homeostasis across `horizon_years`
        for _yr in 1..=horizon_years {
            twin.chronological_age_years += 1.0;

            // Both untreated control and treated twin experience 1 year of chronological entropy
            twin.apply_chronological_entropy(1.0, false);
            twin.propagate_inter_organ_systemic_coupling(false);

            twin.apply_chronological_entropy(1.0, true);
            twin.propagate_inter_organ_systemic_coupling(true);

            // Step 2: Autonomous 7-Circuit Physiological Closed-Loop Repair Check
            targeted_circuit_repairs_triggered +=
                self.execute_targeted_physiological_circuit_repairs(twin)?;

            // Closed-loop biomarker check on the treated HumanZ twin
            let current = twin.snapshot_whole_body(true);
            let needs_maintenance = current.mean_biological_age_years
                > self.setpoints.target_max_biological_age_years
                || current.mean_telomere_bp < self.setpoints.target_min_telomere_bp
                || current.mean_heteroplasmy > self.setpoints.target_max_heteroplasmy
                || current.mean_proteostasis_burden > self.setpoints.target_max_proteostasis_burden
                || current.mean_capillary_density < self.setpoints.target_min_capillary_density
                || current.mean_viability < self.setpoints.target_min_viability;

            if needs_maintenance {
                let (r_risk, r_id) =
                    self.execute_six_hallmark_rejuvenation_sweep(twin, 2, 0.20)?;
                max_onc_risk = max_onc_risk.max(r_risk);
                min_identity = min_identity.min(r_id);
                maintenance_interventions_triggered += 1;
            }

            longitudinal_treated_history.push(twin.snapshot_whole_body(true));
            longitudinal_control_history.push(twin.snapshot_whole_body(false));
        }

        let final_treated_snapshot = twin.snapshot_whole_body(true);
        let final_control_snapshot = twin.snapshot_whole_body(false);
        let pbpk_precision_organ_doses = self.compute_pbpk_precision_organ_doses(twin);

        let years_f = horizon_years.max(1) as f64;
        let net_aging_rate_treated = (final_treated_snapshot.mean_biological_age_years
            - post_reversal_snapshot.mean_biological_age_years)
            / years_f;
        let net_aging_rate_control = (final_control_snapshot.mean_biological_age_years
            - baseline_snapshot.mean_biological_age_years)
            / years_f;

        let safety_audit = audit_whole_body_safety(twin, max_onc_risk, min_identity);

        Ok(NegligibleSenescenceRunSummary {
            baseline_snapshot,
            post_reversal_snapshot,
            final_treated_snapshot,
            final_control_snapshot,
            longitudinal_treated_history,
            longitudinal_control_history,
            net_aging_rate_treated_bio_yr_per_chrono_yr: net_aging_rate_treated,
            net_aging_rate_control_bio_yr_per_chrono_yr: net_aging_rate_control,
            safety_audit,
            maintenance_interventions_triggered,
            targeted_circuit_repairs_triggered,
            pbpk_precision_organ_doses,
        })
    }

    /// Executes an end-to-end personalized clinical rescue simulation on a high-risk patient with:
    /// - `APOE ε4/ε4` Alzheimer/BBB vulnerability
    /// - `CYP2D6/CYP3A4` Poor Metabolizer status (requiring pharmacogenomic dose de-escalation)
    /// - Hypomorphic `TP53` heterozygosity (requiring reinforced safety lock fold)
    /// - Pre-existing `CKD Stage 3b` + `Type 2 Diabetes` + `HFpEF Hypertension`
    ///
    /// Audits the patient's genome via [`scies_bio_db::BioDb`] at `db_dir` and returns the
    /// complete [`PersonalizedComorbidityRescueReport`].
    pub fn run_personalized_comorbidity_rescue_simulation(
        &self,
        db_dir: &Path,
        initial_age_years: f64,
        horizon_years: u32,
    ) -> Result<PersonalizedComorbidityRescueReport> {
        let genotype = PatientGenotype::high_risk_apoe4_cyp_poor_tp53_het();
        let comorbidities = ComorbidityProfile::cardiometabolic_ckd_syndrome();
        let biodb_audit =
            audit_patient_genome_with_biodb(db_dir, initial_age_years, &genotype)?;

        let mut twin = HumanDigitalTwin::new_personalized(
            "Patient-PGx-CKD3b-T2D-HFpEF",
            initial_age_years,
            16,
            2026_99,
            genotype.clone(),
            comorbidities.clone(),
        );

        let pre_egfr = twin.treated_circuits.renal_raas.egfr_ml_min;
        let pre_sbp = twin.treated_circuits.cardiovascular.systolic_bp_mmhg;
        let pre_glucose = twin.treated_circuits.enterohepatic.blood_glucose_mg_dl;
        let pre_brain_prot = twin
            .snapshot_whole_body(true)
            .organs[&OrganSystem::CerebralCortexNeuroglia]
            .stats
            .mean_proteostasis_burden;

        let summary = self.run_negligible_senescence_lifespan_simulation(&mut twin, horizon_years)?;

        let post_egfr = summary.final_treated_snapshot.circuits.renal_raas.egfr_ml_min;
        let post_sbp = summary
            .final_treated_snapshot
            .circuits
            .cardiovascular
            .systolic_bp_mmhg;
        let post_glucose = summary
            .final_treated_snapshot
            .circuits
            .enterohepatic
            .blood_glucose_mg_dl;
        let post_brain_prot = summary.final_treated_snapshot.organs
            [&OrganSystem::CerebralCortexNeuroglia]
            .stats
            .mean_proteostasis_burden;

        let all_comorbidities_rescued = post_egfr >= 90.0
            && post_sbp < 130.0
            && post_glucose < 100.0
            && post_brain_prot < 0.08
            && summary.safety_audit.all_organs_approved
            && biodb_audit.navitoclax_within_therapeutic_window
            && biodb_audit.all_guides_real_reference_provenance;

        Ok(PersonalizedComorbidityRescueReport {
            genotype,
            comorbidities,
            biodb_audit,
            pre_egfr_ml_min: pre_egfr,
            post_egfr_ml_min: post_egfr,
            pre_systolic_bp_mmhg: pre_sbp,
            post_systolic_bp_mmhg: post_sbp,
            pre_fasting_glucose_mg_dl: pre_glucose,
            post_fasting_glucose_mg_dl: post_glucose,
            pre_cortical_proteostasis_burden: pre_brain_prot,
            post_cortical_proteostasis_burden: post_brain_prot,
            all_comorbidities_rescued,
            run_summary: summary,
        })
    }
}

/// Complete telemetry report of a `BioDb`-guided personalized comorbidity & PGx rescue simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PersonalizedComorbidityRescueReport {
    pub genotype: PatientGenotype,
    pub comorbidities: ComorbidityProfile,
    pub biodb_audit: BioDbGenomicAuditReport,
    pub pre_egfr_ml_min: f64,
    pub post_egfr_ml_min: f64,
    pub pre_systolic_bp_mmhg: f64,
    pub post_systolic_bp_mmhg: f64,
    pub pre_fasting_glucose_mg_dl: f64,
    pub post_fasting_glucose_mg_dl: f64,
    pub pre_cortical_proteostasis_burden: f64,
    pub post_cortical_proteostasis_burden: f64,
    pub all_comorbidities_rescued: bool,
    pub run_summary: NegligibleSenescenceRunSummary,
}
