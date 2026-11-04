//! # HumanZ (`humanz`)
//!
//! **Whole-Body Multi-Organ Negligible Senescence & Autonomous Rejuvenation Digital Twin Platform.**
//!
//! Built atop [`scies_bio_th`], [`scies_bio_db`], and [`scies_bio_sim`], `humanz` models 6 interconnected
//! human organ compartments (`CardiovascularEndothelium`, `CerebralCortexNeuroglia`,
//! `HepaticMetabolicParenchyma`, `HematopoieticBoneMarrowNiche`, `ThymicLymphoidOrgan`,
//! `DermalEpithelialTissue`) coupled via systemic circulatory SASP, vascular oxygenation,
//! and immune surveillance.

pub mod body;
pub mod export;
pub mod failure;
pub mod fm_index;
pub mod genome_store;
pub mod genomics;
pub mod protocol;
pub mod safety_guard;
pub mod sais;
pub mod sha256;
pub mod stress_test;
pub mod telemetry;
pub mod validation;
pub mod vcf;
pub mod wgs_pipeline;

#[cfg(test)]
mod tests {
    use super::body::{HumanDigitalTwin, OrganSystem};
    use super::protocol::{HomeostasisSetpoints, NegligibleSenescenceController};
    use super::stress_test::run_clinical_stress_challenge_suite;
    use super::telemetry::evaluate_whole_body_six_axis_statistics;
    use super::validation::{
        run_monte_carlo_cohort_validation, validate_against_clinical_reference_ranges,
    };

    #[test]
    fn test_humanz_multi_organ_negligible_senescence_and_safety_guardrails() {
        let mut twin = HumanDigitalTwin::new("Subject-HZ-001", 72.0, 25, 2026);
        let controller = NegligibleSenescenceController::new(HomeostasisSetpoints::default());

        let summary = controller
            .run_negligible_senescence_lifespan_simulation(&mut twin, 15)
            .unwrap();

        // 1. Deep systemic age reversal achieved (72.0 yrs -> < 38.0 biological yrs, respecting Gompertz-Horvath floor)
        assert!(summary.post_reversal_snapshot.mean_biological_age_years < 38.0);

        // 2. Negligible senescence maintained across 15 chronological years (Age 72 -> 87, bio age <= 32.0)
        assert!(summary.final_treated_snapshot.mean_biological_age_years < 32.0);
        assert!(
            summary.net_aging_rate_treated_bio_yr_per_chrono_yr <= 0.10,
            "Treated net aging rate must be near zero (got {})",
            summary.net_aging_rate_treated_bio_yr_per_chrono_yr
        );
        assert!(summary.net_aging_rate_control_bio_yr_per_chrono_yr > 0.75);

        // 3. Post-senolytic & stem niche cellular viability >= 95% across all 6 organs
        assert!(summary.final_treated_snapshot.mean_viability >= 0.95);
        for organ in OrganSystem::all() {
            let org_stat = &summary.final_treated_snapshot.organs[&organ].stats;
            assert!(org_stat.mean_viability >= 0.95);
            assert!(org_stat.stem_cell_fraction <= 0.15 + 1e-9);
            assert!(summary.pbpk_precision_organ_doses.contains_key(&organ));
        }

        // 4. Whole-body safety guardrails approved across all 6 organs
        assert!(summary.safety_audit.all_organs_approved);
        assert!(summary.safety_audit.max_oncogenic_risk < 0.15);

        // 5. All 6 Mann-Whitney U axes pass p < 0.05
        let stats = evaluate_whole_body_six_axis_statistics(&twin).unwrap();
        assert!(stats.all_axes_significant());

        // 6. Clinical Stress & Acute Challenge Suite passes all resilience criteria for HumanZ
        let stress_report = run_clinical_stress_challenge_suite(&twin);
        assert!(stress_report.all_humanz_resilience_criteria_passed);
        assert!(
            stress_report.humanz_exercise.vo2_max_ml_kg_min
                > stress_report.control_exercise.vo2_max_ml_kg_min * 2.0
        );
        assert!(stress_report.humanz_cold.hypothermia_avoided);
        assert!(stress_report.humanz_endotoxemia.acute_kidney_injury_avoided);

        // 7. Clinical Reference-Range Validation (Framingham, CKD-EPI, ADA, ASE, ACSM)
        let ref_val = validate_against_clinical_reference_ranges(&twin);
        for check in &ref_val.checks {
            assert!(
                check.within_clinical_bounds,
                "Benchmark '{}' ({}) failed: observed {:.2} {} outside [{:.2}, {:.2}]",
                check.biomarker_name,
                check.cohort_group,
                check.observed_value,
                check.units,
                check.reference_min,
                check.reference_max
            );
        }
        assert!(ref_val.all_benchmarks_compliant);

        // 8. Multi-Seed Stochastic Monte Carlo Cohort & Sensitivity Analysis (N = 6 virtual subjects)
        let mc_report = run_monte_carlo_cohort_validation(6, 12).unwrap();
        assert_eq!(mc_report.cohort_size_n, 6);
        assert!(mc_report.fraction_passing_all_safety_guards >= 1.0);
        assert!(mc_report.fraction_passing_all_stress_tests >= 1.0);
        // Verify non-zero standard deviation across perturbed cohort subjects
        assert!(mc_report.distributions[0].control_sd > 0.1);
        assert!(mc_report.distributions[0].humanz_sd > 0.1);

        // 9. BioDb Genomic Personalization, PGx Dosing & Cardiometabolic-CKD Comorbidity Rescue
        let tmp_dir = std::env::temp_dir().join("humanz_test_export");
        let _ = std::fs::create_dir_all(&tmp_dir);
        let rescue = controller
            .run_personalized_comorbidity_rescue_simulation(&tmp_dir, 72.0, 10)
            .unwrap();
        assert_eq!(rescue.biodb_audit.chromosomes_indexed, 7);
        assert!(rescue.biodb_audit.all_guides_real_reference_provenance);
        assert!(rescue.biodb_audit.navitoclax_within_therapeutic_window);
        assert!(rescue.all_comorbidities_rescued);
        assert!(rescue.post_egfr_ml_min > rescue.pre_egfr_ml_min + 25.0);
        assert!(rescue.post_systolic_bp_mmhg < rescue.pre_systolic_bp_mmhg - 15.0);

        // 10. JSON & CSV Telemetry Export Verification
        let campaign = super::failure::run_failure_campaign().unwrap();
        let bundle = super::export::HumanzCompleteSimulationBundle {
            summary: summary.clone(),
            stress_challenge: stress_report,
            reference_validation: ref_val,
            monte_carlo_cohort: mc_report,
            failure_campaign: campaign.clone(),
            personalized_rescue: rescue,
        };
        let json_path = tmp_dir.join("humanz_bundle.json");
        let traj_csv = tmp_dir.join("trajectory.csv");
        let fail_csv = tmp_dir.join("failures.csv");
        super::export::export_simulation_bundle_json(&bundle, &json_path).unwrap();
        super::export::export_longitudinal_trajectory_csv(&summary, &traj_csv).unwrap();
        super::export::export_failure_campaign_csv(&campaign, &fail_csv).unwrap();
        assert!(json_path.exists() && std::fs::metadata(&json_path).unwrap().len() > 1000);
        assert!(traj_csv.exists() && std::fs::metadata(&traj_csv).unwrap().len() > 100);
        assert!(fail_csv.exists() && std::fs::metadata(&fail_csv).unwrap().len() > 100);
        let _ = std::fs::remove_dir_all(&tmp_dir);
    }
}
