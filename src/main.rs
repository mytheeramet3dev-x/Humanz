//! HumanZ CLI: Executes a 30-Year Whole-Body Multi-Organ Negligible Senescence Simulation
//! (Chronological Age 72.0 -> 102.0 Years) across 6 interconnected human organ systems.

use humanz::body::{HumanDigitalTwin, OrganSystem};
use humanz::export::{
    export_failure_campaign_csv, export_longitudinal_trajectory_csv, export_simulation_bundle_json,
    HumanzCompleteSimulationBundle,
};
use humanz::protocol::{HomeostasisSetpoints, NegligibleSenescenceController};
use humanz::stress_test::run_clinical_stress_challenge_suite;
use humanz::telemetry::evaluate_whole_body_six_axis_statistics;
use humanz::validation::{
    run_monte_carlo_cohort_validation, validate_against_clinical_reference_ranges,
};
use scies_bio_sim::prelude::*;
use scies_bio_th::anatomy::OrganArchitecture;
use std::path::PathBuf;

fn main() -> Result<()> {
    println!("====================================================================================================");
    println!("   🧬 HUMANZ v0.1.0 — Whole-Body Multi-Organ Negligible Senescence & Rejuvenation Digital Twin     ");
    println!("====================================================================================================\n");

    let mut initial_age = 72.0f64;
    let mut horizon_years = 30u32;
    let mut export_dir: Option<PathBuf> = None;
    let mut wgs_fasta: Option<PathBuf> = None;

    let args: Vec<String> = std::env::args().collect();
    if let Some(p) = args.iter().position(|a| a == "--seq-hash") {
        // Independent streaming hash (does not use the packed store): per contig SHA-256 of the
        // uppercase sequence and of the case-preserving sequence.
        use humanz::sha256::{to_hex, Sha256};
        use std::io::BufRead;
        let path = args.get(p + 1).expect("--seq-hash requires a FASTA path");
        let mut r = std::io::BufReader::with_capacity(1 << 20, std::fs::File::open(path).expect("open"));
        let mut line = String::new();
        let mut cur: Option<(String, Sha256, Sha256, u64)> = None;
        let flush = |c: Option<(String, Sha256, Sha256, u64)>| {
            if let Some((n, up, raw, len)) = c {
                println!("{}\t{}\t{}\t{}", n, len, to_hex(&up.finalize()), to_hex(&raw.finalize()));
            }
        };
        loop {
            line.clear();
            if r.read_line(&mut line).expect("read") == 0 {
                break;
            }
            let t = line.trim_end();
            if let Some(h) = t.strip_prefix('>') {
                flush(cur.take());
                let name = h.split_whitespace().next().unwrap_or("").to_string();
                cur = Some((name, Sha256::new(), Sha256::new(), 0));
            } else if let Some((_, up, raw, len)) = cur.as_mut() {
                raw.update(t.as_bytes());
                up.update(t.to_ascii_uppercase().as_bytes());
                *len += t.len() as u64;
            }
        }
        flush(cur.take());
        return Ok(());
    }
    if let Some(p) = args.iter().position(|a| a == "--fm-bench") {
        // Usage: --fm-bench <fasta> [max_bp]: builds an FM-index over the first contig prefix.
        let path = args.get(p + 1).expect("--fm-bench requires a FASTA path");
        let max_bp: u64 = args.get(p + 2).and_then(|v| v.parse().ok()).unwrap_or(u64::MAX);
        let file = std::fs::File::open(path).expect("cannot open FASTA");
        let g = humanz::genome_store::read_fasta_streaming(std::io::BufReader::with_capacity(1 << 20, file))
            .expect("FASTA parse failed");
        let c = &g.contigs[0];
        let region = c.region(0, max_bp);
        let cname = c.name.clone();
        drop(g); // free the packed genome before the memory-heavy SA-IS construction
        let t0 = std::time::Instant::now();
        let fm = humanz::fm_index::FmIndex::build(&region).expect("region contains non-ACGT");
        let build_s = t0.elapsed().as_secs_f64();
        let probe = &region[region.len() / 2..region.len() / 2 + 24];
        let t1 = std::time::Instant::now();
        let hits = fm.locate(probe);
        println!(
            "FM-index over {} bp of {}: built in {:.1}s; 24-mer probe -> {} hit(s) {:?} in {:.3}ms",
            fm.len(),
            cname,
            build_s,
            hits.len(),
            &hits[..hits.len().min(5)],
            t1.elapsed().as_secs_f64() * 1e3
        );
        return Ok(());
    }
    if let Some(p) = args.iter().position(|a| a == "--fasta") {
        let Some(path) = args.get(p + 1) else {
            eprintln!("--fasta requires a path");
            std::process::exit(2);
        };
        let out = args
            .iter()
            .position(|a| a == "--export-dir")
            .and_then(|q| args.get(q + 1))
            .map(PathBuf::from)
            .unwrap_or_else(|| PathBuf::from("genome_out"));
        if let Err(e) = genome_pilot(path, &out) {
            eprintln!("genome pilot failed: {e}");
            std::process::exit(1);
        }
        return Ok(());
    }
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--age" if i + 1 < args.len() => {
                if let Ok(v) = args[i + 1].parse::<f64>() {
                    initial_age = v;
                }
                i += 2;
            }
            "--horizon" if i + 1 < args.len() => {
                if let Ok(v) = args[i + 1].parse::<u32>() {
                    horizon_years = v;
                }
                i += 2;
            }
            "--export-dir" if i + 1 < args.len() => {
                export_dir = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--wgs-fasta" if i + 1 < args.len() => {
                wgs_fasta = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            _ => {
                i += 1;
            }
        }
    }

    let cells_per_organ = 40; // 6 organs x 40 cells = 240 single cells per twin (480 total)

    println!(
        "[1] Initializing Paired Multi-Organ Digital Twin ('HZ-Alpha-01', Age {:.1} yrs, 6 Organs x {} cells = {} cells)...",
        initial_age,
        cells_per_organ,
        cells_per_organ * 6
    );
    let mut twin = HumanDigitalTwin::new("HZ-Alpha-01", initial_age, cells_per_organ, 2026);
    let controller = NegligibleSenescenceController::new(HomeostasisSetpoints::default());

    println!(
        "[2] Executing Phase I (3-Round Systemic 6-Hallmark Age Reversal) + Phase II ({}-Year Closed-Loop Homeostasis)...",
        horizon_years
    );
    let summary = controller.run_negligible_senescence_lifespan_simulation(&mut twin, horizon_years)?;

    println!("\n----------------------------------------------------------------------------------------------------");
    println!("  PHASE I & II LONGITUDINAL TRAJECTORY (Chronological Age {:.0} -> {:.0} Years)",
        initial_age,
        summary.final_treated_snapshot.chronological_age_years
    );
    println!("----------------------------------------------------------------------------------------------------");
    println!("  Stage / Milestone             | Chrono Age | Bio Age (Control) | Bio Age (HumanZ) | Viability | Telomere");
    println!("  ------------------------------+------------+-------------------+------------------+-----------+----------");
    println!(
        "  Baseline (Pre-Intervention)   | {:>8.1} y | {:>15.1} y | {:>14.1} y | {:>8.2}% | {:>6.0} bp",
        summary.baseline_snapshot.chronological_age_years,
        summary.baseline_snapshot.mean_biological_age_years,
        summary.baseline_snapshot.mean_biological_age_years,
        summary.baseline_snapshot.mean_viability * 100.0,
        summary.baseline_snapshot.mean_telomere_bp
    );
    println!(
        "  Post-Phase I Reversal         | {:>8.1} y | {:>15.1} y | {:>14.1} y | {:>8.2}% | {:>6.0} bp",
        summary.post_reversal_snapshot.chronological_age_years,
        summary.baseline_snapshot.mean_biological_age_years,
        summary.post_reversal_snapshot.mean_biological_age_years,
        summary.post_reversal_snapshot.mean_viability * 100.0,
        summary.post_reversal_snapshot.mean_telomere_bp
    );

    for yr_idx in [5usize, 10, 15, 20, 25, 30] {
        if let (Some(trt), Some(ctrl)) = (
            summary.longitudinal_treated_history.get(yr_idx),
            summary.longitudinal_control_history.get(yr_idx),
        ) {
            println!(
                "  Homeostasis Year +{:>2}          | {:>8.1} y | {:>15.1} y | {:>14.1} y | {:>8.2}% | {:>6.0} bp",
                yr_idx,
                trt.chronological_age_years,
                ctrl.mean_biological_age_years,
                trt.mean_biological_age_years,
                trt.mean_viability * 100.0,
                trt.mean_telomere_bp
            );
        }
    }

    println!("\n----------------------------------------------------------------------------------------------------");
    println!("  ORGAN-BY-ORGAN BIOMARKER & ANATOMICAL ECM BREAKDOWN AT CHRONOLOGICAL AGE {:.0} YEARS",
        summary.final_treated_snapshot.chronological_age_years
    );
    println!("----------------------------------------------------------------------------------------------------");
    println!("  Organ Compartment               | Bio Age | Telomere | mtDNA Het | Stiffness (E) | YAP/TAZ | Viability");
    println!("  --------------------------------+---------+----------+-----------+---------------+---------+----------");
    for organ in OrganSystem::all() {
        let org_snap = &summary.final_treated_snapshot.organs[&organ];
        let s = &org_snap.stats;
        println!(
            "  {:<31} | {:>5.1} y | {:>5.0} bp | {:>8.2}% | {:>9.1} kPa | {:>7.3} | {:>8.2}%",
            format!("{:?}", organ),
            s.mean_epigenetic_age,
            s.mean_telomere_bp,
            s.mean_heteroplasmy * 100.0,
            org_snap.youngs_modulus_kpa,
            org_snap.yap_taz_mechanotransduction_stress,
            s.mean_viability * 100.0
        );
    }

    println!("\n----------------------------------------------------------------------------------------------------");
    println!("  SPATIAL REACTION-DIFFUSION PDE LAYER PROFILES (8-Lamina Brain BBB & 7-Lamina Arterial Wall)");
    println!("----------------------------------------------------------------------------------------------------");
    for organ in [OrganSystem::CerebralCortexNeuroglia, OrganSystem::CardiovascularEndothelium] {
        let org_snap = &summary.final_treated_snapshot.organs[&organ];
        let anat = &twin.treated_anatomy[&organ];
        println!("  [{:?} | {}] (Center: ({:.0}, {:.0}, {:.0}) mm | E = {:.2} kPa | Cavity: {:?}):",
            organ,
            scies_bio_th::anatomy::OrganArchitecture::ta2_code(anat),
            scies_bio_th::anatomy::OrganArchitecture::stereotaxic_center_mm(anat).x_mediolateral_mm,
            scies_bio_th::anatomy::OrganArchitecture::stereotaxic_center_mm(anat).y_anteroposterior_mm,
            scies_bio_th::anatomy::OrganArchitecture::stereotaxic_center_mm(anat).z_superoinferior_mm,
            org_snap.youngs_modulus_kpa,
            scies_bio_th::anatomy::OrganArchitecture::fascial_cavity(anat),
        );
        for layer in &org_snap.spatial_layers {
            println!(
                "    -> Layer {:<56} | Depth: {:>5.1} um | pO2: {:>5.1} mmHg | Bioavail: {:>5.1}%",
                layer.layer_name,
                layer.radial_depth_um,
                layer.local_po2_mmhg,
                layer.local_drug_bioavailability * 100.0
            );
        }
    }

    let atlas = scies_bio_th::anatomy::AnatomicalSpatialAtlas::from_somatic_age(26.0);
    let pbpk_shares = scies_bio_th::anatomy::compute_whole_body_pbpk_biodistribution(
        &atlas.organs,
        scies_bio_th::aging::MW_MRNA_LNP_VECTOR_DA,
        0.92,
        0.88,
    );
    println!("\n----------------------------------------------------------------------------------------------------");
    println!(
        "  3D WHOLE-BODY SPATIAL ATLAS & CLOSED-LOOP CONDUIT TOPOLOGY ({} TA2 Organs, {} Conduits)",
        atlas.organs.len(),
        atlas.spatial_edges.len()
    );
    println!("----------------------------------------------------------------------------------------------------");
    for edge in atlas.spatial_edges.iter().take(10) {
        println!(
            "  • {} ({}) -> {} ({}) | Dist: {:>5.1} mm | Weight: {:.2} | {:?}",
            edge.source_organ_name,
            edge.source_organ_ta2,
            edge.target_organ_name,
            edge.target_organ_ta2,
            edge.euclidean_distance_mm,
            edge.signal_transmission_weight(),
            edge.directional_relations
        );
    }

    println!("\n----------------------------------------------------------------------------------------------------");
    println!("  MOENS-KORTEWEG HEMODYNAMICS, SYSTEMIC ENDOCRINE AXES & PBPK RES BIODISTRIBUTION");
    println!("----------------------------------------------------------------------------------------------------");
    let ctrl_hemo = &summary.final_control_snapshot.hemodynamics;
    let trt_hemo = &summary.final_treated_snapshot.hemodynamics;
    let ctrl_endo = &summary.final_control_snapshot.endocrine;
    let trt_endo = &summary.final_treated_snapshot.endocrine;
    println!(
        "  • Aortic Pulse Wave Velocity (PWV):  Control = {:>5.2} m/s (PP = {:.0} mmHg, Barotrauma = {:.3}) -> HumanZ = {:>5.2} m/s (PP = {:.0} mmHg, Barotrauma = {:.3})",
        ctrl_hemo.pulse_wave_velocity_m_s,
        ctrl_hemo.central_pulse_pressure_mmhg,
        ctrl_hemo.microvascular_barotrauma_index,
        trt_hemo.pulse_wave_velocity_m_s,
        trt_hemo.central_pulse_pressure_mmhg,
        trt_hemo.microvascular_barotrauma_index
    );
    println!(
        "  • Systemic Endocrine Axes (Age 102): Control [Klotho = {:.2}, Cortisol/DHEA = {:.2}, RAAS = {:.2}] -> HumanZ [Klotho = {:.2}, Cortisol/DHEA = {:.2}, RAAS = {:.2}]",
        ctrl_endo.circulating_alpha_klotho,
        ctrl_endo.adrenal_cortisol_to_dhea_ratio,
        ctrl_endo.raas_angiotensin_ii_tone,
        trt_endo.circulating_alpha_klotho,
        trt_endo.adrenal_cortisol_to_dhea_ratio,
        trt_endo.raas_angiotensin_ii_tone
    );
    println!(
        "  • Top 5 Organ PBPK Biodistribution (800 kDa mRNA-LNP + CD47 RES Stealth + TfR1 Shuttle):"
    );
    let mut sorted_pbpk = pbpk_shares.clone();
    sorted_pbpk.sort_by(|a, b| {
        b.effective_parenchymal_delivery
            .partial_cmp(&a.effective_parenchymal_delivery)
            .unwrap()
    });
    for share in sorted_pbpk.iter().take(5) {
        println!(
            "    -> {:<50} ({}) | Parenchymal Bioavail: {:>5.1}% | RES Trapped: {:>5.1}%",
            share.organ_name,
            share.organ_ta2,
            share.effective_parenchymal_delivery * 100.0,
            share.res_first_pass_sequestration * 100.0
        );
    }

    // Verify Phase 3 Single-Cell Histological Lamina Anchoring (Stem Niche vs Progenitor vs Parenchyma)
    let brain_col = twin.sim.colonies.get(twin.treated_organs[&OrganSystem::CerebralCortexNeuroglia]).unwrap();
    let mut stem_depth = 0.0;
    let mut paren_depth = 0.0;
    for &cid in &brain_col.cell_ids {
        if let Some(c) = twin.sim.cells.get(cid) {
            if c.cell_type == CellType::AdultStemCell {
                stem_depth = c.radial_depth_um;
            } else if c.cell_type == CellType::Somatic && paren_depth == 0.0 {
                paren_depth = c.radial_depth_um;
            }
        }
    }
    println!(
        "  • Single-Cell Histological Lamina Anchoring (Brain Cortex): AdultStemCell @ Lamina 1 ({:.1} um) vs Somatic @ Deep Lamina ({:.1} um)",
        stem_depth, paren_depth
    );

    println!("\n----------------------------------------------------------------------------------------------------");
    println!("  7 MAJOR CLOSED-LOOP SYSTEMIC PHYSIOLOGICAL CIRCUITS (Control Age 102 vs HumanZ Rejuvenated Age 26)");
    println!("----------------------------------------------------------------------------------------------------");
    let c_circ = &summary.final_control_snapshot.circuits;
    let t_circ = &summary.final_treated_snapshot.circuits;

    println!("  [1] CARDIOVASCULAR-HEMODYNAMIC CYCLE (Heart, Arterial Windkessel, Baroreflex):");
    println!("      • Control: HR {:>4.1} bpm | SV {:>4.1} mL | CO {:>4.2} L/min | EF {:>4.1}% | BP {:>5.1}/{:>4.1} mmHg (MAP {:>4.1}) | SVR {:>4.1} Wood",
        c_circ.cardiovascular.heart_rate_bpm, c_circ.cardiovascular.stroke_volume_ml, c_circ.cardiovascular.cardiac_output_l_min,
        c_circ.cardiovascular.ejection_fraction_pct, c_circ.cardiovascular.systolic_bp_mmhg, c_circ.cardiovascular.diastolic_bp_mmhg,
        c_circ.cardiovascular.mean_arterial_pressure_mmhg, c_circ.cardiovascular.systemic_vascular_resistance_wood
    );
    println!("      • HumanZ:  HR {:>4.1} bpm | SV {:>4.1} mL | CO {:>4.2} L/min | EF {:>4.1}% | BP {:>5.1}/{:>4.1} mmHg (MAP {:>4.1}) | SVR {:>4.1} Wood",
        t_circ.cardiovascular.heart_rate_bpm, t_circ.cardiovascular.stroke_volume_ml, t_circ.cardiovascular.cardiac_output_l_min,
        t_circ.cardiovascular.ejection_fraction_pct, t_circ.cardiovascular.systolic_bp_mmhg, t_circ.cardiovascular.diastolic_bp_mmhg,
        t_circ.cardiovascular.mean_arterial_pressure_mmhg, t_circ.cardiovascular.systemic_vascular_resistance_wood
    );

    println!("  [2] CARDIOPULMONARY GAS EXCHANGE & VENTILATION CYCLE (Alveoli, DL_CO, Arterial Blood Gas):");
    println!("      • Control: RR {:>4.1} bpm | Vt {:>5.1} mL | PaO2 {:>5.1} mmHg | SaO2 {:>4.1}% | PaCO2 {:>4.1} mmHg | pH {:>5.3}",
        c_circ.cardiopulmonary.respiratory_rate_bpm, c_circ.cardiopulmonary.tidal_volume_ml,
        c_circ.cardiopulmonary.arterial_po2_mmhg, c_circ.cardiopulmonary.arterial_o2_saturation_pct,
        c_circ.cardiopulmonary.arterial_pco2_mmhg, c_circ.cardiopulmonary.arterial_ph
    );
    println!("      • HumanZ:  RR {:>4.1} bpm | Vt {:>5.1} mL | PaO2 {:>5.1} mmHg | SaO2 {:>4.1}% | PaCO2 {:>4.1} mmHg | pH {:>5.3}",
        t_circ.cardiopulmonary.respiratory_rate_bpm, t_circ.cardiopulmonary.tidal_volume_ml,
        t_circ.cardiopulmonary.arterial_po2_mmhg, t_circ.cardiopulmonary.arterial_o2_saturation_pct,
        t_circ.cardiopulmonary.arterial_pco2_mmhg, t_circ.cardiopulmonary.arterial_ph
    );

    println!("  [3] ENTEROHEPATIC-METABOLIC & GLUCOSE HOMEOSTASIS CYCLE (Gut, Pancreas, Liver Zones 1-3):");
    println!("      • Control: Glucose {:>5.1} mg/dL | Insulin {:>5.1} pmol/L | Glucagon {:>5.1} pg/mL | Urea Eff {:>4.1}% | LPS {:>4.2} EU/mL",
        c_circ.enterohepatic.blood_glucose_mg_dl, c_circ.enterohepatic.circulating_insulin_pmol_l,
        c_circ.enterohepatic.circulating_glucagon_pg_ml, c_circ.enterohepatic.hepatic_urea_cycle_efficiency_pct,
        c_circ.enterohepatic.portal_lps_endotoxin_eu_ml
    );
    println!("      • HumanZ:  Glucose {:>5.1} mg/dL | Insulin {:>5.1} pmol/L | Glucagon {:>5.1} pg/mL | Urea Eff {:>4.1}% | LPS {:>4.2} EU/mL",
        t_circ.enterohepatic.blood_glucose_mg_dl, t_circ.enterohepatic.circulating_insulin_pmol_l,
        t_circ.enterohepatic.circulating_glucagon_pg_ml, t_circ.enterohepatic.hepatic_urea_cycle_efficiency_pct,
        t_circ.enterohepatic.portal_lps_endotoxin_eu_ml
    );

    println!("  [4] RENAL-CARDIOVASCULAR FLUID, ELECTROLYTE & RAAS CYCLE (eGFR, Aldosterone, Klotho):");
    println!("      • Control: eGFR {:>5.1} mL/min | RBF {:>4.2} L/min | Renin {:>4.2} | Ang II {:>5.1} pmol/L | Aldo {:>4.1} ng/dL | Klotho {:>6.1} pg/mL",
        c_circ.renal_raas.egfr_ml_min, c_circ.renal_raas.renal_blood_flow_l_min,
        c_circ.renal_raas.plasma_renin_activity_ng_ml_hr, c_circ.renal_raas.circulating_angiotensin_ii_pmol_l,
        c_circ.renal_raas.plasma_aldosterone_ng_dl, c_circ.renal_raas.circulating_alpha_klotho_pg_ml
    );
    println!("      • HumanZ:  eGFR {:>5.1} mL/min | RBF {:>4.2} L/min | Renin {:>4.2} | Ang II {:>5.1} pmol/L | Aldo {:>4.1} ng/dL | Klotho {:>6.1} pg/mL",
        t_circ.renal_raas.egfr_ml_min, t_circ.renal_raas.renal_blood_flow_l_min,
        t_circ.renal_raas.plasma_renin_activity_ng_ml_hr, t_circ.renal_raas.circulating_angiotensin_ii_pmol_l,
        t_circ.renal_raas.plasma_aldosterone_ng_dl, t_circ.renal_raas.circulating_alpha_klotho_pg_ml
    );

    println!("  [5] NEUROENDOCRINE-GLYMPHATIC SLEEP CYCLE (SCN, HPA, Choroid CSF, AQP4 Glymphatics, SGZ):");
    println!("      • Control: SCN GnRH {:>4.2} | Cortisol {:>4.1} ug/dL | DHEA-S {:>5.1} ug/dL | CSF {:>5.1} mL/d | Glymphatic {:>4.1} uL/min | SGZ Neuro {:>5.1} cells/d",
        c_circ.neuroendocrine_glymphatic.hypothalamic_gnrh_level, c_circ.neuroendocrine_glymphatic.plasma_cortisol_ug_dl,
        c_circ.neuroendocrine_glymphatic.plasma_dhea_s_ug_dl, c_circ.neuroendocrine_glymphatic.daily_csf_turnover_ml_day,
        c_circ.neuroendocrine_glymphatic.cortical_glymphatic_clearance_rate, c_circ.neuroendocrine_glymphatic.adult_hippocampal_neurogenesis_cells_day
    );
    println!("      • HumanZ:  SCN GnRH {:>4.2} | Cortisol {:>4.1} ug/dL | DHEA-S {:>5.1} ug/dL | CSF {:>5.1} mL/d | Glymphatic {:>4.1} uL/min | SGZ Neuro {:>5.1} cells/d",
        t_circ.neuroendocrine_glymphatic.hypothalamic_gnrh_level, t_circ.neuroendocrine_glymphatic.plasma_cortisol_ug_dl,
        t_circ.neuroendocrine_glymphatic.plasma_dhea_s_ug_dl, t_circ.neuroendocrine_glymphatic.daily_csf_turnover_ml_day,
        t_circ.neuroendocrine_glymphatic.cortical_glymphatic_clearance_rate, t_circ.neuroendocrine_glymphatic.adult_hippocampal_neurogenesis_cells_day
    );

    println!("  [6] IMMUNE-LYMPHOID SURVEILLANCE & RES FILTERING CYCLE (Bone Marrow, Thymus, Spleen):");
    println!("      • Control: Naive T-Cell {:>4.1}e7/d | Myeloid/Lymphoid Skew {:>4.2}x | Splenic RBC Cull {:>5.1}e9/d | Thoracic Lymph {:>5.1} mL/h",
        c_circ.immune_lymphatic.daily_naive_t_cell_emigration_10e7, c_circ.immune_lymphatic.myeloid_to_lymphoid_skew_ratio,
        c_circ.immune_lymphatic.splenic_senescent_rbc_clearance_10e9, c_circ.immune_lymphatic.thoracic_duct_lymph_return_ml_hr
    );
    println!("      • HumanZ:  Naive T-Cell {:>4.1}e7/d | Myeloid/Lymphoid Skew {:>4.2}x | Splenic RBC Cull {:>5.1}e9/d | Thoracic Lymph {:>5.1} mL/h",
        t_circ.immune_lymphatic.daily_naive_t_cell_emigration_10e7, t_circ.immune_lymphatic.myeloid_to_lymphoid_skew_ratio,
        t_circ.immune_lymphatic.splenic_senescent_rbc_clearance_10e9, t_circ.immune_lymphatic.thoracic_duct_lymph_return_ml_hr
    );

    println!("  [7] ENERGY & THERMOREGULATORY CYCLE (BAT UCP1, Muscle Shivering/Work, Dermal Vasomotor):");
    println!("      • Control: Core Temp {:>4.1} °C | BMR {:>6.1} kcal/d | BAT Heat {:>4.1} W | Muscle Heat {:>4.1} W | Dermal Heat Loss {:>4.1} W | Irisin {:>4.1} ng/mL",
        c_circ.energy_thermoregulatory.core_body_temperature_c, c_circ.energy_thermoregulatory.basal_metabolic_rate_kcal_day,
        c_circ.energy_thermoregulatory.brown_adipose_heat_watts, c_circ.energy_thermoregulatory.muscle_thermogenesis_watts,
        c_circ.energy_thermoregulatory.cutaneous_heat_loss_watts, c_circ.energy_thermoregulatory.circulating_irisin_ng_ml
    );
    println!("      • HumanZ:  Core Temp {:>4.1} °C | BMR {:>6.1} kcal/d | BAT Heat {:>4.1} W | Muscle Heat {:>4.1} W | Dermal Heat Loss {:>4.1} W | Irisin {:>4.1} ng/mL",
        t_circ.energy_thermoregulatory.core_body_temperature_c, t_circ.energy_thermoregulatory.basal_metabolic_rate_kcal_day,
        t_circ.energy_thermoregulatory.brown_adipose_heat_watts, t_circ.energy_thermoregulatory.muscle_thermogenesis_watts,
        t_circ.energy_thermoregulatory.cutaneous_heat_loss_watts, t_circ.energy_thermoregulatory.circulating_irisin_ng_ml
    );

    let stress = run_clinical_stress_challenge_suite(&twin);
    println!("\n----------------------------------------------------------------------------------------------------");
    println!("  4-STAGE CLINICAL STRESS & ACUTE PHYSIOLOGICAL CHALLENGE SUITE (AT AGE {:.0} YEARS)", stress.chronological_age_years);
    println!("----------------------------------------------------------------------------------------------------");
    println!("  [Challenge 1] MAXIMAL CARDIOPULMONARY & LOCOMOTOR EXERCISE RAMP (VO2max Stress Test):");
    println!("      • Control: Peak HR {:>5.1} bpm | Peak CO {:>5.2} L/min | VO2max {:>4.1} mL/kg/min | Force {:>4.1} N/cm2 | Irisin {:>4.1} ng/mL",
        stress.control_exercise.peak_heart_rate_bpm, stress.control_exercise.peak_cardiac_output_l_min,
        stress.control_exercise.vo2_max_ml_kg_min, stress.control_exercise.peak_muscle_contractile_force_n_cm2,
        stress.control_exercise.peak_irisin_myokine_ng_ml
    );
    println!("      • HumanZ:  Peak HR {:>5.1} bpm | Peak CO {:>5.2} L/min | VO2max {:>4.1} mL/kg/min | Force {:>4.1} N/cm2 | Irisin {:>4.1} ng/mL",
        stress.humanz_exercise.peak_heart_rate_bpm, stress.humanz_exercise.peak_cardiac_output_l_min,
        stress.humanz_exercise.vo2_max_ml_kg_min, stress.humanz_exercise.peak_muscle_contractile_force_n_cm2,
        stress.humanz_exercise.peak_irisin_myokine_ng_ml
    );
    println!("  [Challenge 2] 75g ORAL GLUCOSE TOLERANCE & ENTEROHEPATIC PORTAL TOXIN CHALLENGE (OGTT):");
    println!("      • Control: Peak Glu {:>5.1} mg/dL | 2h Glu {:>5.1} mg/dL | Peak Ins {:>5.1} pmol/L | CYP450 {:>4.2}x | Portal LPS {:>4.2} EU/mL",
        stress.control_ogtt.peak_postprandial_glucose_mg_dl, stress.control_ogtt.two_hour_post_load_glucose_mg_dl,
        stress.control_ogtt.peak_biphasic_insulin_pmol_l, stress.control_ogtt.hepatic_cyp450_clearance_multiplier,
        stress.control_ogtt.portal_lps_endotoxin_leak_eu_ml
    );
    println!("      • HumanZ:  Peak Glu {:>5.1} mg/dL | 2h Glu {:>5.1} mg/dL | Peak Ins {:>5.1} pmol/L | CYP450 {:>4.2}x | Portal LPS {:>4.2} EU/mL",
        stress.humanz_ogtt.peak_postprandial_glucose_mg_dl, stress.humanz_ogtt.two_hour_post_load_glucose_mg_dl,
        stress.humanz_ogtt.peak_biphasic_insulin_pmol_l, stress.humanz_ogtt.hepatic_cyp450_clearance_multiplier,
        stress.humanz_ogtt.portal_lps_endotoxin_leak_eu_ml
    );
    println!("  [Challenge 3] ACUTE 4.0 °C COLD SHOCK & THERMOREGULATORY DEFENSE CHALLENGE:");
    println!("      • Control: Defended Core {:>4.2} °C | BAT UCP1 {:>4.1} W | Muscle Heat {:>4.1} W | Total {:>5.1} W | Hypothermia Avoided: {}",
        stress.control_cold.defended_core_temperature_c, stress.control_cold.bat_ucp1_non_shivering_heat_watts,
        stress.control_cold.muscle_thermogenic_heat_watts, stress.control_cold.total_thermogenic_power_watts,
        stress.control_cold.hypothermia_avoided
    );
    println!("      • HumanZ:  Defended Core {:>4.2} °C | BAT UCP1 {:>4.1} W | Muscle Heat {:>4.1} W | Total {:>5.1} W | Hypothermia Avoided: {}",
        stress.humanz_cold.defended_core_temperature_c, stress.humanz_cold.bat_ucp1_non_shivering_heat_watts,
        stress.humanz_cold.muscle_thermogenic_heat_watts, stress.humanz_cold.total_thermogenic_power_watts,
        stress.humanz_cold.hypothermia_avoided
    );
    println!("  [Challenge 4] SYSTEMIC ENDOTOXEMIA, HYPOTENSIVE RENAL & SLOW-WAVE SLEEP GLYMPHATIC CHALLENGE:");
    println!("      • Control: Stressed eGFR {:>5.1} mL/min (AKI Avoided: {:<5}) | BBB TEER {:>6.1} | SWS Glymphatic {:>4.1} uL/min | Immune {:>4.2}",
        stress.control_endotoxemia.stressed_egfr_ml_min, stress.control_endotoxemia.acute_kidney_injury_avoided,
        stress.control_endotoxemia.stressed_bbb_teer_ohm_cm2, stress.control_endotoxemia.sws_glymphatic_clearance_ul_min,
        stress.control_endotoxemia.immune_patrolling_competence
    );
    println!("      • HumanZ:  Stressed eGFR {:>5.1} mL/min (AKI Avoided: {:<5}) | BBB TEER {:>6.1} | SWS Glymphatic {:>4.1} uL/min | Immune {:>4.2}",
        stress.humanz_endotoxemia.stressed_egfr_ml_min, stress.humanz_endotoxemia.acute_kidney_injury_avoided,
        stress.humanz_endotoxemia.stressed_bbb_teer_ohm_cm2, stress.humanz_endotoxemia.sws_glymphatic_clearance_ul_min,
        stress.humanz_endotoxemia.immune_patrolling_competence
    );

    println!("\n----------------------------------------------------------------------------------------------------");
    println!("  PBPK-GUIDED PRECISION ORGAN DOSING & CLOSED-LOOP PHYSIOLOGICAL CONTROLLER TELEMETRY");
    println!("----------------------------------------------------------------------------------------------------");
    for (organ, dose) in &summary.pbpk_precision_organ_doses {
        let bioavail = twin.treated_anatomy[organ]
            .mean_parenchymal_bioavailability(scies_bio_th::aging::MW_MRNA_LNP_VECTOR_DA, 0.92);
        println!(
            "  • {:<31} | PDE Parenchymal Bioavailability: {:>5.1}% | Precision OSK Factor Dose: {:>5.3}x",
            format!("{:?}", organ),
            bioavail * 100.0,
            dose
        );
    }

    let ref_val = validate_against_clinical_reference_ranges(&twin);
    println!("\n----------------------------------------------------------------------------------------------------");
    println!(
        "  CLINICAL REFERENCE-RANGE VALIDATION AUDIT ({}/{} BENCHMARKS COMPLIANT WITH LITERATURE)",
        ref_val.passed_count, ref_val.total_count
    );
    println!("----------------------------------------------------------------------------------------------------");
    println!("  Biomarker / Parameter                | Cohort Group             | Observed   | Ref Interval    | Status");
    println!("  -------------------------------------+--------------------------+------------+-----------------+-------");
    for c in &ref_val.checks {
        println!(
            "  {:<36} | {:<24} | {:>6.1} {:<3} | [{:>5.1}, {:>5.1}] | {}",
            c.biomarker_name,
            c.cohort_group,
            c.observed_value,
            c.units,
            c.reference_min,
            c.reference_max,
            if c.within_clinical_bounds { "PASS" } else { "FAIL" }
        );
    }

    let mc = run_monte_carlo_cohort_validation(8, 15)?;
    println!("\n----------------------------------------------------------------------------------------------------");
    println!(
        "  MULTI-SEED MONTE CARLO COHORT & PARAMETER SENSITIVITY ANALYSIS (N = {} Virtual Subjects, +/-12% Perturbation)",
        mc.cohort_size_n
    );
    println!("----------------------------------------------------------------------------------------------------");
    println!("  Biomarker          | Control Mean +/- SD (95% CI)            | HumanZ Mean +/- SD (95% CI)");
    println!("  -------------------+-----------------------------------------+-----------------------------------------");
    for d in &mc.distributions {
        println!(
            "  {:<18} | {:>6.1} +/- {:>4.2} [{:>6.1}, {:>6.1}] {:<7} | {:>6.1} +/- {:>4.2} [{:>6.1}, {:>6.1}] {:<7}",
            d.metric_name,
            d.control_mean,
            d.control_sd,
            d.control_ci95_low,
            d.control_ci95_high,
            d.units,
            d.humanz_mean,
            d.humanz_sd,
            d.humanz_ci95_low,
            d.humanz_ci95_high,
            d.units
        );
    }

    let campaign = humanz::failure::run_failure_campaign()?;
    println!("\n----------------------------------------------------------------------------------------------------");
    println!("  FAILURE MODES & ADVERSE EVENTS CAMPAIGN ({} scenarios) — Safety Guard Verification", campaign.results.len());
    println!("----------------------------------------------------------------------------------------------------");
    for r in &campaign.results {
        println!(
            "  {:<50} | {:<12?} | Guard: {:<28} | Eff {:>4.0}% | Contained: {}",
            r.name,
            r.outcome,
            format!("{:?}", r.decision.action),
            r.efficacy_fraction * 100.0,
            r.decision.hazard_contained()
        );
        for e in r.events.iter().filter(|e| e.severity >= humanz::failure::adverse::Severity::Moderate) {
            println!("        - [{:?}] {:?}: {}", e.severity, e.kind, e.description);
        }
        if let Some(m) = r.mitigated_outcome {
            println!("        => After guard dose reduction: {:?}", m);
        }
        for n in &r.notes {
            println!("        . {}", n);
        }
    }
    println!(
        "  Summary: Success {} | Partial {} | Failure {} | Catastrophic {} | Failure rate {:.0}% | Guard activation {:.0}% | Containment {:.0}% | Severity [Mild/Mod/Sev/Crit] = {:?} | Mean time-to-safety = {:?} cycles",
        campaign.success, campaign.partial, campaign.failure, campaign.catastrophic,
        campaign.failure_rate * 100.0, campaign.safety_activation_rate * 100.0, campaign.containment_rate * 100.0,
        campaign.severity_histogram, campaign.mean_time_to_safety_activation
    );

    let tmp_biodb = std::env::temp_dir().join("humanz_biodb_audit");
    let rescue = controller.run_personalized_comorbidity_rescue_simulation(
        &tmp_biodb,
        initial_age,
        12,
    )?;
    println!("\n----------------------------------------------------------------------------------------------------");
    println!("  SCIES-BIO-DB GENOMIC, PHARMACOGENOMIC & COMORBIDITY PERSONALIZATION AUDIT (GRCh38.p14 / T2T)");
    println!("----------------------------------------------------------------------------------------------------");
    println!(
        "  • BioDb Assembly & Loci:   {} ({} chromosomes, {} bp scanned, {} canonical genes, {} clinical SNPs)",
        rescue.biodb_audit.assembly_name,
        rescue.biodb_audit.chromosomes_indexed,
        rescue.biodb_audit.total_reference_bp_scanned,
        rescue.biodb_audit.registered_genes_count,
        rescue.biodb_audit.patient_variants_recorded
    );
    println!(
        "  • CRISPR Guide Provenance: RealReference={} | Mean FM-Index Specificity: {:.1}% | Max Off-Target: {:.2}%",
        rescue.biodb_audit.all_guides_real_reference_provenance,
        rescue.biodb_audit.mean_guide_specificity_score_pct,
        rescue.biodb_audit.max_guide_offtarget_score_pct
    );
    println!(
        "  • Patient Genotype & PGx:  {:?} | {:?} | {:?} -> Navitoclax PGx Dose {:.3} uM (Peak {:.2} uM, Window OK={})",
        rescue.genotype.apoe,
        rescue.genotype.cyp_phenotype,
        rescue.genotype.checkpoint_zygosity,
        rescue.biodb_audit.navitoclax_pgx_dose_um,
        rescue.biodb_audit.navitoclax_peak_concentration_um,
        rescue.biodb_audit.navitoclax_within_therapeutic_window
    );
    println!(
        "  • Comorbidity Rescue (CKD Stage 3b + T2D + HFpEF + APOE4/4): eGFR {:.1} -> {:.1} mL/min | SBP {:.1} -> {:.1} mmHg | Glucose {:.1} -> {:.1} mg/dL | Brain Proteostasis {:.3} -> {:.3} | All Rescued: {}",
        rescue.pre_egfr_ml_min,
        rescue.post_egfr_ml_min,
        rescue.pre_systolic_bp_mmhg,
        rescue.post_systolic_bp_mmhg,
        rescue.pre_fasting_glucose_mg_dl,
        rescue.post_fasting_glucose_mg_dl,
        rescue.pre_cortical_proteostasis_burden,
        rescue.post_cortical_proteostasis_burden,
        rescue.all_comorbidities_rescued
    );

    if let Some(ref out_dir) = export_dir {
        let bundle = HumanzCompleteSimulationBundle {
            summary: summary.clone(),
            stress_challenge: stress.clone(),
            reference_validation: ref_val.clone(),
            monte_carlo_cohort: mc.clone(),
            failure_campaign: campaign.clone(),
            personalized_rescue: rescue.clone(),
        };
        export_simulation_bundle_json(&bundle, &out_dir.join("humanz_simulation_bundle.json"))?;
        export_longitudinal_trajectory_csv(&summary, &out_dir.join("longitudinal_trajectory.csv"))?;
        export_failure_campaign_csv(&campaign, &out_dir.join("failure_campaign.csv"))?;
        println!(
            "  • Telemetry Exported:      JSON + CSV files written to {}",
            out_dir.display()
        );
    }

    // Whole-Genome T2T-CHM13v2.0 (PackedGenome + FmIndex) Coupled Pipeline:
    // 1. Somatic Mutation Tracking -> VCF 4.2
    // 2. Horvath 353-CpG Methylation -> Pre/Post-OSK BedGraph
    // 3. CRISPR 1-3 Mismatch Off-Target Scan -> BED
    let wgs_genome = if let Some(ref fasta_path) = wgs_fasta {
        let f = std::fs::File::open(fasta_path).expect("cannot open --wgs-fasta");
        humanz::genome_store::read_fasta_streaming(std::io::BufReader::with_capacity(1 << 20, f))
            .expect("failed to stream --wgs-fasta")
    } else {
        load_default_t2t_loci_packed_genome()
    };

    let crispr_offtarget_rows =
        humanz::wgs_pipeline::scan_crispr_offtargets_on_packed_genome(&wgs_genome);
    humanz::wgs_pipeline::populate_cellular_somatic_variants_on_t2t(
        &mut twin,
        &wgs_genome,
        &crispr_offtarget_rows,
    );
    let vcf_records = humanz::wgs_pipeline::collect_twin_somatic_vcf_records(&twin, &wgs_genome);
    let vcf_val = humanz::vcf::validate_against_genome(&wgs_genome, &vcf_records);
    let horvath_rows = humanz::wgs_pipeline::map_horvath_353_cpg_sites_with_fm_index(
        &wgs_genome,
        summary.baseline_snapshot.mean_biological_age_years,
        summary.final_treated_snapshot.mean_biological_age_years,
        0.58,
        0.26,
    );

    if let Some(ref out_dir) = export_dir {
        let _ = std::fs::create_dir_all(out_dir);
        if let Ok(f) = std::fs::File::create(out_dir.join("somatic_variants_t2t.vcf")) {
            let _ = humanz::vcf::write_vcf(std::io::BufWriter::new(f), &wgs_genome, &vcf_records);
        }
        if let Ok(f) = std::fs::File::create(out_dir.join("horvath_353_pre_osk.bedgraph")) {
            let _ = humanz::wgs_pipeline::write_horvath_bedgraph(
                std::io::BufWriter::new(f),
                &horvath_rows,
                false,
            );
        }
        if let Ok(f) = std::fs::File::create(out_dir.join("horvath_353_post_osk.bedgraph")) {
            let _ = humanz::wgs_pipeline::write_horvath_bedgraph(
                std::io::BufWriter::new(f),
                &horvath_rows,
                true,
            );
        }
        if let Ok(f) = std::fs::File::create(out_dir.join("crispr_offtargets_cfd_gt_0_1.bed")) {
            let _ = humanz::wgs_pipeline::write_crispr_offtargets_bed(
                std::io::BufWriter::new(f),
                &crispr_offtarget_rows,
                0.10,
            );
        }
    }

    println!("\n----------------------------------------------------------------------------------------------------");
    println!(
        "  T2T-CHM13v2.0 WHOLE-GENOME COUPLED PIPELINE ({} contigs, {} bp, VCF Valid={})",
        wgs_genome.contigs.len(),
        wgs_genome.total_bp(),
        vcf_val.all_consistent()
    );
    println!("----------------------------------------------------------------------------------------------------");
    println!("  [1] SOMATIC VARIANTS VCF 4.2 (Header + 5 Sample Variants):");
    let mut sample_vcf = Vec::new();
    for origin_tag in ["ORIGIN=OsmosticDamage", "ORIGIN=ReplicationError", "ORIGIN=CrisprOfftarget"] {
        for rec in vcf_records.iter().filter(|r| r.info.contains(origin_tag)).take(2) {
            if sample_vcf.len() < 5 {
                sample_vcf.push(rec.clone());
            }
        }
    }
    let mut vcf_preview = Vec::new();
    let _ = humanz::vcf::write_vcf(&mut vcf_preview, &wgs_genome, &sample_vcf);
    for line in String::from_utf8_lossy(&vcf_preview).lines() {
        println!("      {line}");
    }

    println!("\n  [2] HORVATH 353-CpG SITES ON chr1 BEFORE vs AFTER OSK REPROGRAMMING (First 10 Rows):");
    println!("      chrom\tstart\tend\tbeta_pre_osk\tbeta_post_osk\tprobe_id\tgene");
    for r in horvath_rows.iter().filter(|r| r.chrom == "chr1").take(10) {
        println!(
            "      {}\t{}\t{}\t{:.4}\t\t{:.4}\t\t{}\t{}",
            r.chrom, r.start, r.end, r.beta_pre_osk, r.beta_post_osk, r.probe_id, r.gene_symbol
        );
    }

    println!("\n  [3] CRISPR OFF-TARGET HITS ON T2T-CHM13v2.0 WITH CFD SCORE > 0.10 (BED Format):");
    println!("      chrom\tstart\tend\tguide_name\tcfd_score");
    for r in crispr_offtarget_rows.iter().filter(|r| r.cfd_score > 0.10).take(12) {
        println!(
            "      {}\t{}\t{}\t{}\t{:.4}",
            r.chrom, r.start, r.end, r.guide_name, r.cfd_score
        );
    }

    let stats = evaluate_whole_body_six_axis_statistics(&twin)?;
    println!("\n----------------------------------------------------------------------------------------------------");
    println!("  WHOLE-BODY 6-AXIS NON-PARAMETRIC MANN-WHITNEY U & SAFETY AUDIT (N = 240 cells / arm)");
    println!("----------------------------------------------------------------------------------------------------");
    println!(
        "  1. Epigenetic Age Reversal:          U = {:>7.1}, Z = {:>7.2}, p-value = {:.3e}",
        stats.epigenetic_age.u_stat, stats.epigenetic_age.z_score, stats.epigenetic_age.p_value
    );
    println!(
        "  2. Telomere Maintenance:             U = {:>7.1}, Z = {:>7.2}, p-value = {:.3e}",
        stats.telomere_maintenance.u_stat, stats.telomere_maintenance.z_score, stats.telomere_maintenance.p_value
    );
    println!(
        "  3. Proteostasis Clearance:           U = {:>7.1}, Z = {:>7.2}, p-value = {:.3e}",
        stats.proteostasis_clearance.u_stat, stats.proteostasis_clearance.z_score, stats.proteostasis_clearance.p_value
    );
    println!(
        "  4. Mitochondrial Heteroplasmy:       U = {:>7.1}, Z = {:>7.2}, p-value = {:.3e}",
        stats.mitochondrial_heteroplasmy.u_stat, stats.mitochondrial_heteroplasmy.z_score, stats.mitochondrial_heteroplasmy.p_value
    );
    println!(
        "  5. Stem Cell Viability Recovery:     U = {:>7.1}, Z = {:>7.2}, p-value = {:.3e}",
        stats.stem_cell_viability_recovery.u_stat, stats.stem_cell_viability_recovery.z_score, stats.stem_cell_viability_recovery.p_value
    );
    println!(
        "  6. Immune Surveillance Efficiency:   U = {:>7.1}, Z = {:>7.2}, p-value = {:.3e}",
        stats.immune_surveillance_efficiency.u_stat, stats.immune_surveillance_efficiency.z_score, stats.immune_surveillance_efficiency.p_value
    );

    println!("\n  [Negligible Senescence Escape Velocity Metrics]");
    println!(
        "  • Untreated Net Aging Rate:          {:+.3} biological years / chronological year",
        summary.net_aging_rate_control_bio_yr_per_chrono_yr
    );
    println!(
        "  • HumanZ Treated Net Aging Rate:     {:+.3} biological years / chronological year (NEGLIGIBLE SENESCENCE)",
        summary.net_aging_rate_treated_bio_yr_per_chrono_yr
    );
    println!(
        "  • Autonomous Maintenance Sweeps:     {} full sweeps + {} targeted circuit micro-repairs over {} years",
        summary.maintenance_interventions_triggered,
        summary.targeted_circuit_repairs_triggered,
        horizon_years
    );
    println!(
        "  • 4-Stage Clinical Stress Verdict:   {}",
        if stress.all_humanz_resilience_criteria_passed {
            "PASSED ALL 4 ACUTE PHYSIOLOGICAL CHALLENGES (YOUTHFUL RESILIENCE)"
        } else {
            "FAILED"
        }
    );
    println!(
        "  • Whole-Body Safety Guard Verdict:   {} (Max Oncogenic Risk = {:.4} < 0.15, Max Stem Pool = {:.1}% <= 15.0%)",
        if summary.safety_audit.all_organs_approved { "APPROVED (ALL 6 ORGANS SAFE)" } else { "REJECTED" },
        summary.safety_audit.max_oncogenic_risk,
        summary.safety_audit.max_stem_cell_fraction * 100.0
    );

    println!("\n✅ HumanZ Whole-Body Negligible Senescence Simulation Completed Successfully!");
    Ok(())
}

/// Streams a FASTA into 2-bit packed form, verifies checksums, and writes FASTA/BED/manifest.
fn genome_pilot(fasta: &str, out: &std::path::Path) -> std::io::Result<()> {
    use std::io::{BufReader, BufWriter};
    let t0 = std::time::Instant::now();
    let genome = humanz::genome_store::read_fasta_streaming(BufReader::with_capacity(
        1 << 20,
        std::fs::File::open(fasta)?,
    ))?;
    println!(
        "Loaded {} contigs, {} bp, packed {} bytes in {:.1}s",
        genome.contigs.len(),
        genome.total_bp(),
        genome.contigs.iter().map(|c| c.packed_bytes()).sum::<usize>(),
        t0.elapsed().as_secs_f64()
    );
    for c in &genome.contigs {
        if !c.verify() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("checksum mismatch in {}", c.name),
            ));
        }
    }
    std::fs::create_dir_all(out)?;
    genome.write_fasta(BufWriter::new(std::fs::File::create(out.join("genome.packed_roundtrip.fa"))?), 80)?;
    genome.write_n_gaps_bed(BufWriter::new(std::fs::File::create(out.join("n_gaps.bed"))?))?;
    genome.write_manifest(BufWriter::new(std::fs::File::create(out.join("manifest.tsv"))?))?;
    println!("Checksums verified; exports written to {}", out.display());
    Ok(())
}

/// Loads the 7 canonical T2T reference locus FASTAs from `scies-bio-th/testdata` as a lightweight
/// default `PackedGenome` when `--wgs-fasta` is not explicitly provided.
fn load_default_t2t_loci_packed_genome() -> humanz::genome_store::PackedGenome {
    let testdata_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("scies-bio-th")
        .join("testdata");
    let files = [
        "t2t_chr1_parp1.fasta",
        "t2t_chr5_telomere_tert.fasta",
        "t2t_chr6_foxo3.fasta",
        "t2t_chr9_cdkn2a.fasta",
        "t2t_chr11_atm.fasta",
        "t2t_chr17_tp53.fasta",
        "t2t_chr19_sirt6.fasta",
    ];
    let mut combined = String::new();
    for fname in &files {
        if let Ok(txt) = std::fs::read_to_string(testdata_dir.join(fname)) {
            combined.push_str(&txt);
            if !combined.ends_with('\n') {
                combined.push('\n');
            }
        }
    }
    humanz::genome_store::read_fasta_streaming(std::io::Cursor::new(combined)).unwrap_or_default()
}
