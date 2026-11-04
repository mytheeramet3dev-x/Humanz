//! Whole-Body Statistical Hypothesis Testing (6-Axis Pooled Mann-Whitney U & Multi-Organ Wilcoxon)
//! and Longitudinal Telemetry Reporting for `HumanZ`.

use crate::body::HumanDigitalTwin;
use scies_bio_sim::prelude::*;
use scies_bio_th::stats::test::{mann_whitney_u, wilcoxon_signed_rank, MannWhitneyResult, WilcoxonResult};

/// Pooled whole-body Mann-Whitney U statistical test results across all 6 Hallmarks of Aging.
#[derive(Debug, Clone)]
pub struct WholeBodySixAxisMannWhitney {
    pub epigenetic_age: MannWhitneyResult,
    pub telomere_maintenance: MannWhitneyResult,
    pub proteostasis_clearance: MannWhitneyResult,
    pub mitochondrial_heteroplasmy: MannWhitneyResult,
    pub stem_cell_viability_recovery: MannWhitneyResult,
    pub immune_surveillance_efficiency: MannWhitneyResult,
    pub multi_organ_paired_wilcoxon_age: WilcoxonResult,
}

impl WholeBodySixAxisMannWhitney {
    /// Returns `true` if all 6 Mann-Whitney U axes achieve statistical significance ($p < 0.05$).
    pub fn all_axes_significant(&self) -> bool {
        self.epigenetic_age.p_value < 0.05
            && self.telomere_maintenance.p_value < 0.05
            && self.proteostasis_clearance.p_value < 0.05
            && self.mitochondrial_heteroplasmy.p_value < 0.05
            && self.stem_cell_viability_recovery.p_value < 0.05
            && self.immune_surveillance_efficiency.p_value < 0.05
    }
}

/// Computes pooled whole-body 6-axis Mann-Whitney U tests across all cells in all 6 organs
/// (`control_organs` vs `treated_organs`) plus a paired Wilcoxon test across organ means.
pub fn evaluate_whole_body_six_axis_statistics(
    twin: &HumanDigitalTwin,
) -> Result<WholeBodySixAxisMannWhitney> {
    let mut ctrl_age = Vec::new();
    let mut trt_age = Vec::new();
    let mut ctrl_telo = Vec::new();
    let mut trt_telo = Vec::new();
    let mut ctrl_prot = Vec::new();
    let mut trt_prot = Vec::new();
    let mut ctrl_het = Vec::new();
    let mut trt_het = Vec::new();
    let mut ctrl_viab = Vec::new();
    let mut trt_viab = Vec::new();
    let mut ctrl_imm = Vec::new();
    let mut trt_imm = Vec::new();

    let mut organ_mean_ctrl = Vec::new();
    let mut organ_mean_trt = Vec::new();

    for (&organ, &ctrl_col_id) in &twin.control_organs {
        let trt_col_id = twin.treated_organs[&organ];

        let ctrl_snap = twin.sim.snapshot(ctrl_col_id).unwrap();
        let trt_snap = twin.sim.snapshot(trt_col_id).unwrap();
        organ_mean_ctrl.push(ctrl_snap.stats.mean_epigenetic_age);
        organ_mean_trt.push(trt_snap.stats.mean_epigenetic_age);

        if let Some(col) = twin.sim.colonies.get(ctrl_col_id) {
            for &cid in &col.cell_ids {
                if let Some(c) = twin.sim.cells.get(cid) {
                    ctrl_viab.push(c.viability);
                    if c.is_alive() {
                        ctrl_age.push(c.epigenetic_age);
                        ctrl_telo.push(c.telomere_length_bp);
                        ctrl_prot.push(c.proteostasis.total_burden());
                        ctrl_het.push(c.mitochondria.heteroplasmy_fraction);
                        ctrl_imm.push(c.immune.surveillance_efficiency());
                    }
                }
            }
        }

        if let Some(col) = twin.sim.colonies.get(trt_col_id) {
            for &cid in &col.cell_ids {
                if let Some(c) = twin.sim.cells.get(cid) {
                    trt_viab.push(c.viability);
                    if c.is_alive() {
                        trt_age.push(c.epigenetic_age);
                        trt_telo.push(c.telomere_length_bp);
                        trt_prot.push(c.proteostasis.total_burden());
                        trt_het.push(c.mitochondria.heteroplasmy_fraction);
                        trt_imm.push(c.immune.surveillance_efficiency());
                    }
                }
            }
        }
    }

    Ok(WholeBodySixAxisMannWhitney {
        epigenetic_age: mann_whitney_u(&ctrl_age, &trt_age)?,
        telomere_maintenance: mann_whitney_u(&ctrl_telo, &trt_telo)?,
        proteostasis_clearance: mann_whitney_u(&ctrl_prot, &trt_prot)?,
        mitochondrial_heteroplasmy: mann_whitney_u(&ctrl_het, &trt_het)?,
        stem_cell_viability_recovery: mann_whitney_u(&ctrl_viab, &trt_viab)?,
        immune_surveillance_efficiency: mann_whitney_u(&ctrl_imm, &trt_imm)?,
        multi_organ_paired_wilcoxon_age: wilcoxon_signed_rank(&organ_mean_ctrl, &organ_mean_trt)?,
    })
}
