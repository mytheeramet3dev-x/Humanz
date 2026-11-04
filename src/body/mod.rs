//! Multi-Organ Human Digital Twin (`HumanDigitalTwin`) with systemic circulatory SASP,
//! microvascular oxygenation, and immune surveillance coupling across 6 vital organ systems.

use crate::genomics::{ApoeGenotype, ComorbidityProfile, PatientGenotype};
use scies_bio_sim::prelude::*;
use scies_bio_th::anatomy::{
    AnatomicalOrgan, AnatomicalSpatialAtlas, ArterialVascularOrgan, BoneMarrowNicheOrgan,
    CerebralCortexOrgan, DermalCutaneousOrgan, HemodynamicBiomechanicsState, HepaticLobuleOrgan,
    OrganArchitecture, SpatialLayerProfile, SystemicEndocrineAxes, ThymicLymphoidOrgan,
    WholeBodySystemicCircuits,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Vital organ compartments modeled in the `HumanZ` whole-body digital twin.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum OrganSystem {
    /// Vascular endothelium and myocardial microcirculation (governs systemic perfusion & drug delivery).
    CardiovascularEndothelium,
    /// Cerebral cortex neurons, astrocytes, and microglia (highly sensitive to amyloid/tau/lipofuscin & mtDNA).
    CerebralCortexNeuroglia,
    /// Hepatic metabolic parenchyma (central metabolic & xenobiotic clearance hub).
    HepaticMetabolicParenchyma,
    /// Bone marrow hematopoietic & mesenchymal stem cell niche (primary stem cell reservoir).
    HematopoieticBoneMarrowNiche,
    /// Thymic cortex/medulla and peripheral lymphoid organs (governs naive T-cell output & NK/macrophage surveillance).
    ThymicLymphoidOrgan,
    /// Dermal fibroblast and epidermal barrier tissue.
    DermalEpithelialTissue,
}

impl OrganSystem {
    /// Returns all 6 canonical organ systems in deterministic order.
    pub const fn all() -> [OrganSystem; 6] {
        [
            OrganSystem::CardiovascularEndothelium,
            OrganSystem::CerebralCortexNeuroglia,
            OrganSystem::HepaticMetabolicParenchyma,
            OrganSystem::HematopoieticBoneMarrowNiche,
            OrganSystem::ThymicLymphoidOrgan,
            OrganSystem::DermalEpithelialTissue,
        ]
    }

    /// Human-readable lineage name for `OskPulseRegimen`.
    pub const fn lineage_name(&self) -> &'static str {
        match self {
            OrganSystem::CardiovascularEndothelium => "CoronaryMicrovascularEndothelium",
            OrganSystem::CerebralCortexNeuroglia => "CorticalNeuroglialParenchyma",
            OrganSystem::HepaticMetabolicParenchyma => "HepaticSinusoidalHepatocyte",
            OrganSystem::HematopoieticBoneMarrowNiche => "BoneMarrowOsteoblasticNiche",
            OrganSystem::ThymicLymphoidOrgan => "ThymicEpithelialMedulla",
            OrganSystem::DermalEpithelialTissue => "AdultDermalFibroblast",
        }
    }

    /// Instantiates the concrete OOP [`AnatomicalOrgan`] calibrated to `age_years`.
    pub fn build_anatomical_organ(&self, age_years: f64) -> AnatomicalOrgan {
        match self {
            OrganSystem::CardiovascularEndothelium => {
                AnatomicalOrgan::ArterialVascular(ArterialVascularOrgan::from_somatic_age(age_years))
            }
            OrganSystem::CerebralCortexNeuroglia => {
                AnatomicalOrgan::CerebralCortex(CerebralCortexOrgan::from_somatic_age(age_years))
            }
            OrganSystem::HepaticMetabolicParenchyma => {
                AnatomicalOrgan::HepaticLobule(HepaticLobuleOrgan::from_somatic_age(age_years))
            }
            OrganSystem::HematopoieticBoneMarrowNiche => {
                AnatomicalOrgan::BoneMarrowNiche(BoneMarrowNicheOrgan::from_somatic_age(age_years))
            }
            OrganSystem::ThymicLymphoidOrgan => {
                AnatomicalOrgan::ThymicLymphoid(ThymicLymphoidOrgan::from_somatic_age(age_years))
            }
            OrganSystem::DermalEpithelialTissue => {
                AnatomicalOrgan::DermalCutaneous(DermalCutaneousOrgan::from_somatic_age(age_years))
            }
        }
    }
}

/// Snapshot of a single organ's multi-hallmark biological and anatomical metrics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrganBiomarkerSnapshot {
    pub organ: OrganSystem,
    pub colony_id: EntityId,
    pub stats: PopulationStats,
    pub youngs_modulus_kpa: f64,
    pub ages_crosslink_burden: f64,
    pub yap_taz_mechanotransduction_stress: f64,
    pub spatial_layers: Vec<SpatialLayerProfile>,
}

/// Aggregated whole-body systemic biomarkers across all 6 organ compartments.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WholeBodySnapshot {
    pub chronological_age_years: f64,
    pub mean_biological_age_years: f64,
    pub mean_telomere_bp: f64,
    pub mean_heteroplasmy: f64,
    pub mean_proteostasis_burden: f64,
    pub mean_viability: f64,
    pub mean_stem_cell_fraction: f64,
    pub max_stem_cell_fraction: f64,
    pub mean_capillary_density: f64,
    pub mean_po2_mmhg: f64,
    pub mean_immune_surveillance: f64,
    pub systemic_inflammaging_sasp: f64,
    pub max_pluripotency: f64,
    pub hemodynamics: HemodynamicBiomechanicsState,
    pub endocrine: SystemicEndocrineAxes,
    pub circuits: WholeBodySystemicCircuits,
    pub organs: BTreeMap<OrganSystem, OrganBiomarkerSnapshot>,
}

/// Multi-Organ Human Digital Twin orchestrating 6 organ colonies, their [`AnatomicalOrgan`]
/// architectures, closed-loop Moens-Korteweg hemodynamics, multi-organ endocrine axes,
/// and 7 systemic physiological feedback circuits.
pub struct HumanDigitalTwin {
    pub subject_name: String,
    pub chronological_age_years: f64,
    pub cells_per_organ: usize,
    pub genotype: PatientGenotype,
    pub comorbidities: ComorbidityProfile,
    pub sim: SimulationRunner,
    pub control_organs: BTreeMap<OrganSystem, EntityId>,
    pub treated_organs: BTreeMap<OrganSystem, EntityId>,
    pub control_anatomy: BTreeMap<OrganSystem, AnatomicalOrgan>,
    pub treated_anatomy: BTreeMap<OrganSystem, AnatomicalOrgan>,
    pub control_atlas: AnatomicalSpatialAtlas,
    pub treated_atlas: AnatomicalSpatialAtlas,
    pub control_circuits: WholeBodySystemicCircuits,
    pub treated_circuits: WholeBodySystemicCircuits,
    pub control_hemodynamics: HemodynamicBiomechanicsState,
    pub treated_hemodynamics: HemodynamicBiomechanicsState,
    pub control_endocrine: SystemicEndocrineAxes,
    pub treated_endocrine: SystemicEndocrineAxes,
}

impl HumanDigitalTwin {
    /// Initializes a paired Multi-Organ Human Digital Twin (`control_organs` vs `treated_organs`)
    /// at a given `initial_age_years` with population reference genotype and healthy-aged baseline.
    pub fn new(
        subject_name: &str,
        initial_age_years: f64,
        cells_per_organ: usize,
        seed: u64,
    ) -> Self {
        Self::new_personalized(
            subject_name,
            initial_age_years,
            cells_per_organ,
            seed,
            PatientGenotype::wild_type_reference(),
            ComorbidityProfile::healthy_aged(),
        )
    }

    /// Initializes a personalized Multi-Organ Human Digital Twin parameterized by `genotype`
    /// (`APOE`, `CYP2D6/CYP3A4`, `TP53/CDKN2A`, `KL-VS`) and pre-existing `comorbidities`
    /// (`CKD Stage 3b`, `Type 2 Diabetes`, `HFpEF Hypertension`).
    pub fn new_personalized(
        subject_name: &str,
        initial_age_years: f64,
        cells_per_organ: usize,
        seed: u64,
        genotype: PatientGenotype,
        comorbidities: ComorbidityProfile,
    ) -> Self {
        let mut sim = SimulationRunner::new(seed);
        let mut control_organs = BTreeMap::new();
        let mut treated_organs = BTreeMap::new();
        let mut control_anatomy = BTreeMap::new();
        let mut treated_anatomy = BTreeMap::new();
        let human_tax_id = 9606;

        for organ in OrganSystem::all() {
            // Organ-specific baseline senescent burden at advanced age
            let sen_frac = match organ {
                OrganSystem::CardiovascularEndothelium => 0.24,
                OrganSystem::CerebralCortexNeuroglia => 0.18,
                OrganSystem::HepaticMetabolicParenchyma => 0.20,
                OrganSystem::HematopoieticBoneMarrowNiche => 0.22,
                OrganSystem::ThymicLymphoidOrgan => 0.26,
                OrganSystem::DermalEpithelialTissue => 0.22,
            };

            let ctrl_id = sim.create_heterogeneous_somatic_colony(
                human_tax_id,
                cells_per_organ,
                initial_age_years,
                sen_frac,
            );
            let trt_id = sim.create_heterogeneous_somatic_colony(
                human_tax_id,
                cells_per_organ,
                initial_age_years,
                sen_frac,
            );

            // Cerebral cortex has elevated baseline amyloid/tau sensitivity with age (amplified by APOE e4)
            if organ == OrganSystem::CerebralCortexNeuroglia {
                let apoe_amyloid_boost = match genotype.apoe {
                    ApoeGenotype::E4E4 => 0.22,
                    ApoeGenotype::E3E4 => 0.14,
                    ApoeGenotype::E3E3 => 0.08,
                    ApoeGenotype::E2E3 => 0.04,
                };
                let apoe_tau_boost = match genotype.apoe {
                    ApoeGenotype::E4E4 => 0.18,
                    ApoeGenotype::E3E4 => 0.11,
                    ApoeGenotype::E3E3 => 0.07,
                    ApoeGenotype::E2E3 => 0.04,
                };
                for col_id in [ctrl_id, trt_id] {
                    let cids = sim.colonies.get(col_id).unwrap().cell_ids.clone();
                    for cid in cids {
                        if let Some(cell) = sim.cells.get_mut(cid) {
                            cell.proteostasis.amyloid_burden = (cell.proteostasis.amyloid_burden
                                + apoe_amyloid_boost)
                                .clamp(0.0, 0.95);
                            cell.proteostasis.tau_burden = (cell.proteostasis.tau_burden
                                + apoe_tau_boost)
                                .clamp(0.0, 0.95);
                        }
                    }
                }
            }

            let mut ctrl_anat = organ.build_anatomical_organ(initial_age_years);
            let mut trt_anat = organ.build_anatomical_organ(initial_age_years);
            if organ == OrganSystem::CardiovascularEndothelium && comorbidities.hfpef_hypertension {
                ctrl_anat.apply_chronological_aging(16.0);
                trt_anat.apply_chronological_aging(16.0);
            }
            if organ == OrganSystem::HepaticMetabolicParenchyma
                && comorbidities.type2_diabetes_metabolic_syndrome
            {
                ctrl_anat.apply_chronological_aging(15.0);
                trt_anat.apply_chronological_aging(15.0);
            }

            // Anchor each single cell to its histological lamina (Stem Niche vs Transit-Amplifying vs Parenchyma)
            for (col_id, anat) in [(ctrl_id, &ctrl_anat), (trt_id, &trt_anat)] {
                let layers = anat.layers();
                let n_layers = layers.len().max(1);
                let cids = sim.colonies.get(col_id).unwrap().cell_ids.clone();
                for (idx, cid) in cids.into_iter().enumerate() {
                    if let Some(cell) = sim.cells.get_mut(cid) {
                        let lamina_idx = match cell.cell_type {
                            CellType::AdultStemCell => 1.min(n_layers - 1),
                            CellType::TransitAmplifyingProgenitor => 2.min(n_layers - 1),
                            _ => (3 + (idx % (n_layers.saturating_sub(3).max(1)))).min(n_layers - 1),
                        };
                        cell.lamina_index = lamina_idx;
                        let mut depth = 0.0;
                        for l in layers.iter().take(lamina_idx) {
                            depth += l.thickness_um;
                        }
                        depth += 0.5 * layers[lamina_idx].thickness_um;
                        cell.radial_depth_um = depth;
                    }
                }
            }

            control_organs.insert(organ, ctrl_id);
            treated_organs.insert(organ, trt_id);
            control_anatomy.insert(organ, ctrl_anat);
            treated_anatomy.insert(organ, trt_anat);
        }

        let art_e = control_anatomy[&OrganSystem::CardiovascularEndothelium]
            .ecm()
            .youngs_modulus_kpa;
        let control_hemodynamics = HemodynamicBiomechanicsState::from_arterial_modulus(art_e);
        let treated_hemodynamics = control_hemodynamics.clone();
        let mut control_endocrine = SystemicEndocrineAxes::from_somatic_age(initial_age_years);
        if genotype.klotho_vs_heterozygote {
            control_endocrine.circulating_alpha_klotho =
                (control_endocrine.circulating_alpha_klotho + 0.14).clamp(0.12, 1.25);
        }
        if comorbidities.ckd_stage >= 3 {
            control_endocrine.circulating_alpha_klotho =
                (control_endocrine.circulating_alpha_klotho * 0.65).clamp(0.12, 1.25);
        }
        let treated_endocrine = control_endocrine.clone();

        let mut control_atlas = AnatomicalSpatialAtlas::from_somatic_age(initial_age_years);
        if comorbidities.ckd_stage >= 3 {
            let extra = (comorbidities.ckd_stage as f64) * 6.5;
            control_atlas.organs[7].apply_chronological_aging(extra);
            control_atlas.organs[16].apply_chronological_aging(extra * 0.7);
        }
        if comorbidities.type2_diabetes_metabolic_syndrome {
            control_atlas.organs[9].apply_chronological_aging(18.0);
            control_atlas.organs[10].apply_chronological_aging(15.0);
            control_atlas.organs[2].apply_chronological_aging(15.0);
            control_atlas.organs[11].apply_chronological_aging(14.0);
        }
        if comorbidities.hfpef_hypertension {
            control_atlas.organs[0].apply_chronological_aging(16.0);
            control_atlas.organs[6].apply_chronological_aging(16.0);
        }
        if genotype.apoe == ApoeGenotype::E4E4 {
            control_atlas.organs[1].apply_chronological_aging(14.0);
            control_atlas.organs[13].apply_chronological_aging(14.0);
        }
        let treated_atlas = control_atlas.clone();

        let mut control_circuits = WholeBodySystemicCircuits::new();
        control_circuits.step_synchronized(&control_atlas, 1.0, 22.0, 75.0, 0.20);
        let mut treated_circuits = WholeBodySystemicCircuits::new();
        treated_circuits.step_synchronized(&treated_atlas, 1.0, 22.0, 75.0, 0.20);

        Self {
            subject_name: subject_name.to_string(),
            chronological_age_years: initial_age_years,
            cells_per_organ,
            genotype,
            comorbidities,
            sim,
            control_organs,
            treated_organs,
            control_anatomy,
            treated_anatomy,
            control_atlas,
            treated_atlas,
            control_circuits,
            treated_circuits,
            control_hemodynamics,
            treated_hemodynamics,
            control_endocrine,
            treated_endocrine,
        }
    }

    /// Propagates inter-organ circulatory SASP inflammaging, systemic vascular oxygenation,
    /// and circulating thymic/NK immune surveillance across the specified organ map (`control` or `treated`).
    pub fn propagate_inter_organ_systemic_coupling(&mut self, is_treated: bool) {
        let organ_map = if is_treated {
            self.treated_organs.clone()
        } else {
            self.control_organs.clone()
        };

        // 1. Measure circulating endocrine SASP load, cardiovascular capillary density, and thymic immune competence
        let mut total_sasp = 0.0;
        let mut total_live = 0usize;
        let mut cardio_capillary = 0.75;
        let mut thymic_immune = scies_bio_th::aging::ImmunosenescenceState::youthful();

        for (&organ, &col_id) in &organ_map {
            if let Some(col) = self.sim.colonies.get(col_id) {
                if organ == OrganSystem::CardiovascularEndothelium {
                    cardio_capillary = col.vascular.capillary_density;
                }
                if organ == OrganSystem::ThymicLymphoidOrgan {
                    thymic_immune = col.immune.clone();
                }
                for &cid in &col.cell_ids {
                    if let Some(cell) = self.sim.cells.get(cid) {
                        if cell.is_alive() {
                            total_sasp += cell.sasp_secretion;
                            total_live += 1;
                        }
                    }
                }
            }
        }

        let systemic_sasp = if total_live > 0 {
            total_sasp / total_live as f64
        } else {
            0.0
        };

        // Update Moens-Korteweg hemodynamic biomechanics from arterial Young's modulus
        let art_e = if is_treated {
            self.treated_anatomy[&OrganSystem::CardiovascularEndothelium]
                .ecm()
                .youngs_modulus_kpa
        } else {
            self.control_anatomy[&OrganSystem::CardiovascularEndothelium]
                .ecm()
                .youngs_modulus_kpa
        };
        let hemo = HemodynamicBiomechanicsState::from_arterial_modulus(art_e);
        let barotrauma = hemo.microvascular_barotrauma_index;
        if is_treated {
            self.treated_hemodynamics = hemo;
        } else {
            self.control_hemodynamics = hemo;
        }

        let klotho = if is_treated {
            self.treated_endocrine.circulating_alpha_klotho
        } else {
            self.control_endocrine.circulating_alpha_klotho
        };

        // Advance closed-loop whole-body physiological circuits across all 20 organs
        let (glymphatic_clearance, portal_lps) = if is_treated {
            self.treated_circuits.step_synchronized(
                &self.treated_atlas,
                1.0,
                22.0,
                75.0,
                systemic_sasp,
            );
            (
                self.treated_circuits
                    .neuroendocrine_glymphatic
                    .cortical_glymphatic_clearance_rate,
                self.treated_circuits
                    .enterohepatic
                    .portal_lps_endotoxin_eu_ml,
            )
        } else {
            self.control_circuits.step_synchronized(
                &self.control_atlas,
                1.0,
                22.0,
                75.0,
                systemic_sasp,
            );
            (
                self.control_circuits
                    .neuroendocrine_glymphatic
                    .cortical_glymphatic_clearance_rate,
                self.control_circuits
                    .enterohepatic
                    .portal_lps_endotoxin_eu_ml,
            )
        };

        // 2. Apply systemic coupling across all organs (Vascular + Thymic + Barotrauma + Klotho + Circuit Feedback + Lamina Anchoring)
        let anatomy_snapshot = if is_treated {
            self.treated_anatomy.clone()
        } else {
            self.control_anatomy.clone()
        };

        for (&organ, &col_id) in &organ_map {
            let anat = &anatomy_snapshot[&organ];
            let layers = anat.layers();
            let n_layers = layers.len().max(1);

            let cids = if let Some(col) = self.sim.colonies.get_mut(col_id) {
                // Cardiovascular perfusion couples to peripheral organ capillary beds (20% systemic blend)
                col.vascular.capillary_density = (0.80 * col.vascular.capillary_density
                    + 0.20 * cardio_capillary
                    - 0.015 * barotrauma
                    + 0.012 * klotho)
                    .clamp(0.20, 0.98);
                // Thymic lymphoid output couples to peripheral immune surveillance (35% systemic blend)
                col.immune.thymic_output_naive_t = (0.65 * col.immune.thymic_output_naive_t
                    + 0.35 * thymic_immune.thymic_output_naive_t)
                    .clamp(0.05, 1.0);
                col.immune.nk_cytotoxicity = (0.65 * col.immune.nk_cytotoxicity
                    + 0.35 * thymic_immune.nk_cytotoxicity)
                    .clamp(0.10, 1.0);
                col.cell_ids.clone()
            } else {
                continue;
            };

            for (idx, cid) in cids.into_iter().enumerate() {
                if let Some(cell) = self.sim.cells.get_mut(cid) {
                    if !cell.is_alive() {
                        continue;
                    }
                    // Anchor cell to its histological lamina (Stem Niche = 1, Progenitor = 2, Parenchyma >= 3)
                    let lamina_idx = match cell.cell_type {
                        CellType::AdultStemCell => 1.min(n_layers - 1),
                        CellType::TransitAmplifyingProgenitor => 2.min(n_layers - 1),
                        _ => (3 + (idx % (n_layers.saturating_sub(3).max(1)))).min(n_layers - 1),
                    };
                    cell.lamina_index = lamina_idx;
                    let mut depth = 0.0;
                    for l in layers.iter().take(lamina_idx) {
                        depth += l.thickness_um;
                    }
                    depth += 0.5 * layers[lamina_idx].thickness_um;
                    cell.radial_depth_um = depth;

                    // Bi-directional feedback from Neuro-Glymphatic circuit to Cerebral Cortex proteostasis
                    if organ == OrganSystem::CerebralCortexNeuroglia {
                        let clearance_boost = ((glymphatic_clearance - 12.0).max(0.0) / 25.0)
                            .clamp(0.0, 0.015);
                        cell.proteostasis.amyloid_burden =
                            (cell.proteostasis.amyloid_burden - clearance_boost).clamp(0.005, 0.95);
                        cell.proteostasis.tau_burden =
                            (cell.proteostasis.tau_burden - clearance_boost * 0.85)
                                .clamp(0.005, 0.95);
                    }

                    // Bi-directional feedback from Enterohepatic portal LPS endotoxin to Hepatic Zone inflammation
                    if organ == OrganSystem::HepaticMetabolicParenchyma && portal_lps > 0.35 {
                        let lps_stress = (portal_lps - 0.35) * 0.015;
                        cell.sasp_secretion =
                            (cell.sasp_secretion + lps_stress).clamp(0.01, 0.95);
                    }

                    // If systemic SASP inflammaging is high (> 0.18), endocrine cytokines accelerate multi-organ aging
                    if systemic_sasp > 0.18 {
                        cell.adjust_epigenetic_age(systemic_sasp * 0.35);
                        cell.dna_damage =
                            (cell.dna_damage + systemic_sasp * 0.012).clamp(0.0, 1.0);
                    }
                }
            }
        }
    }

    /// Simulates `delta_years` of chronological time and natural biological entropy across an organ map.
    pub fn apply_chronological_entropy(&mut self, delta_years: f64, is_treated: bool) {
        let organ_map = if is_treated {
            self.treated_organs.clone()
        } else {
            self.control_organs.clone()
        };

        // Advance chronological ECM stiffening, microvascular rarefaction & endocrine decay
        let endo = if is_treated {
            for anat in self.treated_anatomy.values_mut() {
                anat.apply_chronological_aging(delta_years);
            }
            self.treated_atlas.apply_chronological_aging(delta_years);
            &mut self.treated_endocrine
        } else {
            for anat in self.control_anatomy.values_mut() {
                anat.apply_chronological_aging(delta_years);
            }
            self.control_atlas.apply_chronological_aging(delta_years);
            &mut self.control_endocrine
        };
        endo.circulating_alpha_klotho =
            (endo.circulating_alpha_klotho - 0.006 * delta_years).clamp(0.12, 1.25);
        endo.hepatic_igf1_level =
            (endo.hepatic_igf1_level - 0.005 * delta_years).clamp(0.15, 1.20);
        endo.adrenal_cortisol_to_dhea_ratio =
            (endo.adrenal_cortisol_to_dhea_ratio + 0.012 * delta_years).clamp(0.80, 4.50);
        endo.hypothalamic_gnrh_scn_tone =
            (endo.hypothalamic_gnrh_scn_tone - 0.0055 * delta_years).clamp(0.15, 1.0);

        for (&organ, &col_id) in &organ_map {
            let yap_taz = if is_treated {
                self.treated_anatomy[&organ]
                    .ecm()
                    .yap_taz_mechanotransduction_stress()
            } else {
                self.control_anatomy[&organ]
                    .ecm()
                    .yap_taz_mechanotransduction_stress()
            };
            let cids = if let Some(col) = self.sim.colonies.get_mut(col_id) {
                col.vascular.capillary_density =
                    (col.vascular.capillary_density - 0.0055 * delta_years).clamp(0.22, 0.96);
                col.vascular.enos_activity =
                    (col.vascular.enos_activity - 0.0060 * delta_years).clamp(0.20, 0.96);
                col.vascular.po2_mmhg = (14.0
                    + 31.0 * (0.65 * col.vascular.capillary_density + 0.35 * col.vascular.enos_activity))
                    .clamp(15.0, 46.0);
                col.niche.wnt_signal =
                    (col.niche.wnt_signal - 0.005 * delta_years).clamp(0.15, 0.96);
                col.niche.notch_signal =
                    (col.niche.notch_signal - 0.006 * delta_years).clamp(0.15, 0.96);
                col.immune.thymic_output_naive_t =
                    (col.immune.thymic_output_naive_t - 0.008 * delta_years).clamp(0.04, 0.98);
                col.immune.nk_cytotoxicity =
                    (col.immune.nk_cytotoxicity - 0.005 * delta_years).clamp(0.15, 0.98);
                col.cell_ids.clone()
            } else {
                continue;
            };

            let mut sasp_sum = 0.0;
            let mut live_n = 0usize;
            for &cid in &cids {
                if let Some(c) = self.sim.cells.get(cid) {
                    if c.is_alive() {
                        sasp_sum += c.sasp_secretion;
                        live_n += 1;
                    }
                }
            }
            let colony_mean_sasp = if live_n > 0 { sasp_sum / live_n as f64 } else { 0.0 };

            for cid in cids {
                if let Some(cell) = self.sim.cells.get_mut(cid) {
                    if !cell.is_alive() {
                        continue;
                    }
                    let plasticity = ((cell.pluripotency - 0.01).max(0.0) / 0.15).clamp(0.0, 2.0);
                    if is_treated && plasticity > 0.05 {
                        // Paracrine bystander senescence & unbuffered YAP/TAZ mechanotransduction
                        // when recently reprogrammed cells sit in an unconsolidated stiff/oxidative niche
                        let bystander_drive = (colony_mean_sasp - 0.06).max(0.0) * 0.28
                            + (yap_taz - 0.12).max(0.0) * 0.22
                            + (cell.mitochondria.ros_level - 0.18).max(0.0) * 0.18;
                        if bystander_drive > 0.0 {
                            cell.sasp_secretion = (cell.sasp_secretion
                                + bystander_drive * delta_years)
                                .clamp(0.01, 0.95);
                        }
                    }

                    // Chromatin plasticity + microenvironmental discordance rebound:
                    // Newly reprogrammed cells (elevated pluripotency) placed in an unconsolidated,
                    // stiff (YAP/TAZ), oxidative (ROS), or inflammatory (SASP) niche experience
                    // accelerated epigenetic re-methylation drift upon treatment withdrawal.
                    let microenv_stress = 0.38 * yap_taz
                        + 0.32 * cell.sasp_secretion
                        + 0.18 * (cell.mitochondria.ros_level - 0.10).max(0.0)
                        + 0.12 * (cell.mitochondria.heteroplasmy_fraction - 0.08).max(0.0);
                    let rebound_drift = if is_treated && plasticity > 0.05 && microenv_stress > 0.08 {
                        plasticity * (microenv_stress - 0.08) * 2.2 * delta_years
                    } else {
                        0.0
                    };

                    cell.adjust_epigenetic_age(delta_years + rebound_drift);
                    cell.adjust_telomere(-(52.0 + 35.0 * rebound_drift) * delta_years);
                    cell.dna_damage =
                        (cell.dna_damage + (0.006 + 0.015 * rebound_drift) * delta_years)
                            .clamp(0.005, 0.98);
                    cell.mitochondria.heteroplasmy_fraction =
                        (cell.mitochondria.heteroplasmy_fraction + 0.0022 * delta_years)
                            .clamp(0.005, 0.90);
                    cell.mitochondria.ros_level =
                        (cell.mitochondria.ros_level + 0.004 * delta_years).clamp(0.05, 0.95);
                    cell.proteostasis.amyloid_burden =
                        (cell.proteostasis.amyloid_burden + 0.003 * delta_years).clamp(0.0, 0.95);
                    cell.proteostasis.tau_burden =
                        (cell.proteostasis.tau_burden + 0.0025 * delta_years).clamp(0.0, 0.95);
                    cell.proteostasis.lipofuscin_burden =
                        (cell.proteostasis.lipofuscin_burden + 0.0035 * delta_years).clamp(0.0, 0.95);
                    cell.immune.thymic_output_naive_t =
                        (cell.immune.thymic_output_naive_t - 0.008 * delta_years).clamp(0.04, 0.98);
                    cell.immune.nk_cytotoxicity =
                        (cell.immune.nk_cytotoxicity - 0.005 * delta_years).clamp(0.15, 0.98);
                    cell.viability = (cell.viability - 0.0035 * delta_years).clamp(0.55, 0.99);
                    if cell.telomere_length_bp <= 3_800.0 || cell.dna_damage > 0.70 {
                        cell.sasp_secretion = (cell.sasp_secretion + 0.15).clamp(0.65, 0.95);
                    }
                }
            }
        }
    }

    /// Captures a complete `WholeBodySnapshot` across all 6 organs for either the `treated` or `control` twin.
    pub fn snapshot_whole_body(&self, is_treated: bool) -> WholeBodySnapshot {
        let (organ_map, anatomy_map) = if is_treated {
            (&self.treated_organs, &self.treated_anatomy)
        } else {
            (&self.control_organs, &self.control_anatomy)
        };

        let mut organs = BTreeMap::new();
        let mut sum_age = 0.0;
        let mut sum_telo = 0.0;
        let mut sum_het = 0.0;
        let mut sum_prot = 0.0;
        let mut sum_viab = 0.0;
        let mut sum_stem = 0.0;
        let mut max_stem: f64 = 0.0;
        let mut sum_cap = 0.0;
        let mut sum_po2 = 0.0;
        let mut sum_imm = 0.0;
        let mut max_pluri: f64 = 0.0;
        let mut sum_sasp = 0.0;
        let mut live_cells = 0usize;

        for (&organ, &col_id) in organ_map {
            let snap = self.sim.snapshot(col_id).unwrap();
            let anat = &anatomy_map[&organ];
            let shuttle_eff = if is_treated { 0.92 } else { 0.0 };
            let spatial_layers = anat.solve_spatial_pde_profiles(
                scies_bio_th::aging::MW_MRNA_LNP_VECTOR_DA,
                shuttle_eff,
            );

            sum_age += snap.stats.mean_epigenetic_age;
            sum_telo += snap.stats.mean_telomere_bp;
            sum_het += snap.stats.mean_heteroplasmy;
            sum_prot += snap.stats.mean_proteostasis_burden;
            sum_viab += snap.stats.mean_viability;
            sum_stem += snap.stats.stem_cell_fraction;
            max_stem = max_stem.max(snap.stats.stem_cell_fraction);
            sum_cap += snap.stats.capillary_density;
            sum_po2 += snap.stats.po2_mmhg;
            sum_imm += snap.stats.mean_immune_surveillance;
            max_pluri = max_pluri.max(snap.stats.mean_pluripotency);

            if let Some(col) = self.sim.colonies.get(col_id) {
                for &cid in &col.cell_ids {
                    if let Some(c) = self.sim.cells.get(cid) {
                        if c.is_alive() {
                            sum_sasp += c.sasp_secretion;
                            live_cells += 1;
                        }
                    }
                }
            }

            organs.insert(
                organ,
                OrganBiomarkerSnapshot {
                    organ,
                    colony_id: col_id,
                    stats: snap.stats,
                    youngs_modulus_kpa: anat.ecm().youngs_modulus_kpa,
                    ages_crosslink_burden: anat.ecm().ages_crosslink_burden,
                    yap_taz_mechanotransduction_stress: anat
                        .ecm()
                        .yap_taz_mechanotransduction_stress(),
                    spatial_layers,
                },
            );
        }

        let n = organs.len().max(1) as f64;
        WholeBodySnapshot {
            chronological_age_years: self.chronological_age_years,
            mean_biological_age_years: sum_age / n,
            mean_telomere_bp: sum_telo / n,
            mean_heteroplasmy: sum_het / n,
            mean_proteostasis_burden: sum_prot / n,
            mean_viability: sum_viab / n,
            mean_stem_cell_fraction: sum_stem / n,
            max_stem_cell_fraction: max_stem,
            mean_capillary_density: sum_cap / n,
            mean_po2_mmhg: sum_po2 / n,
            mean_immune_surveillance: sum_imm / n,
            systemic_inflammaging_sasp: if live_cells > 0 {
                sum_sasp / live_cells as f64
            } else {
                0.0
            },
            max_pluripotency: max_pluri,
            hemodynamics: if is_treated {
                self.treated_hemodynamics.clone()
            } else {
                self.control_hemodynamics.clone()
            },
            endocrine: if is_treated {
                self.treated_endocrine.clone()
            } else {
                self.control_endocrine.clone()
            },
            circuits: if is_treated {
                self.treated_circuits.clone()
            } else {
                self.control_circuits.clone()
            },
            organs,
        }
    }
}
