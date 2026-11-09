# Bunny 10,000-Year Computational Experiment

> Status: experiment specification only. No 10,000-year simulation result is claimed by this document.

## Goal

Compare whether the current HumanZ intervention strategy and a one-time-intervention strategy maintain modeled health over a 10,000-year horizon. Bunny is the name of the simulated subject.

This is a computational stress test of model assumptions, not a prediction that a human could live for 10,000 years.

## Experimental arms

Run each arm from the same initial age, genotype, comorbidity profile, cell count, and random seed wherever possible.

1. **Control** — no HumanZ rejuvenation intervention. Keep the model's ordinary baseline biology.
2. **Bunny A: one-shot** — apply the initial induction/rejuvenation phase once. Disable subsequent *additional anti-aging intervention sweeps* and targeted HumanZ rescue actions. Do not disable ordinary physiological functions (for example, baseline DNA repair, immune function, cell turnover, and wound repair) if those are represented by the underlying simulator. First audit the simulator to distinguish those baseline processes from the added HumanZ interventions; do not assume they are currently modeled.
3. **Bunny B: closed-loop maintenance** — use the current HumanZ strategy: initial induction plus threshold-triggered rejuvenation and targeted circuit repairs.

The existing `run_negligible_senescence_lifespan_simulation` currently performs four induction rounds, then checks yearly thresholds and may perform further rejuvenation sweeps. It is therefore not a one-shot model.

## Initial conditions

- Starting chronological age: 72 years (current CLI default).
- Horizon: 10,000 chronological years.
- Current paired digital-twin design: six modeled organ compartments and 40 cells per organ per twin (240 cells per arm).
- Use identical seeds and inputs across arms for paired comparisons. Repeat with multiple seeds once the deterministic single-seed path is validated.

## Required yearly outputs

Record a snapshot at years 0, 1, 5, 10, 25, 50, 100, 250, 500, 1,000, 2,500, 5,000, 7,500, and 10,000:

- chronological age and mean biological age
- telomere length
- mtDNA heteroplasmy
- proteostasis burden
- cell viability and live/dead cell counts
- stem-cell fraction and maximum stem-cell fraction
- capillary density and oxygenation
- immune surveillance and systemic SASP/inflammation
- organ-specific biomarker values
- number and type of extra interventions
- safety flags, terminal events, and reason for termination

Do not rely only on averages: retain per-organ and distributional summaries so a single failing organ or high-risk subpopulation is not hidden by whole-body means.

## Endpoints and failure handling

A simulation must not label Bunny as alive solely because a biomarker is clamped to a nonzero floor. Define explicit terminal criteria and report the first year and reason a criterion is met. Candidate model-level endpoints include irreversible loss of a required organ function, depletion of viable cells below a declared threshold, uncontrolled modeled oncogenic expansion, or a fatal stochastic event if such events are implemented.

Keep accident, infection, and other external mortality separate from senescence. If they are not modeled, state that they are excluded rather than implying immunity to them.

## Numerical and model checks

1. Run short horizons (30, 100, 500, and 1,000 years) before attempting 10,000.
2. Check that every state variable remains finite and within biologically meaningful bounds; document every clamp and its consequences.
3. Verify cell birth, death, replacement, and population-size behavior. A fixed 40-cell sample must not silently stand in for an indefinitely maintained human cell population.
4. Check time-step stability and computational cost at 10,000 annual steps.
5. Add tests for the one-shot arm proving that extra rejuvenation/repair interventions are not invoked after induction, while baseline simulator functions remain enabled.
6. Add tests that the maintenance arm still invokes interventions when thresholds are crossed.
7. Export machine-readable longitudinal data (CSV/JSON), not only a few hard-coded display years.

## Interpretation

Report outcomes as conditional on the model and parameter choices. A passing test means the code behaves as specified; it does not establish biological efficacy, safety, or human longevity. Do not convert a simulated 10,000-year trajectory into a clinical claim.