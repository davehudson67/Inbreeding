# Single target model — Joint infection and survival (17 November 2026)

**Decision:** stop adding standalone SEI model variants after Stages 80, 82 and 83. Retain those as exploratory benchmarks. Develop, test and fit **one joint infection-survival model** in a single codebase.

## Scientific estimands

Estimate the separate associations of standardized multilocus homozygosity with (i) susceptible-to-exposed infection acquisition, (ii) exposed-to-infectious progression, and (iii) mortality, while allowing the current infection state to affect mortality. A claim of mediation or causality would require additional assumptions beyond this observational association model.

## Target joint model

- **Latent quarterly states:** S (susceptible), E (exposed), I (infectious), D (dead, absorbing). Transitions S→E, E→I, and S/E/I→D. No forced observed disease state after the last test.
- **Genetic, sex and age covariates:** both infection transitions and mortality; standardize homozygosity once in the original genotyped cohort.
- **Historical infection pressure:** one regularized, smoothly varying calendar baseline from the outset. Avoid a sequence of linear/quadratic/period scientific models. Age effects can be flexible where warranted, but one fixed primary specification must be selected before seeing new genetic results.
- **Mortality:** Gompertz–Makeham age-dependent baseline; add sex, homozygosity and current infection-state effects. Model survival jointly with disease, not by assigning each animal a fixed inferred lifetime infection status.
- **Entry:** explicitly model S/E/I probabilities at first confirmed live capture, including age, sex, calendar year and homozygosity, with left truncation. Acknowledge survivorship/genotyping selection.
- **Observation:** distinguish a known live capture, monitored but uncaptured, and not monitored (survival history codes 1/0/-1). A non-capture is not a known death. Model imperfect capture probability, observing the actual study windows and effort information. Retain repeated captures within quarters and stack usable assay results.
- **Known death recoveries:** integrate recovery information and censoring; resolve dates where death and live capture occur within the same quarter before declaring an exact death interval.
- **Diagnostics:** six distinct assays (BELT, indirect ELISA, culture, DPP, IGRA, StatPak). Fixed original sensitivity, specificity and exposed-state multipliers; each assay's properties constant through time, while availability/test mix can change. Ambiguous outcomes stay missing; not negative. Document provisional mixed-culture rules.
- **Primary exclusions:** the 317 non-confirmed-live original-only badger-quarter occasions remain excluded; never delete the source database records.
- **Cohort:** start with the full 1,603 matched genotyped survival-input individuals, not just the 1,224 with two or more diagnostic quarters. Confirm eligibility/identity and preserve one-test/no-test animals who supply encounter or survival information. Do not impute unobserved genotypes or assume all 1,603 pass checks.

## One implementation, not many models

The next code deliverable should be one **compiled four-state forward likelihood** and its input-preparation pipeline. It must integrate over the unobserved disease/death states and observed capture/diagnostic process. Validate transition ordering against small enumerated synthetic histories, then recover known parameters in simulations generated from **the same full joint model**. Pilot runs and alternate starts test implementation, not new biological model specifications.

Before fitting, check that genetic IDs, q indices, death times, first-capture entry and detection codes agree across the private existing inputs. Do not silently overwrite or repair records. Next fit the one model on the complete qualifying cohort, quantify parameter uncertainty and posterior/predictive fit to diagnostic tests, recaptures and deaths.

## Existing local source inputs

- data/model_inputs/integrated_gm_inputs_group_proxy_2026.rds — 1,603-genotyped survival data, fields observed in existing source scripts: id, y (-1/0/1), T, birth_q, first_q, death_q, hom, sex.
- outputs/sei_transitions/stage76_PRIVATE_CANDIDATE_NOT_VALIDATED.csv — conservative confirmed-live diagnostic candidate matrix (time,id,captured,test1..test6).
- outputs/sei_transitions/stage76_private_id_crosswalk.csv — SEI numeric ID to tattoo mapping.
- data/infection_posteriors/sei_2026/fits/sei_inferred_accuracy_.rds — original six assay characteristics and exposed-state multipliers.
- Stage 80/82/83 output — exploratory reference points, not proof of validity for joint survival.

**Definition of success:** one identified and computationally validated joint model providing genetic effect estimates and credible uncertainty for infection acquisition, progression and mortality, without extrapolating prior disease trajectories as observed truth.