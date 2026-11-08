# HRS–ELSA disease multistate analysis contract

## Scientific estimand

- Population: HRS and ELSA participants aged at least 50 years with a valid cardiovascular state, a complete three-symptom insomnia assessment, and at least one later state or death record.
- Exposure: latest valid three-symptom insomnia score at each interval start, analyzed per one-symptom increment and as 0, 1, or at least 2 symptoms.
- Outcome: transition-specific intensity among S0 (no target disease), S1 (hypertension only), S2 (one of heart disease or stroke without hypertension), S3 (at least two target diseases), and S4 (death).
- Disease components: cumulative doctor-diagnosed hypertension, broad heart disease, and stroke. The HRS definitions are unchanged; HRS stroke includes transient ischaemic attack, whereas ELSA explicitly records stroke only.
- ELSA heart phenotype: broad heart disease is the union of angina, myocardial infarction/heart attack, congestive heart failure, heart murmur, abnormal heart rhythm, and other heart trouble. For W4-W9, the Harmonized ELSA source-wide file is the component-level source of record and supplies the W5 `angine` values missing from the `GlobalAgingR` package data; package missingness is not interpreted as absence of disease.
- ELSA W10 confirmation and dispute rule: reconstruction is component specific. `HENvr*=1` confirms a preloaded prior diagnosis and is positive evidence; `HENvr*=2/3/4` disputes that component and retrospectively corrects its earlier reports. Unlabelled raw codes 85/86 remain unknown rather than being recoded as negative.
- ELSA heart ambiguity: when an earlier aggregate `hearte=1` cannot be reconciled after observable component disputes because coexistence of unobserved other heart trouble is not identifiable, the primary phenotype remains tri-state and records the aggregate as unknown. A lower-bound sensitivity treats these ambiguous records as negative; an upper-bound sensitivity retains the aggregate-positive history. These bounds describe the identification range and are not asserted to recover the true diagnosis.
- ELSA myocardial-infarction limitation: W10 has no archived direct `HEEverMI` or `HENvrMI`. Substantive routed responses in `HEAgeMI`, `HEOthHA`, or `HEMedHA` can confirm a prior myocardial infarction, but the available fields cannot support an analogous retrospective MI dispute correction.
- Time origin: first exposure-valid eligible visit; primary time scale is years since that visit.
- Cohort role: HRS is the discovery cohort and ELSA is the external-validation cohort. Models are fit separately; fixed-effect synthesis is secondary and never replaces cohort-specific estimates.
- Association language only. The design does not identify a causal effect of insomnia.

## Allowed transition graph

S0→S1, S0→S2, S0→S4, S1→S3, S1→S4, S2→S3, S2→S4, and S3→S4. Disease histories are cumulative after applying the documented ELSA component-level retrospective corrections, so reverse transitions are not allowed. Unknown component/aggregate values are not forced into a disease-free state. The Q matrix has no separate S0→S3 intensity in either the primary implementation or the supplied sensitivity implementations: an observed S0-to-S3 panel jump is treated as an interval-censored, unobserved two-step path through S1 or S2.

## Model hierarchy

1. Empirical transition counts and person-time rates.
2. Primary analysis: panel-observed continuous-time Markov model using the GlobalAgingR multi-state workflow and `msm::msm` for the protocol-required transition-specific effects. The disease onset times remain interval censored.
3. Approximate-time sensitivity analysis: transition-specific cause-specific Cox models using visit-detection times, with participant-clustered robust standard errors. These estimates are not presented as a replacement for the panel model.
4. Sequential adjustment: M1 age and sex; M2 additionally race/ethnicity, comparable education, wealth quintile and marital status; M3 additionally smoking, drinking, physical inactivity and BMI; M4 (primary) additionally diabetes, lung disease and frequent pain. The available harmonized CES-D cannot be guaranteed to exclude its sleep item, so it is added only in `M4_CESD_sensitivity`.
5. Secondary exposure models: each symptom, cumulative mean, previous valid score, dynamic pattern, and HRS-only extended score/sleep medication.
6. Robustness: complete-case, baseline-weighted, baseline-fixed exposure, directly observed sleep intervals, two-year landmark, alternative Q structure, attained-age/non-homogeneous time, long-interval exclusion, proxy exclusion, semi-Markov approximation and alternative death timing where the necessary source information exists.
7. Sex and age-group heterogeneity, proportional-hazards diagnostics, event-per-parameter checks, and HRS–ELSA heterogeneity statistics.

## Missing-data and weighting contract

- Primary adjusted estimates use 20 multiple imputations of baseline covariates; exposure, state, event time, cohort, and weight are never imputed.
- The primary Markov model is unweighted. Positive baseline respondent/cross-sectional weights carried over follow-up are used for cohort-specific weighted sensitivity analyses only; HRS and ELSA weights are never pooled.
- Complete-case estimates are reported as sensitivity analyses.
- The available harmonized CES-D summary cannot be guaranteed to exclude a sleep-related item in both cohorts. It is therefore omitted from the harmonized primary M4 and added only in an explicitly labelled depression-adjusted sensitivity analysis; this departure from the ideal protocol is disclosed.
- ELSA W4-W10 is used for the non-fatal S0-S3 validation. Because the available Harmonized ELSA EOL A2 death year ends in 2012, ELSA five-state/death results are restricted to that covered period and are labelled limited mortality validation, not a full W10 death replication. Retrospective disease-history correction from W10 does not extend mortality ascertainment beyond 2012.

## Figure contract

- Provisional core conclusion: higher time-varying insomnia burden may mark faster progression through cardiovascular disease states; the figures must show whether this association is transition-specific and whether it replicates in ELSA.
- Results-level question: where in the S0–S4 process does insomnia burden predict progression, and how stable is the signal across cohorts and specifications?
- Archetype: quantitative evidence sequence, with one dominant inferential role per figure.
- Backend: R only (`ggplot2`, `patchwork`, `svglite`, `cairo_pdf`, `ragg`).
- Target/output: generic high-impact-journal double-column figures, 183 mm wide, editable SVG/PDF plus 600-dpi TIFF and review PNG.
- Main-figure sequence:
  - Figure 1: prespecified state structure and all eight allowed intensities.
  - Figure 2: primary non-fatal transition intensity-ratio forest plot.
  - Figure 3: fitted 2/5/10-year state probabilities from S0.
  - Figure 4: fitted state-occupation trajectories.
  - Figure 5: mean sojourn times.
  - Figure 6: direct HRS–ELSA comparison of observed rates and adjusted intensities.
- Extended figures retain death-transition estimates, approximate-time Cox results, complete transition-probability matrices, observed-versus-expected calibration, the full Cox robustness heatmap, and the 16-family protocol sensitivity/coverage audit.
- Statistics: participant/event counts, person-years, hazard ratios, 95% confidence intervals, two-sided p values, interaction p values, and heterogeneity Q/I² where applicable.
- Source data: every quantitative figure has a dedicated CSV generated from the fitted objects.
- Reviewer risks: interval-censored disease onset, residual confounding, non-proportional hazards, sparse transitions, different symptom wording, the HRS-versus-ELSA stroke/TIA difference, W10 MI disputes that are not observable, unidentified coexistence of other heart trouble behind some aggregate-positive ELSA histories, unlabelled W10 codes 85/86 retained as unknown, incomplete ELSA death ascertainment after 2012, and possible CES-D sleep-item overlap.
