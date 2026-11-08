# HRS–ELSA secondary extension: frozen questions and reporting rules

This extension follows inspection of the original multistate results. All new
domain and disease comparisons are post-hoc secondary analyses. The questions
below are fixed before inspecting their newly fitted estimates. They are not
described as preregistered confirmatory analyses.

## Active sequence

1. Formal onset-versus-progression, mortality-origin, and cross-cohort
   coefficient comparisons using the saved primary M4-MI20 models.
2. Joint comparison of difficulty initiating sleep, difficulty maintaining
   sleep, and non-restorative sleep within a single multistate model per
   imputed data set.
3. Disease-specific extension distinguishing broad heart disease from stroke
   after a disease-free state and after hypertension alone.
4. Independent model and figure checks, complete result tables, source data,
   reproducible scripts, and an integrated report distinguishing established
   primary findings from the new secondary findings.

Per the user's latest instruction, no further S0-to-S1 inclusion/exclusion
experiments are active. Previously generated files are retained as records and
are not used to select the analysis population. The S0-to-S1 edge remains part
of the original state graph and of the joint symptom-domain model.

## Authoritative inputs

Use the frozen delivery under
`outputs/hrs_elsa_multistate_final_20260909/data_multistate/`, not the older
root-level data files. Disease-component states are prepared by
`R/11a_prepare_disease_specific_states.R` into
`outputs/hrs_elsa_multistate_nextsteps_20260909/data_disease_specific/`.

HRS is the five-state discovery model; ELSA W4–W10 is the four-state non-fatal
validation model. ELSA mortality ascertainment ends in 2012 and must not be
extended to W10 by assumption. CHARLS is outside the active analysis.

## Joint symptom domains

- Update the three symptom indicators together at a complete observed
  three-item assessment, then carry the triplet forward. Do not combine item
  values updated at different assessments without labelling this as a separate
  sensitivity analysis. Assess changes from the existing asynchronous item
  histories and the primary score.
- Rebuild 20 imputations for missing baseline adjustment covariates with the
  three domains represented in the predictor set. Do not impute exposure,
  outcome state, event time, identifiers, or weights. Avoid including a
  deterministic sum of the three baseline domains alongside those domains.
  Report the remaining limitations of baseline-covariate imputation for a
  time-varying exposure model.
- Fit all three domains jointly. Each domain effect is transition-specific;
  shared nuisance effects follow the primary M4 convention.
- For the four common non-fatal edges, compare equality of all three domain
  effects with an eight-degree-of-freedom omnibus test. Per-edge equality
  tests have two degrees of freedom. Gate specific claims by the global test
  and apply BH to the planned 12 pairwise domain comparisons; retain all
  estimates regardless of the resulting p values.
- Per-domain four-edge association tests address whether a domain is
  associated with any transition. They do not establish differences between
  domains. HRS death-edge domain contrasts are exploratory.
- A death-censored HRS four-state model is a structural-comparability
  sensitivity, not a replacement for the HRS five-state model.

## Disease specificity

- Use eight living combinations of hypertension, broad heart disease, and
  stroke. HRS adds death as a ninth state. The intensity matrix allows only
  single-disease increments, plus HRS death from each living state.
- The four target transitions are no disease to heart disease only, no
  disease to stroke only, hypertension only to hypertension plus heart
  disease, and hypertension only to hypertension plus stroke.
- Preserve observed multi-disease jumps as panel observations. Integrate
  over compatible unobserved paths; do not assign an observed ordering.
- Test heart-versus-stroke coefficient equality separately at the two origin
  states, using covariance from the same fitted model. Pool contrasts within
  and across imputations. Report the two-contrast joint test and multiplicity
  correction for the two planned pairwise comparisons. Cross-cohort
  comparisons retain cohort-specific estimates and phenotype limitations.
- Evaluate sparse non-target edges before estimation. Declare any equality
  constraints and their rationale before inspecting their insomnia effects.
  Do not remove an edge because its association is weak or non-significant.
- Check whether the existing imputation model includes disease-specific
  outcome auxiliaries. Reusing a combined-outcome imputation model requires
  explicit qualification; where feasible, augment the outcome auxiliaries in
  a fresh imputation model and preserve a sensitivity using the original MI.
- HRS heart disease cannot be reliably separated into myocardial infarction,
  angina, or heart failure with these harmonized longitudinal variables. HRS
  stroke includes TIA; ELSA stroke does not explicitly include TIA. Neither
  limitation is resolved by the new state graph.

## Estimation and acceptance

The participant is the independent unit. Retain interval censoring of disease
onset in the primary continuous-time panel Markov likelihood. Time-varying
exposures predict the subsequent interval from its start. Keep approximate
visit-time Cox estimates as secondary sensitivity estimates if used.

Pool estimates and their full within-imputation covariance, then add
between-imputation uncertainty by Rubin's rules. Report the reference
distribution: finite-imputation t for scalar contrasts; any asymptotic joint
Wald chi-square must be explicitly labelled as such. A difference between
significant and non-significant p values is not a test of coefficient
difference. Do not interpret a non-significant contrast as equivalence.

Each accepted fit must have convergence and finite parameter/covariance
diagnostics. A full MI20 result requires 20 accepted fits; partial results
must be labelled partial. Save cache signatures and model inventories so a
restart cannot silently mix specifications. Preserve failed fits and reasons.

Figures use R only, with editable PDF/SVG, raster TIFF/PNG, uncertainty labels,
and matched source-data tables. Retain original automated figure audits and
document any Cairo text-encoding false positives with rendered geometry and
visual evidence. Statistical results are never changed to satisfy graphics
checks.

## Method references

- [msm model, observation-scheme, and coefficient-constraint documentation](https://chjackson.github.io/msm/reference/msm.html)
- [mice pooling rules and variance diagnostics](https://amices.org/mice/reference/pool.html)
- [mice multivariate D1 testing documentation](https://amices.org/mice/reference/D1.html)

The first two references describe the underlying estimation and pooling
operations; the D1 reference does not imply that a large-sample Wald test in
these scripts is the finite-sample D1 statistic.
