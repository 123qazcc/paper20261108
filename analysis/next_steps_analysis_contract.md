# HRS–ELSA next-step inference and figure contract

## Analysis scope

- These are post-hoc secondary analyses developed after inspection of the
  primary results. The contrasts and multiplicity families below were fixed
  before their evaluation in this extension; they are not presented as
  preregistered confirmatory tests from the original protocol.
- Reuse the frozen HRS five-state and ELSA Wave 4–10 four-state M4-MI20
  continuous-time Markov fits. No disease states, exposure values, covariates,
  event times, or imputations are changed in this first extension.
- The independent unit remains the participant. Repeated panel observations are
  handled inside the fitted continuous-time Markov likelihood; HRS and ELSA are
  independent cohorts.
- Every linear contrast is evaluated inside each of the 20 imputed-data fits and
  then pooled with Rubin's rules. Scalar contrasts use Rubin total variance and
  finite-imputation degrees of freedom. Joint tests use the Rubin total
  covariance with a large-sample Wald chi-square reference, matching the
  existing primary joint-test convention.
- All tests are two-sided. A ratio of transition intensity ratios above 1 means
  that the insomnia association is larger for the transition in the numerator.
- The two path-specific onset-versus-progression contrasts form one
  planned family for this extension and receive BH-FDR adjustment. HRS death pairwise
  contrasts are exploratory and receive BH-FDR adjustment across all six
  state-pair comparisons.
- Cross-cohort comparisons use independent-cohort covariance. Per-transition
  Cochran Q, I-squared, a fixed-effect synthesis, and an exploratory
  DerSimonian-Laird random-effects synthesis are reported. With only two
  cohorts, tau-squared and I-squared are explicitly labelled unstable; a
  modified Hartung-Knapp interval is supplied for the random-effects summary.

## Figure contracts

### Figure 8 — disease stage discrimination

- Core conclusion: formal coefficient contrasts determine whether insomnia is
  more strongly associated with later cardiovascular progression than with
  first disease onset.
- Results-level question: do the two progression-edge coefficients exceed their
  corresponding onset-edge coefficients in HRS and ELSA?
- Archetype: quantitative comparison.
- Evidence: two path-specific ratios plus the mean progression-versus-onset
  ratio, with 95% confidence intervals and cohort-specific joint-test p values.

### Extended Data Figure 7 — HRS mortality-origin heterogeneity

- Core conclusion: formal contrasts determine whether the insomnia–mortality
  association differs by the cardiovascular state from which death occurs.
- Results-level question: are the four HRS death-transition coefficients equal?
- Archetype: quantitative comparison.
- Evidence: all six state-pair ratios, 95% confidence intervals, the omnibus
  equality test, and BH-FDR status.

### Figure 9 — cross-cohort replication and synthesis

- Core conclusion: HRS and ELSA are compared with formal heterogeneity tests,
  while pooled estimates remain secondary to cohort-specific estimates.
- Results-level question: which non-fatal transition associations are compatible
  across cohorts, and what is their cautiously synthesized magnitude?
- Archetype: faceted quantitative forest plot.
- Evidence: HRS, ELSA, fixed-effect, and random-effects estimates for each of the
  four harmonized non-fatal transitions, with Q and I-squared annotations.

## Export and integrity contract

- R is the exclusive plotting backend (`ggplot2`).
- Each figure is exported as editable PDF/SVG plus 600-dpi TIFF and PNG.
- Every figure has a dedicated source-data CSV and no synthetic observations.
- Minimum rendered text is 5 pt; uncertainty definitions and multiplicity are
  stated in source data and the generated report.
- Final PDFs undergo the Nature-figure text and collision audits, followed by
  visual inspection at the exported physical size.
