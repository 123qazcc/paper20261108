# Step 4: empirical state flow and absolute-risk uncertainty

This work completes step 4 of the analysis sequence recommended in the task
on 2026-09-09. It is a post-hoc reporting and uncertainty extension. It does
not modify the frozen primary models, the 11c model specification, or any
eligibility rules. CHARLS remains outside the active scope.

## Absolute quantities

Use all 20 saved HRS five-state and all 20 ELSA W4–W10 four-state primary
M4-MI20 fits. Reconstruct their intensity matrices at each fit's nuisance
model-matrix means, fixing insomnia to exactly 0, 1, 2, or 3 symptoms for the
projection. These are model-based conditional projections under a constant
score and constant-intensity assumption, not causal intervention effects or
individual predictions. Exactly 2 symptoms is not the observed >=2 category.

From S0, calculate at 2, 5, and 10 years:

- Endpoint-state probabilities for every state.
- Cumulative probability of reaching S3 before death by making S3 absorbing
  for this derived calculation; this is distinct from being alive in S3 at
  the horizon in HRS.
- Restricted expected years in each state and the sum of years in S0–S2.
  The HRS sum is alive-and-S3-free time. ELSA omits death and cannot estimate
  alive-and-S3-free time or be used to compare healthy life expectancy.
- Paired differences between fixed score 2 and fixed score 0, with the
  covariance between their predictions retained.

Propagate the full relevant parameter covariance with the delta method using
analytic Frechet derivatives of matrix exponentials within each imputation,
cross-checked against two fixed-step central-difference grids. Combine within- and between-imputation
uncertainty with Rubin's rules. Transform probabilities and restricted years
divided by the horizon to the logit scale for bounded marginal intervals;
use the identity scale for paired differences. Report finite-imputation t
degrees of freedom. Treat empirical nuisance means as fixed conditioning
values; uncertainty in their sampling distribution is not included.

Before inference, verify reconstructed Q and probabilities against msm and
existing per-fit predictions, compare derivative step sizes, check all
probability/occupation sums, and require finite covariance and all 20 fits.
Structural zeros receive deterministic intervals and are identified as such.
Intervals are pointwise 95% model-based confidence intervals, not simultaneous
bands or distribution-free bootstrap intervals.

## Empirical state flow

Show observed state composition by wave, with denominators and missingness.
The empirical display must not use the fitted q matrix as flow counts.
Use a clearly defined entry cohort for participant flows. Preserve a separate
not-observed state for missing visits and do not fill unobserved disease
states for the purpose of drawing a complete path. HRS death may be carried
as an absorbing state after recorded death; ELSA death after 2012 is not
ascertained and cannot be inferred from loss to follow-up. Document delayed
entry, return visits, dropout, observed jumps, and any waves omitted.

Export R-generated PDF/SVG/600-dpi TIFF/PNG with matched source data. Full
tables retain all projected scores, horizons, states, and uncertainty. Main
plots may emphasize the fixed-score 2 versus 0 comparison for readability.
