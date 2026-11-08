# Empirical state-flow display contract

Question: how do observed disease composition and participant observation
availability evolve across the frozen HRS and ELSA analysis waves?

R is the exclusive plotting backend. Figure 1 compares cohort-specific
observed composition with observed and not-observed counts at each wave.
Figure 2 compares fixed-entry adjacent-wave empirical Sankey flows. Every
initial participant is retained at every wave. These are unweighted exact
sample counts, not model-based transition probabilities; no test or CI is used.

Use final msm participants and analysis entry waves. Real post-entry visits
and disease-item completeness come from augmented exposure panels. At visits
with all three disease items known, preserve the frozen cumulative-history
state; otherwise use Not observed. Never fill living disease states into
unobserved waves. HRS recorded death is placed at a missing-visit wave once
death precedes the wave median actual visit date and remains absorbing. A real
living visit takes precedence at its own wave. ELSA mortality is incompletely
ascertained: no death state is shown and missingness never creates death.

Fixed entry is HRS W8 or ELSA W4. Delayed entrants participate in observed-wave
composition after entry but are outside the fixed-entry flow definition.
Retain all available waves, returns after missingness and direct endpoint
jumps; do not invent disease order between visits. Survey waves are equally
spaced categorically, not by elapsed time. Use shared state colours and 0–100%
axes across comparable panels. Export PDF/SVG/600-dpi TIFF/PNG plus long/wide
states, trajectory counts, adjacent-flow counts, denominators, missing reasons,
captions, Chinese report and numerical/text/collision/alignment QA.
