# Sleep symptoms and cardiovascular disease transitions in HRS, ELSA and CHARLS

**Public-release candidate, prepared 2026-11-08 (not a deposit date).** This repository accompanies the three-cohort multistate report revised on 2026-10-18. It contains the study's analysis scripts, selected aggregate results and a short numerical check. The date identifies the snapshot; it is not a public release date.

## Start here (no individual data needed)

Unpack the whole archive, open a terminal in the directory containing this README, and run:

~~~sh
Rscript --vanilla check_results.R
~~~

Use R 4.5.1 (tested); no add-on packages, network access or special hardware are needed. The check should take seconds and print **PASS**. A mismatch stops the script with a nonzero exit code. It does not fit models or modify the supplied tables.

The script recalculates HRs, pointwise 95% CIs, two-sided P values and the specified BH corrections from pooled coefficients, SEs and degrees of freedom. It also recalculates P values from supplied joint statistics and pairs selected figure tables with the full-precision outputs by analysis, transition and contrast. This checks numerical reporting; it is not an independent reconstruction of estimates from raw records.

Two small checks also invoke the actual archived functions: the 27 possible hypertension/heart-disease/stroke input combinations (including unknown values) and a synthetic two-coefficient example for Rubin pooling. For the ELSA mortality extension the audit coefficient export leaves its BH column empty; the four nonfatal BH values are checked against the saved Fig. 2 table instead. Death-edge P values are not included in that correction.

## Where each reported result comes from

| Report result | Upstream code | Aggregate reference |
|---|---|---|
| HRS/ELSA: three symptom items, disease states and follow-up | R/01–05c; scripts/publication_audit_20261017.R | reference/followup_summary.csv |
| Main Markov networks and per-symptom HRs (Figs 1–2) | R/07_fit_markov_models.R, 08_collate_multistate_outputs.R, 19a_markov_network_hr.R | figure_data/Figure1*, Figure2* |
| Predicted state probabilities and absolute differences (Figs 3–4) | R/12a_absolute_risk_uncertainty.R, 19b_markov_probabilities.R | figure_data/Figure3*, Figure4* |
| Disease-specific paths and direct heart vs stroke contrasts (Fig 5) | R/11a, 11c_fit; CHARLS R/13a–14e; R/18_integrated_publication_figures.R | figure_data/Figure5a*, Figure5b* |
| CHARLS: added CESD-9 and duration adjustment (Fig 6) | R/16a–16e; R/18_integrated_publication_figures.R | figure_data/Figure6* |
| W2-landmark sleep patterns, direct comparisons and baseline (Figs 7–8) | R/17a–17f, 17i–17o, 18b_integrated_dynamic_figures.R | figure_data/Figure7*, Figure8*; reference/dynamic* |
| Measured-sleep sensitivity and supplementary methods | R/06, 07b; scripts/publication_*_20261017.R | See original report Table 12; this is a Cox sensitivity, not a measured-visits-only Markov refit |

File numbers refer to the original stage labels. CORE is the CHARLS baseline adjustment; BOTH adds the nine non-sleep CESD items (as a score) and sleep duration. HRS/ELSA HRs are per additional symptom; CHARLS uses W1 restless-sleep frequency or a specified W1–W2 pattern contrast. There is no pooled three-cohort HR. ELSA W4–W10 is a **nonfatal four-state** model; the separate five-state mortality extension ends in 2012.

## For investigators with authorized inputs

The scripts expose the study calculations; data acquisition is separate. Use a new private working directory; never point runs at frozen results. Review the staged paths below before invoking scripts. Historical QC checks may require signed preparation bundles, prior-stage reviews and auxiliary plotting tools that are **not** distributed. Do not disable checks or reuse old signatures with different inputs. A clean-room full refit is **not** claimed.

1. Obtain the relevant survey releases independently: [HRS Data Portal](https://hrsdata.isr.umich.edu/), [ELSA data access and UK Data Service](https://www.elsa-project.ac.uk/accessing-elsa-data), [CHARLS data access](https://charls.pku.edu.cn/en/). Follow each provider's rules for the specific products. The HRS analysis covers W8–W15; ELSA uses Harmonized G3, Core W4–W10 and EOL A2; CHARLS uses W1–W5 including raw health, tracking and death-status fields. Obtain GlobalAgingR 2.0.0 and any required harmonized/teaching materials separately with permission. The ELSA extension does not use an updated national mortality linkage.
2. Environment: historical R 4.5.1, msm 1.8.2, mice 3.19.0, GlobalAgingR 2.0.0. See environment.tsv for additional packages discovered in the source and known audit-environment versions. This is not a complete dependency lockfile. Windows/PowerShell was used for historical staging; prefer an ASCII-only private path. The aggregate check is plain base R.
3. HRS/ELSA: inspect scripts/run_cleaning.ps1 (ELSA wave filenames and directory layout), run_exposure_baseline.ps1 and run_multistate_analysis.ps1 (Full profile). Run from the package root; supply absolute private stage/source paths to the wrappers. R/07 uses 20 imputations for the main estimates. R/10, 11a–11c and 12a–12b cover formal contrasts, disease splitting, and derived quantities. Presence of a routine in source is not proof that every optional sensitivity was completed.
4. CHARLS: R/13a–13d prepares and audits the source routing; 14a/14c prepares and reviews the bundle; 14b supports pilot, component N (N=1…20), then pool; 14e independently checks pooling. The 16 and 17 series are separate modules: run preparation and bundle verification before fitting, and preserve the required prior-stage signatures. Both use 20 imputations under each of CORE and BOTH. For example, after the dynamic preconditions are met, the original CLI is `Rscript R/17c_fit_charls_dynamic_models.R component CORE 1`; its `pool` mode requires the full accepted set. 17d independently reconstructs the pooled contrasts. No full-run time estimate is asserted for other hardware.
5. Layout: private/stage/ is the default equivalent of the archived staging area; private/source/ is the licensed upstream data/materials tree; private/model_store/ is for large model objects; private/render_tools/ is for separately obtained production graphics tooling. Keep the original relative hierarchy of outputs/ and analysis/ within your private project. The dated final-output directories are not distributed. Historical figure and publication-audit scripts expect those prior outputs; the aggregate check uses only reference/ and figure_data/. Under private/source/, use reference_scripts/msm_example.R for the licensed reference script, and harmonized/ELSA/ and harmonized/CHARLS/ for the respective harmonized products. The files must be the original authorized versions; a different file with the same name is not a substitute. Documentation, non-scientific configuration labels and path aliases were edited for distribution. Historical and candidate source hashes are listed in provenance.tsv. They do not change the packaged reference estimates. Archived acceptance hashes identify the historical run, not a newly signed run.

## Files and rights

- R/, scripts/: study analysis, independent checks and figure sources; no third-party package source or raw data.
- analysis/: selected stage controls and protocol documents used by the study scripts. Office workflows and machine scheduling are not part of this package.
- reference/: selected full-precision **aggregate** coefficients, contrasts, joint results and follow-up summaries. figure_data/: the ten aggregate source tables for the report’s eight figures. No participant rows, imputed individual records, fitted objects or third-party teaching files are distributed.
- MANIFEST.sha256: SHA-256 of every other delivered file. CITATION.txt and RIGHTS.txt record citation and licensing status. Do not describe this snapshot as publicly deposited until a real versioned URL and date are verified.
