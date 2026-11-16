# PRECISE-AI analysis project

This R project reproduces the tables and figures for the PRECISE-AI randomized clinical trial. 

The analysis uses the de-identified participant-level dataset. The variable definitions are provided below. Running the project produces **17 individual PDF figures and 12 CSV tables**.

## Analysis input

The project includes the analysis dataset at `data/PRECISE-AI.csv` and a header-only template at `data/PRECISE-AI_template.csv`.

**Participant-level data are not publicly distributed** because the combination of clinical variables carries a potential re-identification risk. The public repository provides the analysis code, variable definitions, and a header-only input template.

The de-identified analysis dataset may be requested from the corresponding author (Chin Lin, xup6fup0629@gmail.com), with access subject to approval following a simplified ethics review. Full individual-patient data that can be linked to longitudinal electronic health records, and raw ECG signals, require formal review and approval by the ethics committee of Tri-Service General Hospital. These full records and raw signals are not required to run this project.

## Run the analyses

Use R with Cairo PDF support. Install the required packages once:

```r
install.packages(c("logistf", "survey", "survival", "jsonlite"), repos = "https://cloud.r-project.org")
```

Open `PRECISE-AI.Rproj` in RStudio and run:

```r
source("run_all.R")
```

The only retained generated files are the **17 PDFs in `output/figures/`** and the **12 CSVs in `output/tables/`**. Repeated runs overwrite these files. Intermediate calculations are held in a temporary directory and removed when the run finishes.

## Project contents

| Path | Contents |
|---|---|
| `PRECISE-AI.Rproj` | RStudio project |
| `run_all.R` | One-command analysis entry point |
| `data/README.md` | Input instructions |
| `data/PRECISE-AI_template.csv` | Header-only input template |
| `R/01_data.R` | Direct CSV input, validation, and analysis groups |
| `R/02_statistics_core.R` | Regression, tables, subgroup analyses, pathways, and sensitivity analyses |
| `R/03_per_protocol.R` | Main outcome figures and propensity-weighted per-protocol analyses |
| `R/04_recent.R` | Recent-ECG post hoc cohort and baseline table |
| `R/05_timing.R` | Event-time analyses and box plots with beeswarm points |
| `R/06_ecg.R` | Distribution of baseline ECG timing |
| `R/07_utilization.R` | ECG and Holter counts during 180-day follow-up |
| `R/08_exports.R` | Export of the 17 PDFs and 12 CSVs |
| `output/figures/` | 17 individual PDF figures |
| `output/tables/` | 12 CSV tables |

## Figure and table outputs

The project produces main Figures 2–3, Supplemental Figures 1–15, Tables 1–3, and Supplemental Tables 1–9. All figures and tables are calculated from the input CSV.

| Filename | Contents |
|---|---|
| `Fig_02.pdf` | Intention-to-treat primary, secondary, and CIED safety outcomes overall and by screening flag |
| `Fig_03.pdf` | Unweighted and propensity-weighted per-protocol outcomes, including both safety outcomes |
| `Supplemental_Fig_01.pdf` | Timing of the most recent ECG, all ECGs, and the most recent high-risk ECG by screening status |
| `Supplemental_Fig_02.pdf` | Mean examinations per participant and proportions receiving at least one ECG, usual-care Holter, or trial-arranged ambulatory ECG monitoring during 180 days, overall and within three screening/protocol strata |
| `Supplemental_Fig_03.pdf` | Primary outcome across clinical strata, Holter-history categories, age, and sex, all randomized participants |
| `Supplemental_Fig_04.pdf` | Primary outcome across clinical strata, Holter-history categories, age, and sex, participants flagged by the screening algorithm |
| `Supplemental_Fig_05.pdf` | Sick sinus syndrome across clinical strata, Holter-history categories, age, and sex, all randomized participants |
| `Supplemental_Fig_06.pdf` | Sick sinus syndrome across clinical strata, Holter-history categories, age, and sex, participants flagged by the screening algorithm |
| `Supplemental_Fig_07.pdf` | Monitoring pathways, CIED indication, and implantation by treatment arm and screening flag |
| `Supplemental_Fig_08.pdf` | CIED indication across clinical strata, Holter-history categories, age, and sex in the per-protocol intervention cohort |
| `Supplemental_Fig_09.pdf` | CIED implantation across clinical strata, Holter-history categories, age, and sex in the per-protocol intervention cohort |
| `Supplemental_Fig_10.pdf` | Post hoc outcomes using the protocol-defined one-year risk classification |
| `Supplemental_Fig_11.pdf` | Post hoc outcomes across three screening/protocol strata |
| `Supplemental_Fig_12.pdf` | Post hoc outcomes in participants with an ECG within 91.5 days before randomization |
| `Supplemental_Fig_13.pdf` | CIED implantation cumulative incidence and treatment effects at 30-day cutoffs |
| `Supplemental_Fig_14.pdf` | Sick sinus syndrome timing using recorded CIED-indication days and treatment effects at 30-day cutoffs |
| `Supplemental_Fig_15.pdf` | CIED-indication timing and indication-to-implantation intervals, with box plots and beeswarm points |
| `Table_1.csv` | Baseline characteristics by randomized group and screening flag |
| `Table_2.csv` | Per-protocol baseline characteristics before and after weighting, with group-specific effective sample sizes |
| `Table_3.csv` | Exploratory Holter findings in the per-protocol cohort before and after weighting |
| `Supplemental_Table_1.csv` | Propensity model summary and weight diagnostics |
| `Supplemental_Table_2.csv` | Distribution of propensity scores and weights |
| `Supplemental_Table_3.csv` | Baseline characteristics by prespecified clinical stratum |
| `Supplemental_Table_4.csv` | Baseline characteristics by clinical stratum without Holter in the preceding year |
| `Supplemental_Table_5.csv` | Baseline characteristics by clinical stratum with Holter in the preceding year |
| `Supplemental_Table_6.csv` | Characteristics of participants meeting CIED indications |
| `Supplemental_Table_7.csv` | Characteristics of participants undergoing CIED implantation |
| `Supplemental_Table_8.csv` | Baseline characteristics across three screening/protocol strata |
| `Supplemental_Table_9.csv` | Baseline characteristics of the recent-ECG analysis cohort |

## Statistical methods

Randomized participants are analyzed in their assigned groups, with control as the regression reference. Table 1 and analogous baseline tables use the full participant count in each column as the denominator for categorical percentages. Baseline and propensity analyses use the complete input covariates. Exploratory Holter comparisons use observations with available values for the corresponding measure.

`ai_ecg_screening_flag` defines **Flagged by the AI-ECG screening algorithm** and **Not flagged by the AI-ECG screening algorithm**. Screening flags use a three-month high-risk ECG window. Protocol-defined risk strata use the preceding year. Post hoc protocol-defined high risk includes flagged participants or participants with `high_risk_ecg_count_1y > 0`; protocol-defined low risk comprises unflagged participants with a zero count. The three strata retain flagged participants separately and divide unflagged participants into protocol-defined high and low risk.

Unweighted treatment effects use binomial logistic regression with Wald 95% confidence intervals and two-sided Wald P values. Firth regression is selected for an analysis block when an evaluable treatment comparison has events in only one arm. All estimable comparisons in that block then use profile penalized-likelihood intervals and penalized likelihood-ratio P values. Overall comparisons and subgroup blocks are assessed separately. Comparisons with an empty arm or no events in either arm remain not estimable and do not themselves trigger Firth regression. Zero-event safety rows remain displayed.

Global treatment-by-subgroup interactions use the method selected for their block. Ordinary models use a likelihood-ratio comparison of interaction and additive models. Firth models jointly constrain the interaction coefficients to zero using `logistf::logistftest`, retaining the full model's Jeffreys penalty. Interactions with empty treatment-by-subgroup cells or a stratum with no events in either arm are not estimable; estimable individual treatment effects are retained.

The per-protocol analyses use the selected matched monitoring cohort. Every participant with `selected_matched_cohort = 1` completed trial-arranged Holter monitoring before inclusion in this cohort. This indicator therefore defines the trial-arranged ambulatory ECG monitoring group (`Protocol_Holter`) in the utilization analyses. An ordinary logistic model estimates intervention-group membership from **20 clinical-history and medication-history variables**: atrial fibrillation/flutter, hypertension, diabetes, heart failure, chronic kidney disease, chronic obstructive pulmonary disease, peripheral arterial disease, valvular heart disease, cerebrovascular disease, obstructive sleep apnea, thyroid disease, prior myocardial infarction, coronary artery bypass grafting, electrophysiology study or catheter ablation, beta-blocker use, non-dihydropyridine calcium-channel blocker use, class IA, IC, and III antiarrhythmic use, and digoxin use.

Inverse-probability weights are `1/e` for intervention and `1/(1-e)` for control. Ordinary weighted outcome models use `survey::svyglm` with robust variance, normal Wald confidence intervals, and residual-degrees-of-freedom t tests; weighted interactions use joint Wald F tests. Weighted Firth models use case weights normalized to sum to the full per-protocol cohort size, with that normalization retained in subsets. Effective sample size is `(sum(w))^2/sum(w^2)` within each arm. Baseline balance uses absolute standardized differences with the unweighted pooled standard deviation as the fixed denominator. Propensity diagnostics include fitted-probability and weight distributions, and effective sample sizes.

Within-cohort indication and implantation proportions use Wilson 95% confidence intervals. Supplemental Figures 8 and 9 report one global Fisher exact test per prespecified subgroup variable, using the Fisher–Freeman–Halton extension for more than two nonempty categories. There are no pairwise tests. Monitoring routes have been confirmed for all participants and are classified using the protocol Holter cohort and the recorded Holter counter.

For participants with sick sinus syndrome (SSS), data review confirmed that the SSS event date is identical to the date of the first CIED indication. When multiple CIED indications are present, SSS is the earliest indication. Accordingly, `cied_indication_days` is used directly as the SSS event time. Participants without SSS are censored at 180 days in the SSS time-to-event analysis.

The recent-ECG sensitivity analysis includes participants whose most recent ECG was recorded 0–91.5 days before randomization.

## Variable definitions

The following definitions correspond to all 77 input columns. Binary clinical flags use 1 for a recorded qualifying condition and 0 when none is recorded. Blank values remain missing or not applicable where specified. Names and categorical values must match the input template. The propensity model uses only the 20 clinical-history and medication-history variables listed above.

| Variable | Type | Coding | Definition |
|---|---|---|---|
| `case_id` | character | Unique string | Randomly assigned participant identifier independent of source order. |
| `group` | character | Intervention; Control | Randomized treatment group. |
| `selected_matched_cohort` | integer | 0 = no; 1 = yes | Membership in the matched cohort; all members completed trial-arranged Holter monitoring. |
| `clinical_holter_stratum` | character | I-a; I-b; II-a; II-b; III-a; III-b (a: without prior-year Holter; b: with prior-year Holter) | Prespecified clinical stratum and prior-Holter subgroup. |
| `ai_ecg_screening_flag` | integer | 1 = flagged; 0 = not flagged | Flagged by the AI-ECG screening algorithm. |
| `holter_within_1y` | integer | 0 = no; 1 = yes | Holter monitoring within the preceding year. |
| `male` | integer | 1 = male; 0 = female | Male sex. |
| `age_years` | integer | Numeric | Age in years. |
| `age_ge80` | integer | 0 = no; 1 = yes | Age at least 80 years. |
| `clinical_stratum` | character | I; II; III | Prespecified clinical stratum. |
| `af_recency` | character | Long-standing; Recent-onset; Without | Atrial fibrillation/flutter recency category. |
| `history_hypertension` | integer | 0 = no; 1 = yes | History of hypertension. |
| `history_diabetes` | integer | 0 = no; 1 = yes | History of diabetes mellitus. |
| `history_heart_failure` | integer | 0 = no; 1 = yes | History of heart failure. |
| `history_ckd` | integer | 0 = no; 1 = yes | History of chronic kidney disease. |
| `history_copd` | integer | 0 = no; 1 = yes | History of chronic obstructive pulmonary disease. |
| `history_pad` | integer | 0 = no; 1 = yes | History of peripheral arterial disease. |
| `history_valvular_disease` | integer | 0 = no; 1 = yes | History of valvular heart disease. |
| `history_cerebrovascular_disease` | integer | 0 = no; 1 = yes | History of cerebrovascular disease. |
| `history_sleep_apnea` | integer | 0 = no; 1 = yes | History of obstructive sleep apnea. |
| `history_thyroid_disease` | integer | 0 = no; 1 = yes | History of thyroid disease. |
| `history_mi` | integer | 0 = no; 1 = yes | History of myocardial infarction. |
| `history_cabg` | integer | 0 = no; 1 = yes | History of coronary artery bypass grafting. |
| `history_ep_study_or_ablation` | integer | 0 = no; 1 = yes | History of electrophysiology study or catheter ablation. |
| `mobitz_i_av_block` | integer | 0 = no; 1 = yes | Mobitz type I atrioventricular block in ECG or Holter. |
| `borderline_low_heart_rate` | integer | 0 = no; 1 = yes | Borderline low heart rate in ECG or Holter. |
| `rbbb_recency` | character | Long-standing; Recent-onset; Without | Complete right bundle branch block recency category. |
| `lbbb_recency` | character | Long-standing; Recent-onset; Without | Complete left bundle branch block recency category. |
| `left_fascicular_block` | integer | 0 = no; 1 = yes | Left fascicular block in ECG. |
| `pvc_ecg_recency` | character | Long-standing; Recent-onset; Without | Premature ventricular complex in ECG recency category. |
| `holter_pvc_burden_category` | character | <5.0% or without Holter; >=10.0%; 5.0-9.9% | Premature ventricular complex burden category in Holter. |
| `holter_nsvt` | integer | 0 = no; 1 = yes | Non-sustained ventricular tachycardia in Holter. |
| `holter_borderline_pause` | integer | 0 = no; 1 = yes | Borderline long pause in Holter. |
| `ecg_marked_pr_prolongation` | integer | 0 = no; 1 = yes | Marked PR prolongation indicator in ECG. |
| `ecg_mild_pr_prolongation` | integer | 0 = no; 1 = yes | Mild PR prolongation indicator in ECG. |
| `holter_pr_prolongation` | integer | 0 = no; 1 = yes | Documented PR prolongation in Holter. |
| `med_beta_blocker` | integer | 0 = no; 1 = yes | Beta-blocker use. |
| `med_nondhp_ccb` | integer | 0 = no; 1 = yes | Non-dihydropyridine calcium-channel blocker use. |
| `med_class_ia` | integer | 0 = no; 1 = yes | Class IA antiarrhythmic use. |
| `med_class_ic` | integer | 0 = no; 1 = yes | Class IC antiarrhythmic use. |
| `med_class_iii` | integer | 0 = no; 1 = yes | Class III antiarrhythmic use. |
| `med_digoxin` | integer | 0 = no; 1 = yes | Digoxin use. |
| `history_af_flutter` | integer | 0 = no; 1 = yes | Binary history of atrial fibrillation/flutter for the propensity model. |
| `cied_implantation` | integer | 0 = no; 1 = yes | Pacemaker or ICD implantation during 180-day follow-up. |
| `high_grade_av_block` | integer | 0 = no; 1 = yes | High-grade atrioventricular block. |
| `complete_av_block` | integer | 0 = no; 1 = yes | Complete atrioventricular block. |
| `sick_sinus_syndrome` | integer | 0 = no; 1 = yes | Sick sinus syndrome. |
| `ventricular_arrhythmia` | integer | 0 = no; 1 = yes | Ventricular arrhythmia. |
| `new_af` | integer | 0 = no; 1 = yes; blank = unavailable or not applicable | New-onset atrial fibrillation. |
| `holter_average_hr` | integer | Numeric; blank = unavailable or not applicable | Average Holter heart rate, beats per minute. |
| `holter_maximum_hr` | integer | Numeric; blank = unavailable or not applicable | Maximum Holter heart rate, beats per minute. |
| `holter_minimum_hr` | integer | Numeric; blank = unavailable or not applicable | Minimum Holter heart rate, beats per minute. |
| `holter_longest_rr` | numeric | Numeric; blank = unavailable or not applicable | Longest Holter RR interval, seconds. |
| `holter_pac_burden` | numeric | Numeric; blank = unavailable or not applicable | Holter premature atrial complex burden, percent. |
| `holter_pvc_burden` | numeric | Numeric; blank = unavailable or not applicable | Holter premature ventricular complex burden, percent. |
| `new_nsvt` | integer | 0 = no; 1 = yes; blank = unavailable or not applicable | New-onset non-sustained ventricular tachycardia. |
| `new_bundle_branch_block` | integer | 0 = no; 1 = yes; blank = unavailable or not applicable | New-onset bundle branch block. |
| `new_pr_prolongation` | integer | 0 = no; 1 = yes; blank = unavailable or not applicable | New-onset PR prolongation. |
| `new_mobitz_i_av_block` | integer | 0 = no; 1 = yes; blank = unavailable or not applicable | New-onset Mobitz type I atrioventricular block. |
| `cied_adverse_event` | integer | 0 = no; 1 = yes | CIED-related adverse event. |
| `monitor_adverse_event` | integer | 0 = no; 1 = yes | Monitor-related adverse event. |
| `cied_indication` | integer | 0 = no; 1 = yes | A qualifying pacemaker or ICD indication during follow-up. |
| `indication_high_grade_av_block` | integer | 0 = no; 1 = yes | High-grade atrioventricular block recorded as a CIED indication. |
| `indication_complete_av_block` | integer | 0 = no; 1 = yes | Complete atrioventricular block recorded as a CIED indication. |
| `indication_sick_sinus_syndrome` | integer | 0 = no; 1 = yes | Sick sinus syndrome recorded as a CIED indication. |
| `indication_ventricular_arrhythmia` | integer | 0 = no; 1 = yes | Ventricular arrhythmia recorded as a CIED indication. |
| `other_cied_indication` | integer | 0 = no; 1 = yes | Other qualifying pacemaker or ICD indication. |
| `cied_type` | character | Pacemaker; ICD; blank if no qualifying implantation | Type of qualifying implanted device. |
| `cied_implantation_days` | integer | Day offset; blank = no observed event | Days from index to CIED implantation. |
| `cied_indication_days` | integer | Day offset; blank = no observed event | Days from index to the first CIED indication; also the confirmed SSS event time for participants with SSS. |
| `ecg_count_1y` | integer | Numeric | Number of ECGs in the preceding year. |
| `high_risk_ecg_count_1y` | integer | Numeric | Number of high-risk ECGs in the preceding year. |
| `last_ecg_days_before_index` | integer | Numeric | Days from the latest ECG to index. |
| `last_high_risk_ecg_days_before_index` | integer | Day offset; blank = no qualifying high-risk ECG | Days from the latest high-risk ECG to index. |
| `ecg_days_before_index` | character | Nonnegative day offsets separated by semicolons | Semicolon-separated ECG day offsets before index. |
| `ecg_count_180d` | integer | Nonnegative count | Number of ECGs during 180-day follow-up. |
| `recorded_holter_count_180d` | integer | Numeric | Recorded Holter count during 180-day follow-up. | 
