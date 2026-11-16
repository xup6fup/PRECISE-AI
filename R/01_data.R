precise_baseline_variables <- c("ai_ecg_screening_flag", "holter_within_1y", "male", "age_years", "age_ge80", "clinical_stratum", "af_recency", 
"history_hypertension", "history_diabetes", "history_heart_failure", "history_ckd", "history_copd", "history_pad", 
"history_valvular_disease", "history_cerebrovascular_disease", "history_sleep_apnea", "history_thyroid_disease", 
"history_mi", "history_cabg", "history_ep_study_or_ablation", "mobitz_i_av_block", "borderline_low_heart_rate", 
"rbbb_recency", "lbbb_recency", "left_fascicular_block", "pvc_ecg_recency", "holter_pvc_burden_category", "holter_nsvt", 
"holter_borderline_pause", "ecg_marked_pr_prolongation", "ecg_mild_pr_prolongation", "holter_pr_prolongation", 
"med_beta_blocker", "med_nondhp_ccb", "med_class_ia", "med_class_ic", "med_class_iii", "med_digoxin")

precise_propensity_variables <- c("history_af_flutter", "history_hypertension", "history_diabetes", "history_heart_failure", "history_ckd", "history_copd", 
"history_pad", "history_valvular_disease", "history_cerebrovascular_disease", "history_sleep_apnea", "history_thyroid_disease", 
"history_mi", "history_cabg", "history_ep_study_or_ablation", "med_beta_blocker", "med_nondhp_ccb", "med_class_ia", 
"med_class_ic", "med_class_iii", "med_digoxin")

precise_secondary_outcomes <- c("high_grade_av_block", "complete_av_block", "sick_sinus_syndrome", "ventricular_arrhythmia")

precise_exploratory_outcomes <- c("new_af", "holter_average_hr", "holter_maximum_hr", "holter_minimum_hr", "holter_longest_rr", "holter_pac_burden", 
"holter_pvc_burden", "new_nsvt", "new_bundle_branch_block", "new_pr_prolongation", "new_mobitz_i_av_block")

precise_safety_outcomes <- c("cied_adverse_event", "monitor_adverse_event")

read_precise_data <- function(path = file.path("data", "PRECISE-AI.csv")) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"), colClasses = c(case_id = "character"))
  required <- c("case_id", "group", "selected_matched_cohort", "clinical_holter_stratum", 
"ai_ecg_screening_flag", "holter_within_1y", "male", "age_years", 
"age_ge80", "clinical_stratum", "af_recency", "history_hypertension", 
"history_diabetes", "history_heart_failure", "history_ckd", "history_copd", 
"history_pad", "history_valvular_disease", "history_cerebrovascular_disease", 
"history_sleep_apnea", "history_thyroid_disease", "history_mi", 
"history_cabg", "history_ep_study_or_ablation", "mobitz_i_av_block", 
"borderline_low_heart_rate", "rbbb_recency", "lbbb_recency", 
"left_fascicular_block", "pvc_ecg_recency", "holter_pvc_burden_category", 
"holter_nsvt", "holter_borderline_pause", "ecg_marked_pr_prolongation", 
"ecg_mild_pr_prolongation", "holter_pr_prolongation", "med_beta_blocker", 
"med_nondhp_ccb", "med_class_ia", "med_class_ic", "med_class_iii", 
"med_digoxin", "history_af_flutter", "cied_implantation", "high_grade_av_block", 
"complete_av_block", "sick_sinus_syndrome", "ventricular_arrhythmia", 
"new_af", "holter_average_hr", "holter_maximum_hr", "holter_minimum_hr", 
"holter_longest_rr", "holter_pac_burden", "holter_pvc_burden", 
"new_nsvt", "new_bundle_branch_block", "new_pr_prolongation", 
"new_mobitz_i_av_block", "cied_adverse_event", "monitor_adverse_event", 
"cied_indication", "indication_high_grade_av_block", "indication_complete_av_block",
"indication_sick_sinus_syndrome", "indication_ventricular_arrhythmia",
"other_cied_indication", "cied_type", "cied_implantation_days", 
"cied_indication_days", "ecg_count_1y", "high_risk_ecg_count_1y", 
"last_ecg_days_before_index", "last_high_risk_ecg_days_before_index", 
"ecg_days_before_index", "ecg_count_180d", "recorded_holter_count_180d"
)
  nullable <- c("cied_type", "cied_implantation_days", "cied_indication_days", 
"last_high_risk_ecg_days_before_index", "new_af", "holter_average_hr", 
"holter_maximum_hr", "holter_minimum_hr", "holter_longest_rr", 
"holter_pac_burden", "holter_pvc_burden", "new_nsvt", "new_bundle_branch_block", 
"new_pr_prolongation", "new_mobitz_i_av_block")
  stopifnot(nrow(x) > 0L, !anyDuplicated(names(x)), all(required %in% names(x)), !anyDuplicated(x$case_id), !anyNA(x[setdiff(required, nullable)]), all(x$group %in% c("Intervention", "Control")))
  indication_variables <- c("cied_indication", "indication_high_grade_av_block", "indication_complete_av_block", "indication_sick_sinus_syndrome", "indication_ventricular_arrhythmia", "other_cied_indication")
  stopifnot(all(vapply(x[indication_variables], function(v) all(v %in% 0:1), logical(1))))
  d <- x
  d$row_number <- seq_len(nrow(d))
  d$screening_group <- ifelse(d$ai_ecg_screening_flag == 1, "High risk", "Low risk")
  d$sex_group <- ifelse(d$male == 1, "Male", "Female")
  d$age_group <- ifelse(d$age_ge80 == 1, "80-90 years", "65-79 years")
  d$treatment <- as.integer(d$group == "Intervention")
  stopifnot(all(vapply(d[precise_propensity_variables], function(v) all(v %in% 0:1), logical(1))), all(d$ai_ecg_screening_flag %in% 0:1))
  d
}
