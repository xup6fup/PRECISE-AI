precise_core_helpers <- function(envir = parent.frame()) {
    evalq({
        stopifnot(requireNamespace("survey", quietly = TRUE), requireNamespace("logistf", quietly = TRUE))
        logistf <- logistf::logistf
        logistftest <- logistf::logistftest
        ctrl <- logistf::logistf.control(maxit = 1000, collapse = FALSE)
        pctrl <- logistf::logistpl.control(maxit = 1000)
        safe <- function(expr) {
            w <- character()
            v <- tryCatch(withCallingHandlers(expr, warning = function(e) {
                w <<- c(w, conditionMessage(e))
                invokeRestart("muffleWarning")
            }), error = function(e) {
                w <<- c(w, conditionMessage(e))
            })
            list(value = v, warnings = paste(unique(w), collapse = "; "))
        }
        outs <- c("cied_implantation", precise_secondary_outcomes)
        outlabs <- setNames(c("CIED implantation", "High-grade atrioventricular block", "Complete atrioventricular block", "Sick sinus syndrome", 
            "Ventricular arrhythmia"), outs)
        subgroup_levels <- list(group = c("Control", "Intervention"), screening_group = c("High risk", "Low risk"), sex_group = c("Female", 
            "Male"), clinical_stratum = c("I", "II", "III"), clinical_holter_stratum = c("I-a", "I-b", "II-a", "II-b", "III-a", "III-b"), 
            age_group = c("65-79 years", "80-90 years"))
        prepare_analysis_factors <- function(d) {
            for (v in intersect(names(subgroup_levels), names(d))) {
                if (any(!is.na(d[[v]]) & !d[[v]] %in% subgroup_levels[[v]])) 
                  stop("Unexpected category in ", v)
                d[[v]] <- factor(d[[v]], levels = subgroup_levels[[v]])
            }
            d
        }
        format_p <- function(p) ifelse(is.na(p), "NE", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
        counts_for <- function(z, v) {
            z <- z[!is.na(z[[v]]) & !is.na(z$group), , drop = FALSE]
            n <- vapply(c("Control", "Intervention"), function(a) sum(z$group == a), 0L)
            e <- vapply(c("Control", "Intervention"), function(a) sum(z[[v]][z$group == a]), 0)
            setNames(c(e0 = e[1], n0 = n[1], e1 = e[2], n1 = n[2]), c("e0", "n0", "e1", "n1"))
        }
        row_status <- function(ct) {
            if (any(ct[c("n0", "n1")] == 0)) 
                return("Empty treatment arm")
            if (all(ct[c("e0", "e1")] == 0)) 
                return("No events in either arm")
            if (ct["e0"] == ct["n0"] && ct["e1"] == ct["n1"]) 
                return("All participants have event in both arms")
            "OK"
        }
        block_rule <- function(z, v, sg = NULL) {
            slices <- if (is.null(sg)) 
                list(z)
            else lapply(levels(droplevels(z[[sg]])), function(lv) z[!is.na(z[[sg]]) & z[[sg]] == lv, , drop = FALSE])
            ct <- lapply(slices, counts_for, v = v)
            trigger <- vapply(ct, function(x) row_status(x) == "OK" && xor(x["e0"] == 0, x["e1"] == 0), FALSE)
            list(firth = any(trigger), trigger_levels = if (is.null(sg)) if (any(trigger)) "Overall" else "" else paste(levels(droplevels(z[[sg]]))[trigger], 
                collapse = "; "))
        }
        empty_fit <- function(method, status, warnings = "") list(values = rep(NA_real_, 4), method = method, status = status, warnings = warnings, 
            df = NA_real_)
        ordinary_fit <- function(z, adjustment) {
            method <- if (adjustment == "IPTW") 
                "IPTW survey logistic; robust normal Wald 95% CI; residual-df t P"
            else "Ordinary logistic; normal Wald 95% CI; two-sided z P"
            if (adjustment == "IPTW" && any(!is.finite(z$IPTW) | z$IPTW <= 0)) 
                return(empty_fit(method, "Propensity model unavailable; weighted analysis NE"))
            if (adjustment == "Crude") {
                f <- safe(glm(.Y ~ group, data = z, family = binomial(), control = glm.control(maxit = 100, epsilon = 1e-12)))
            }
            else {
                design <- survey::svydesign(ids = ~1, weights = ~IPTW, data = z)
                f <- safe(survey::svyglm(.Y ~ group, design = design, family = quasibinomial(), control = glm.control(maxit = 100, 
                  epsilon = 1e-12)))
            }
            if (is.null(f$value) || !isTRUE(f$value$converged)) 
                return(empty_fit(method, "Ordinary model failed to converge", f$warnings))
            co <- summary(f$value)$coefficients
            if (!"groupIntervention" %in% rownames(co)) 
                return(empty_fit(method, "Treatment coefficient not estimable", f$warnings))
            b <- co["groupIntervention", 1]
            se <- co["groupIntervention", 2]
            vals <- c(exp(b), exp(b - qnorm(0.975) * se), exp(b + qnorm(0.975) * se), co["groupIntervention", 4])
            if (any(!is.finite(vals))) 
                return(empty_fit(method, "Nonfinite ordinary model estimate", f$warnings))
            list(values = vals, method = method, status = "OK", warnings = f$warnings, df = if (adjustment == "IPTW") f$value$df.residual else NA_real_)
        }
        firth_fit <- function(z, adjustment) {
            method <- if (adjustment == "IPTW") 
                "IPTW Firth logistic; full-cohort normalized case weights; profile CI/penalized LR"
            else "Firth logistic; profile CI/penalized LR"
            if (adjustment == "IPTW" && any(!is.finite(z$FIRTH_WEIGHT) | z$FIRTH_WEIGHT <= 0)) 
                return(empty_fit(method, "Propensity model unavailable; weighted analysis NE"))
            z$.W <- if (adjustment == "IPTW") 
                z$FIRTH_WEIGHT
            else rep(1, nrow(z))
            f <- safe(logistf(.Y ~ group, data = z, weights = .W, firth = TRUE, pl = TRUE, plconf = 2, control = ctrl, plcontrol = pctrl))
            if (is.null(f$value) || any(f$value$iter >= ctrl$maxit) || any(f$value$pl.iter >= pctrl$maxit) || grepl("converg|iterations.*exceeded", 
                f$warnings, ignore.case = TRUE)) 
                return(empty_fit(method, "Firth model/profile did not converge", f$warnings))
            term <- "groupIntervention"
            vals <- c(exp(f$value$coefficients[term]), exp(f$value$ci.lower[term]), exp(f$value$ci.upper[term]), f$value$prob[term])
            if (any(!is.finite(vals))) 
                return(empty_fit(method, "Nonfinite Firth estimate", f$warnings))
            list(values = vals, method = method, status = "OK", warnings = f$warnings, df = NA_real_)
        }
        effect <- function(z, v, population, adjustment = "Crude", sg = "Overall", lv = "All", use_firth = NULL, block_id = "", trigger_levels = "") {
            z <- z[!is.na(z[[v]]) & !is.na(z$group), , drop = FALSE]
            z$.Y <- z[[v]]
            ct <- counts_for(z, v)
            status <- row_status(ct)
            if (is.null(use_firth)) 
                use_firth <- block_rule(z, v)$firth
            selected <- empty_fit("Not estimated", status)
            if (status == "OK") {
                selected <- if (use_firth) firth_fit(z, adjustment) else ordinary_fit(z, adjustment)
            }
            vals <- selected$values
            data.frame(population = population, adjustment = adjustment, outcome = v, outcome_label = if (v %in% names(outlabs)) 
                unname(outlabs[[v]])
            else v, subgroup_var = sg, subgroup_level = lv, event1 = unname(ct["e1"]), n1 = unname(ct["n1"]), event0 = unname(ct["e0"]), 
                n0 = unname(ct["n0"]), OR = vals[1], lower = vals[2], upper = vals[3], p = vals[4], method = selected$method, status = selected$status, 
                interaction_p = NA_real_, interaction_method = "", block_id = block_id, 
                block_firth = use_firth, block_trigger_levels = trigger_levels, p_df = selected$df, warnings = selected$warnings, stringsAsFactors = FALSE, 
                row.names = NULL)
        }
        ordinary_interaction <- function(z, adjustment) {
            method <- if (adjustment == "IPTW") 
                "IPTW survey logistic interaction; joint Wald F test"
            else "Ordinary logistic interaction; likelihood-ratio test"
            if (adjustment == "IPTW" && any(!is.finite(z$IPTW) | z$IPTW <= 0)) 
                return(list(p = NA_real_, statistic = NA_real_, df = NA_real_, denominator_df = NA_real_, method = method, status = "Propensity model unavailable; weighted analysis NE", 
                  warnings = ""))
            if (adjustment == "IPTW") {
                nt <- table(z$group, z$.SG)
                ev <- tapply(z$.Y, list(z$group, z$.SG), sum)
                if (any(ev == 0 | ev == nt)) 
                  return(list(p = NA_real_, statistic = NA_real_, df = NA_real_, denominator_df = NA_real_, method = method, status = "Survey Wald interaction NE under separation", 
                    warnings = ""))
            }
            if (adjustment == "Crude") {
                f <- safe(glm(.Y ~ group * .SG, data = z, family = binomial(), control = glm.control(maxit = 100, epsilon = 1e-12)))
                n <- safe(glm(.Y ~ group + .SG, data = z, family = binomial(), control = glm.control(maxit = 100, epsilon = 1e-12)))
                warn <- paste(f$warnings, n$warnings)
                if (is.null(f$value) || is.null(n$value) || !f$value$converged || !n$value$converged) 
                  return(list(p = NA_real_, statistic = NA_real_, df = NA_real_, denominator_df = NA_real_, method = method, status = "Ordinary interaction fit failed", 
                    warnings = warn))
                statistic <- n$value$deviance - f$value$deviance
                df <- f$value$rank - n$value$rank
                p <- pchisq(max(0, statistic), df, lower.tail = FALSE)
                denominator_df <- NA_real_
            }
            else {
                design <- survey::svydesign(ids = ~1, weights = ~IPTW, data = z)
                f <- safe(survey::svyglm(.Y ~ group * .SG, design = design, family = quasibinomial(), control = glm.control(maxit = 100, 
                  epsilon = 1e-12)))
                if (is.null(f$value) || !isTRUE(f$value$converged)) 
                  return(list(p = NA_real_, statistic = NA_real_, df = NA_real_, denominator_df = NA_real_, method = method, status = "Survey interaction fit failed", 
                    warnings = f$warnings))
                test <- safe(survey::regTermTest(f$value, ~group:.SG, method = "Wald"))
                warn <- paste(f$warnings, test$warnings)
                if (is.null(test$value)) 
                  return(list(p = NA_real_, statistic = NA_real_, df = NA_real_, denominator_df = NA_real_, method = method, status = "Survey Wald F test failed", 
                    warnings = warn))
                statistic <- unname(test$value$Ftest)
                df <- test$value$df
                denominator_df <- test$value$ddf
                p <- unname(test$value$p)
            }
            list(p = p, statistic = statistic, df = df, denominator_df = denominator_df, method = method, status = if (is.finite(p)) "OK" else "Nonfinite interaction P", 
                warnings = warn)
        }
        firth_interaction <- function(z, adjustment) {
            method <- if (adjustment == "IPTW") 
                "IPTW Firth interaction; constrained penalized likelihood-ratio test; full-cohort normalized weights"
            else "Firth interaction; constrained penalized likelihood-ratio test"
            fail <- function(status, warnings = "") list(p = NA_real_, statistic = NA_real_, df = NA_real_, denominator_df = NA_real_, 
                method = method, status = status, warnings = warnings)
            if (adjustment == "IPTW" && any(!is.finite(z$FIRTH_WEIGHT) | z$FIRTH_WEIGHT <= 0)) 
                return(fail("Propensity model unavailable; weighted analysis NE"))
            z$.W <- if (adjustment == "IPTW") 
                z$FIRTH_WEIGHT
            else rep(1, nrow(z))
            f <- safe(logistf(.Y ~ group * .SG, data = z, weights = .W, firth = TRUE, pl = FALSE, model = TRUE, control = ctrl))
            if (is.null(f$value) || any(f$value$iter >= ctrl$maxit)) 
                return(fail("Firth interaction fit failed", f$warnings))
            idx <- which(grepl(":", names(f$value$coefficients), fixed = TRUE))
            tt <- safe(logistftest(f$value, test = idx, values = rep(0, length(idx)), control = ctrl))
            warn <- paste(f$warnings, tt$warnings)
            if (is.null(tt$value) || grepl("converg|iterations.*exceeded", warn, ignore.case = TRUE)) 
                return(fail("Firth constrained interaction test failed", warn))
            statistic <- 2 * unname(tt$value$loglik["full"] - tt$value$loglik["null"])
            df <- length(idx)
            p <- pchisq(max(0, statistic), df, lower.tail = FALSE)
            list(p = p, statistic = statistic, df = df, denominator_df = NA_real_, method = method, status = if (is.finite(p)) "OK" else "Nonfinite interaction P", 
                warnings = warn)
        }
        interaction <- function(z, v, sg, population, adjustment = "Crude", use_firth = NULL, block_id = "", trigger_levels = "") {
            z <- z[!is.na(z[[v]]) & !is.na(z[[sg]]) & !is.na(z$group), , drop = FALSE]
            z$.Y <- z[[v]]
            z$.SG <- droplevels(z[[sg]])
            tab <- table(z$group, z$.SG)
            status <- "OK"
            if (is.null(use_firth)) 
                use_firth <- block_rule(z, v, sg)$firth
            if (nlevels(z$.SG) < 2) 
                status <- "Fewer than two subgroup levels"
            else if (any(tab == 0)) 
                status <- "Empty treatment-by-subgroup cell"
            else if (sum(z$.Y) == 0) 
                status <- "No events across subgroup block"
            if (status == "OK") {
                event_cells <- tapply(z$.Y, list(z$group, z$.SG), sum)
                if (any(colSums(event_cells) == 0)) 
                  status <- "No events in either arm within a subgroup"
                else if (!use_firth && any(event_cells == 0 | event_cells == tab)) 
                  status <- "Ordinary interaction NE: nonfinite MLE in treatment-by-subgroup outcome cells"
            }
            fit <- list(p = NA_real_, statistic = NA_real_, df = NA_real_, denominator_df = NA_real_, method = "Not estimated", status = status, 
                warnings = "")
            if (status == "OK") {
                fit <- if (use_firth) firth_interaction(z, adjustment) else ordinary_interaction(z, adjustment)
            }
            data.frame(population = population, adjustment = adjustment, outcome = v, subgroup_var = sg, p = fit$p, method = fit$method, 
                status = fit$status, LR = if (grepl("Wald F", fit$method)) 
                  NA_real_
                else fit$statistic, F = if (grepl("Wald F", fit$method)) 
                  fit$statistic
                else NA_real_, df = fit$df, denominator_df = fit$denominator_df, block_id = block_id, block_firth = use_firth,
                block_trigger_levels = trigger_levels, warnings = fit$warnings, 
                stringsAsFactors = FALSE)
        }
    }, envir)
}

precise_publication_labels <- function(envir = parent.frame()) {
    evalq({
        FLAGGED_LABEL <- "Flagged by the AI-ECG screening algorithm"
        NOT_FLAGGED_LABEL <- "Not flagged by the AI-ECG screening algorithm"
        PROTOCOL_HIGH_LABEL <- "Protocol defined AI-ECG high risk"
        PROTOCOL_LOW_LABEL <- "Protocol defined AI-ECG low risk"
        THREE_GROUP_LABELS <- c(FLAGGED_LABEL, paste0(NOT_FLAGGED_LABEL, " (", PROTOCOL_HIGH_LABEL, ")"), paste0(NOT_FLAGGED_LABEL, 
            " (", PROTOCOL_LOW_LABEL, ")"))
        screening_label <- function(x) {
            y <- as.character(x)
            y[y %in% c("High risk", "High")] <- FLAGGED_LABEL
            y[y %in% c("Low risk", "Low")] <- NOT_FLAGGED_LABEL
            y
        }
        wrap_publication_label <- function(s, width_mm = 46, size = 5.6) {
            words <- strsplit(s, "[[:space:]]+")[[1]]
            lines <- character()
            cur <- ""
            for (w in words) {
                z <- if (nzchar(cur)) 
                  paste(cur, w)
                else w
                if (nzchar(cur) && strwidth(z, units = "inches", cex = size/12) * 25.4 > width_mm) {
                  lines <- c(lines, cur)
                  cur <- w
                }
                else cur <- z
            }
            paste(c(lines, cur), collapse = "\n")
        }
    }, envir)
}

precise_screening_labels <- function(envir = parent.frame()) {
    evalq({
        screening_label <- function(x) {
            map <- c(`High risk` = "Flagged by the AI-ECG screening algorithm", `Low risk` = "Not flagged by the AI-ECG screening algorithm", 
                All = "Overall", Overall = "Overall")
            ans <- unname(map[as.character(x)])
            ans[is.na(ans)] <- as.character(x)[is.na(ans)]
            ans
        }
    }, envir)
}

precise_table2_map <- function(envir = parent.frame()) {
    evalq({
        table2_display_map <- function(tab) {
            canonical <- function(s) tolower(gsub("[[:space:]]+", " ", trimws(s)))
            rows <- vapply(tab$rows, function(r) r$row, 0L)
            labels <- vapply(tab$rows, function(r) r$cells[[1]]$text, "")
            find_row <- function(label, parent = NULL) {
                hit <- which(canonical(labels) == canonical(label))
                if (!is.null(parent)) {
                  p <- which(canonical(labels) == canonical(parent))
                  if (length(p) != 1L) 
                    stop("Ambiguous Table 2 parent: ", parent)
                  hit <- hit[hit > p]
                  if (length(hit)) 
                    hit <- hit[1L]
                }
                if (length(hit) != 1L) 
                  stop("Missing or ambiguous Table 2 row: ", label)
                rows[hit]
            }
            map <- data.frame(row_index = integer(), variable = character(), level = character(), show_smd = logical(), stringsAsFactors = FALSE)
            add <- function(label, variable, level = "1", show_smd = FALSE, parent = NULL) {
                map <<- rbind(map, data.frame(row_index = find_row(label, parent), variable = variable, level = level, show_smd = show_smd, 
                  stringsAsFactors = FALSE))
            }
            add("Flagged by the AI-ECG screening algorithm, n (%)", "ai_ecg_screening_flag", "1")
            add("Not flagged by the AI-ECG screening algorithm, n (%)", "ai_ecg_screening_flag", "0")
            add("Male, n (%)", "male")
            add("Age (years, mean ± SD)", "age_years", "Continuous")
            add("65-79 years, n (%)", "age_ge80", "0")
            add("80-90 years, n (%)", "age_ge80", "1")
            clinical <- c("Atrial fibrillation/flutter, n (%)", "Hypertension, n (%)", "Diabetes mellitus, n (%)", "Heart failure, n (%)", 
                "Chronic kidney disease, n (%)", "Chronic obstructive pulmonary disease, n (%)", "Peripheral arterial disease, n (%)", 
                "Valvular heart disease, n (%)", "Cerebrovascular disease, n (%)", "Obstructive sleep apnea, n (%)", "Thyroid disease, n (%)", 
                "Prior myocardial infarction, n (%)", "History of coronary artery bypass grafting, n (%)", "History of electrophysiology study or catheter ablation, n (%)")
            for (i in seq_along(clinical)) add(clinical[i], precise_propensity_variables[i], show_smd = TRUE)
            add("Mobitz type I AV block in 12-lead ECG or Holter, n (%)", "mobitz_i_av_block")
            add("Borderline low heart rate in 12-lead ECG or Holter, n (%)", "borderline_low_heart_rate")
            for (q in list(c("Complete right bundle branch block in 12-lead ECG", "rbbb_recency"), c("Complete left bundle branch block in 12-lead ECG", 
                "lbbb_recency"), c("Premature ventricular complex in 12-lead ECG", "pvc_ecg_recency"))) {
                add("Recent-onset (≤1 years), n (%)", q[2], "Recent-onset", parent = q[1])
                add("Long-standing (>1 years), n (%)", q[2], "Long-standing", parent = q[1])
            }
            add("Left fascicular block in 12-lead ECG, n (%)", "left_fascicular_block")
            add("Premature ventricular complex burden 5.0-9.9%, n (%)", "holter_pvc_burden_category", "5.0-9.9%")
            add("Premature ventricular complex burden ≥10.0%, n (%)", "holter_pvc_burden_category", ">=10.0%")
            add("Non-sustained ventricular tachycardia in Holter, n (%)", "holter_nsvt")
            add("Borderline pause or asymptomatic pause in Holter, n (%)", "holter_borderline_pause")
            add("Marked PR prolongation (≥250) in 12-lead ECG, n (%)", "ecg_marked_pr_prolongation")
            add("Mild PR prolongation (200-249) in 12-lead ECG, n (%)", "ecg_mild_pr_prolongation")
            add("Documeted PR prolongation (≥200) in Holter, n (%)", "holter_pr_prolongation")
            meds <- c("Beta blocker, n (%)", "Non-dihydropyridine calcium channel blockers, n (%)", "Class IA antiarrhythmia drug, n (%)", 
                "Class IC antiarrhythmia drug, n (%)", "Class III antiarrhythmia drug, n (%)", "Digoxin, n (%)")
            for (i in seq_along(meds)) add(meds[i], precise_propensity_variables[i + 14L], show_smd = TRUE)
            map <- map[order(map$row_index), ]
            stopifnot(!anyDuplicated(map$row_index), nrow(map) == 42L, sum(map$show_smd) == 20L)
            map
        }
    }, envir)
}

precise_table3 <- function(envir = parent.frame()) {
    evalq({
        precise_core_helpers()
        vars <- precise_exploratory_outcomes
        binary <- c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
        labels <- c("New-onset atrial fibrillation", "Average heart rate (bpm)", "Maximal heart rate (bpm)", "Minimal heart rate (bpm)", 
            "Longest RR interval (seconds)", "Premature atrial complex burden (%)", "Premature ventricular complex burden (%)", "New-onset non-sustained ventricular tachycardia", 
            "New-onset bundle branch block", "New-onset PR prolongation", "New-onset Mobitz type I AV block")
        if (exists("outlabs", inherits = FALSE)) 
            outlabs <- c(outlabs, setNames(labels, vars))
        else outlabs <- setNames(labels, vars)
        for (adj in c("Crude", "IPTW")) for (arm in c("Intervention", "Control")) {
            ci <- if (adj == "Crude") 
                1L
            else 5L
            ci <- ci + as.integer(arm == "Control")
            arm_weights <- pp$IPTW[pp$group == arm]
            header_text <- if (adj == "Crude")
                sprintf("%s\n(N = %s)", arm, fmtn(sum(pp$group == arm)))
            else sprintf("%s\n(ESS = %.1f)", arm, sum(arm_weights)^2/sum(arm_weights^2))
            put("main", 3, 1, ci, list(text = header_text, n_total = sum(pp$group == arm),
                arm = arm, adjustment = adj), "header")
        }
        for (i in seq_along(vars)) {
            v <- vars[i]
            ri <- i + 1L
            z <- pp[!is.na(pp[[v]]), , drop = FALSE]
            rule <- if (binary[i]) 
                block_rule(z, v)
            else {
            }
            for (adj in c("Crude", "IPTW")) {
                pv <- NA_real_
                method <- ""
                status <- "OK"
                if (binary[i]) {
                  ee <- effect(z, v, "Monitored", adj, use_firth = rule$firth, block_id = paste("Table3", v, sep = ":"), trigger_levels = rule$trigger_levels)
                  pv <- ee$p
                  method <- ee$method
                  status <- ee$status
                }
                else if (adj == "Crude") {
                  method <- "Two-sided Welch t test"
                  tt <- tryCatch(t.test(z[[v]] ~ z$group), error = function(e) {
                  })
                  if (!is.null(tt)) 
                    pv <- tt$p.value
                  else status <- "Continuous comparison not estimable"
                }
                else {
                  method <- "Treatment-only survey-weighted linear regression; robust variance; residual-df t test"
                  if (any(!is.finite(z$IPTW) | z$IPTW <= 0) || !nrow(z)) {
                    status <- "Invalid propensity estimation; weighted analysis NE"
                  }
                  else {
                    z$.Y <- z[[v]]
                    design <- survey::svydesign(ids = ~1, weights = ~IPTW, data = z)
                    fit <- tryCatch(survey::svyglm(.Y ~ group, design = design, family = gaussian()), error = function(e) {
                    })
                    if (is.null(fit)) {
                      status <- "Weighted linear model failed"
                    }
                    else {
                      co <- summary(fit)$coefficients
                      pv <- co["groupIntervention", 4]
                    }
                  }
                }
                for (arm in c("Intervention", "Control")) {
                  ci <- if (adj == "Crude") 
                    1L
                  else 5L
                  ci <- ci + as.integer(arm == "Control")
                  za <- z[z$group == arm, , drop = FALSE]
                  digits <- if (i %in% c(5, 6, 7)) 
                    2L
                  else 1L
                  val <- summary_cell(za, v, "Yes", adj, digits = digits, show_n = !binary[i], binary_n = binary[i])
                  val$arm <- arm
                  val$population <- "Monitored"
                  val$cohort_n <- sum(pp$group == arm)
                  val$excluded_or_unavailable_n <- val$cohort_n - nrow(za)
                  put("main", 3, ri, ci, val)
                }
                put("main", 3, ri, if (adj == "Crude") 
                  3L
                else 7L, list(text = fmtp(pv), p = pv, method = method, status = status, variable = v, adjustment = adj), "p")
            }
        }
    }, envir)
}

precise_core_analysis <- function(root) {
    out <- file.path(root, "analysis")
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    precise_core_helpers()
    write_csv <- function(x, nm) write.csv(x, file.path(out, paste0(nm, ".csv")), row.names = FALSE, na = "")
    d <- read.csv(file.path(root, "data", "PRECISE-AI.csv"), check.names = FALSE, stringsAsFactors = FALSE, na.strings = "")
    d <- prepare_analysis_factors(d)
    stopifnot(!anyDuplicated(d$row_number), all(vapply(d[outs], function(x) all(!is.na(x) & x %in% 0:1), FALSE)))
    ps_variables <- readLines(file.path(root, "data", "ps_variables.txt"))
    stopifnot(identical(ps_variables, precise_propensity_variables), all(ps_variables %in% names(d)))
    pp <- d[d$selected_matched_cohort == 1, , drop = FALSE]
    complete_ps <- complete.cases(pp[ps_variables])
    ps_formula <- reformulate(ps_variables, response = "treatment")
    ps_ok <- FALSE
    pp$PS <- pp$IPTW <- pp$FIRTH_WEIGHT <- NA_real_
    if (all(complete_ps)) {
        mf <- model.frame(ps_formula, data = pp, na.action = na.fail, drop.unused.levels = FALSE)
        fullx <- model.matrix(ps_formula, data = mf)
        qx <- qr(fullx, tol = 1e-10)
        keep <- sort(qx$pivot[seq_len(qx$rank)])
        design_matrix_summary <- data.frame(column = colnames(fullx), used = seq_len(ncol(fullx)) %in% keep, reason = ifelse(seq_len(ncol(fullx)) %in% 
            keep, "Identifiable design column", "Constant or exactly linearly dependent design column"), stringsAsFactors = FALSE)
        write_csv(design_matrix_summary, "propensity_design")
        X <- fullx[, keep, drop = FALSE]
        fitted <- safe(glm.fit(x = X, y = pp$treatment, family = binomial(), control = glm.control(maxit = 1000, epsilon = 1e-12)))
        ps_fit <- fitted$value
        ps_warn <- fitted$warnings
        ps_ok <- !is.null(ps_fit) && isTRUE(ps_fit$converged) && !isTRUE(ps_fit$boundary) && all(is.finite(ps_fit$coefficients)) && 
            all(is.finite(ps_fit$fitted.values)) && all(ps_fit$fitted.values > 0 & ps_fit$fitted.values < 1) && !grepl("numerically 0 or 1|did not converge|separat", 
            ps_warn, ignore.case = TRUE)
        if (ps_ok) {
            pp$PS <- as.numeric(ps_fit$fitted.values)
            pp$IPTW <- ifelse(pp$treatment == 1, 1/pp$PS, 1/(1 - pp$PS))
            pp$FIRTH_WEIGHT <- nrow(pp) * pp$IPTW/sum(pp$IPTW)
            ps_ok <- all(is.finite(pp$IPTW) & pp$IPTW > 0) & all(is.finite(pp$FIRTH_WEIGHT) & pp$FIRTH_WEIGHT > 0)
        }
        if (!ps_ok) 
            pp$PS <- pp$IPTW <- pp$FIRTH_WEIGHT <- NA_real_
    }
    weights <- data.frame(row_number = pp$row_number, PS = pp$PS, IPTW = pp$IPTW, FIRTH_WEIGHT = pp$FIRTH_WEIGHT)
    write_csv(weights, "propensity_weights")
    saveRDS(list(pp = pp, ps_numeric_usable = ps_ok), file.path(out, "analysis_objects.rds"))
    effects <- ints <- list()
    add_effect <- function(e) {
        effects[[length(effects) + 1]] <<- e
    }
    add_interaction <- function(it) {
        ints[[length(ints) + 1]] <<- it
    }
    attach_interaction <- function(e, it) {
        e$interaction_p <- it$p
        e$interaction_method <- it$method
        e
    }
    for (v in c("cied_implantation", "sick_sinus_syndrome")) {
        for (pop in c("Overall", "High risk")) for (sg in c("clinical_stratum", "clinical_holter_stratum", "age_group", "sex_group")) {
            zz <- if (pop == "Overall") 
                d
            else d[d$screening_group == "High risk", ]
            rule <- block_rule(zz, v, sg)
            id <- paste("ITT", pop, v, sg, sep = ":")
            it <- interaction(zz, v, sg, pop, use_firth = rule$firth, block_id = id, trigger_levels = rule$trigger_levels)
            add_interaction(it)
            for (lv in levels(d[[sg]])) add_effect(attach_interaction(effect(zz[zz[[sg]] == lv, ], v, pop, sg = sg, lv = lv, use_firth = rule$firth, 
                block_id = id, trigger_levels = rule$trigger_levels), it))
        }
    }
    ef <- do.call(rbind, effects)
    inter <- do.call(rbind, ints)
    fmt_event <- function(e, n) if (n == 0) 
        "0/0 (NE)"
    else sprintf("%d/%d (%.1f%%)", e, n, 100 * e/n)
    fmt_effect <- function(v) if (is.na(v[1])) 
        "NE"
    else sprintf("%.2f (%.2f-%.2f)", v[1], v[2], v[3])
    ef$intervention_display <- mapply(fmt_event, ef$event1, ef$n1)
    ef$control_display <- mapply(fmt_event, ef$event0, ef$n0)
    ef$OR_CI <- apply(ef[, c("OR", "lower", "upper")], 1, fmt_effect)
    ef$p_display <- format_p(ef$p)
    ef$interaction_p_display <- format_p(ef$interaction_p)
    write_csv(ef, "effects")
    write_csv(inter, "interactions")
}

precise_core_diagnostics <- function(root) {
    a <- file.path(root, "analysis")
    t <- file.path(root, "tables")
    dir.create(t, showWarnings = FALSE)
    obj <- readRDS(file.path(a, "analysis_objects.rds"))
    pp <- obj$pp
    ps_variables <- readLines(file.path(root, "data", "ps_variables.txt"))
    design <- read.csv(file.path(a, "propensity_design.csv"), check.names = FALSE)
    stopifnot(nrow(pp) == 240L, identical(ps_variables, precise_propensity_variables), !anyNA(pp[ps_variables]), all(is.finite(pp$PS)), 
        all(pp$PS > 0 & pp$PS < 1), all(is.finite(pp$IPTW) & pp$IPTW > 0))
    stopifnot(max(abs(pp$IPTW - ifelse(pp$treatment == 1, 1/pp$PS, 1/(1 - pp$PS)))) < 1e-12, max(abs(pp$FIRTH_WEIGHT - nrow(pp) * pp$IPTW/sum(pp$IPTW))) < 
        1e-12)
    n1 <- sum(pp$treatment == 1)
    n0 <- sum(pp$treatment == 0)
    auc <- (sum(rank(pp$PS, ties.method = "average")[pp$treatment == 1]) - n1 * (n1 + 1)/2)/(n1 * n0)
    brier <- mean((pp$treatment - pp$PS)^2)
    balance <- do.call(rbind, lapply(ps_variables, function(v) {
        x <- pp[[v]]
        ii <- pp$treatment == 1
        m1 <- mean(x[ii])
        m0 <- mean(x[!ii])
        wm1 <- weighted.mean(x[ii], pp$IPTW[ii])
        wm0 <- weighted.mean(x[!ii], pp$IPTW[!ii])
        var1 <- m1 * (1 - m1)
        var0 <- m0 * (1 - m0)
        denom <- sqrt((var1 + var0)/2)
        calc <- function(delta) if (denom == 0) {
            if (abs(delta) < 1e-12) 
                0
            else Inf
        }
        else abs(delta)/denom
        data.frame(smd_before = calc(m1 - m0), smd_after = calc(wm1 - wm0))
    }))
    converged_numeric <- isTRUE(obj$ps_numeric_usable)
    model_values <- c(auc, brier, sum(design$used), length(ps_variables), as.integer(converged_numeric), max(balance$smd_before), max(balance$smd_after), 
        sum(balance$smd_before >= 0.1), sum(balance$smd_after >= 0.1), 0, 100 * max(pp$IPTW)/sum(pp$IPTW))
    model_labels <- c("In-sample c statistic", "In-sample Brier score", "Model rank including intercept", "Nominal candidate baseline indicators", 
        "Numerical IRLS convergence", "Max modeled-covariate SMD before", "Max modeled-covariate SMD after", "Modeled covariates SMD >=0.10 before", 
        "Modeled covariates SMD >=0.10 after", "Propensity values clipped", "Largest weight as percentage of total")
    model_display <- c(sprintf("%.3f", model_values[1:2]), as.character(model_values[3:5]),
        sprintf("%.3f", model_values[6:7]), as.character(model_values[8:10]), sprintf("%.2f%%", model_values[11]))
    model_display_table <- data.frame(metric = model_labels, display = model_display, stringsAsFactors = FALSE)
    cohorts <- list(Overall = pp, Intervention = pp[pp$treatment == 1, , drop = FALSE], Control = pp[pp$treatment == 0, , drop = FALSE])
    distribution <- do.call(rbind, lapply(names(cohorts), function(arm) {
        z <- cohorts[[arm]]
        w <- z$IPTW
        nw <- z$FIRTH_WEIGHT
        data.frame(cohort = arm, participants = nrow(z), minimum_propensity = min(z$PS), median_propensity = median(z$PS), maximum_propensity = max(z$PS), 
            minimum_raw_weight = min(w), median_raw_weight = median(w), p95_raw_weight = unname(quantile(w, 0.95, type = 7)), p99_raw_weight = unname(quantile(w, 
                0.99, type = 7)), maximum_raw_weight = max(w), sum_raw_weights = sum(w), minimum_normalized_weight = min(nw), maximum_normalized_weight = max(nw), 
            sum_normalized_weights = sum(nw), effective_sample_size = sum(w)^2/sum(w^2), stringsAsFactors = FALSE)
    }))
    distribution_labels <- c("Participants", "Minimum propensity", "Median propensity", "Maximum propensity", "Minimum raw weight", 
        "Median raw weight", "95th percentile raw weight", "99th percentile raw weight", "Maximum raw weight", "Sum raw weights", "Minimum normalized weight", 
        "Maximum normalized weight", "Sum normalized weights", "Effective sample size")
    dist_table <- data.frame(metric = distribution_labels, stringsAsFactors = FALSE)
    for (arm in names(cohorts)) {
        z <- distribution[distribution$cohort == arm, -1]
        val <- as.numeric(z[1, ])
        dist_table[[arm]] <- c(as.character(val[1]), sprintf("%.3f", val[-1]))
    }
    write_csv <- function(x, file) write.csv(x, file, row.names = FALSE, na = "")
    write_csv(model_display_table, file.path(t, "Supplemental_Table_01_Propensity_Model_Diagnostics.csv"))
    write_csv(dist_table, file.path(t, "Supplemental_Table_02_Propensity_Distribution.csv"))
}

precise_core_yields <- function(root) {
    out <- file.path(root, "analysis")
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    d <- read.csv(file.path(root, "data", "PRECISE-AI.csv"), check.names = FALSE, stringsAsFactors = FALSE, na.strings = "")
    z <- d[d$selected_matched_cohort == 1 & d$group == "Intervention", , drop = FALSE]
    levels_by <- list(clinical_stratum = c("I", "II", "III"), clinical_holter_stratum = c("I-a", "I-b", "II-a", "II-b", "III-a", "III-b"), 
        age_group = c("65-79 years", "80-90 years"), sex_group = c("Female", "Male"))
    block_labels <- c(clinical_stratum = "Prespecified clinical strata", clinical_holter_stratum = "Strata by recent Holter monitoring", 
        age_group = "Age group", sex_group = "Recorded sex")
    subgroup_labels <- list(clinical_stratum = c(I = "Potentially suboptimal rhythm evaluation (Stratum I)", II = "Risk indicators supporting rhythm monitoring (Stratum II)", 
        III = "No prespecified rhythm risk indicators (Stratum III)"), clinical_holter_stratum = c(`I-a` = "Stratum I (without Holter within 1 year)", 
        `I-b` = "Stratum I (with Holter within 1 year)", `II-a` = "Stratum II (without Holter within 1 year)", `II-b` = "Stratum II (with Holter within 1 year)", 
        `III-a` = "Stratum III (without Holter within 1 year)", `III-b` = "Stratum III (with Holter within 1 year)"))
    label_for <- function(sg, lv) if (!is.null(subgroup_labels[[sg]]) && lv %in% names(subgroup_labels[[sg]])) 
        unname(subgroup_labels[[sg]][[lv]])
    else lv
    wilson <- function(e, n) {
        if (n == 0) 
            return(c(NA_real_, NA_real_))
        q <- qnorm(0.975)
        p <- e/n
        den <- 1 + q * q/n
        mid <- (p + q * q/(2 * n))/den
        half <- q/den * sqrt(p * (1 - p)/n + q * q/(4 * n * n))
        c(max(0, mid - half), min(1, mid + half))
    }
    fmtp <- function(p) if (is.na(p)) 
        "NE"
    else if (p < 0.001) 
        "<0.001"
    else sprintf("%.3f", p)
    safe_fisher <- function(tab) tryCatch(fisher.test(tab, alternative = "two.sided", workspace = 2e+07)$p.value, error = function(e) NA_real_)
    outputs <- global <- list()
    for (v in c("cied_indication", "cied_implantation")) {
        stopifnot(all(!is.na(z[[v]]) & z[[v]] %in% 0:1))
        figure <- if (v == "cied_indication") 
            "indication"
        else "implantation"
        outcome_label <- if (v == "cied_indication") 
            "CIED indication"
        else "CIED implantation"
        order <- 0L
        addrow <- function(sg, lv, zz) {
            e <- sum(zz[[v]])
            n <- nrow(zz)
            ci <- wilson(e, n)
            status <- if (n == 0) 
                "No participants: proportion and confidence interval not estimable"
            else "OK"
            order <<- order + 1L
            outputs[[length(outputs) + 1L]] <<- data.frame(figure = figure, outcome = v, outcome_label = outcome_label, block = if (sg == 
                "Overall") 
                "Overall"
            else block_labels[[sg]], subgroup_var = sg, level = lv, label = if (sg == "Overall") 
                "Overall"
            else label_for(sg, lv), events = e, total = n, proportion = if (n == 0) 
                NA_real_
            else e/n, lower = ci[1], upper = ci[2], percent = if (n == 0) 
                NA_real_
            else 100 * e/n, lower_percent = 100 * ci[1], upper_percent = 100 * ci[2], order = order, status = status, counts_display = sprintf("%d/%d", 
                e, n), yield_display = if (n == 0) 
                "NE"
            else sprintf("%.1f (%.1f-%.1f)", 100 * e/n, 100 * ci[1], 100 * ci[2]), stringsAsFactors = FALSE)
        }
        addrow("Overall", "All", z)
        for (sg in names(levels_by)) {
            stopifnot(all(z[[sg]] %in% levels_by[[sg]]))
            for (lv in levels_by[[sg]]) addrow(sg, lv, z[z[[sg]] == lv, , drop = FALSE])
            n <- vapply(levels_by[[sg]], function(lv) sum(z[[sg]] == lv), 0L)
            e <- vapply(levels_by[[sg]], function(lv) sum(z[[v]][z[[sg]] == lv]), 0)
            keep <- n > 0
            gp <- if (sum(keep) >= 2) 
                safe_fisher(rbind(e[keep], n[keep] - e[keep]))
            else NA_real_
            test <- if (sum(keep) == 2) 
                "Two-sided global Fisher exact test (2 x 2)"
            else "Two-sided global Fisher-Freeman-Halton exact test (2 x K)"
            global[[length(global) + 1L]] <- data.frame(figure = figure, outcome = v, outcome_label = outcome_label, block = block_labels[[sg]], 
                subgroup_var = sg, test = test, p_value = gp, p_display = fmtp(gp), n_levels_tested = sum(keep), tested_levels = paste(levels_by[[sg]][keep], 
                  collapse = ";"), excluded_empty_levels = paste(levels_by[[sg]][!keep], collapse = ";"), status = if (sum(keep) < 
                  2) 
                  "Fewer than two nonempty categories"
                else if (is.na(gp)) 
                  "Exact test failed"
                else "OK")
            stopifnot(sum(n) == nrow(z), sum(e) == sum(z[[v]]))
        }
    }
    res <- do.call(rbind, outputs)
    g <- do.call(rbind, global)
    stopifnot(nrow(res) == 28, all(table(res$figure) == 14), all(res$events <= res$total), all(is.na(res$lower[res$total == 0])), nrow(g) == 
        8L, all(table(g$figure) == 4L), !anyDuplicated(g[c("figure", "subgroup_var")]))
    stopifnot(all(res$lower[res$total > 0] <= res$proportion[res$total > 0] + 1e-12), all(res$upper[res$total > 0] >= res$proportion[res$total > 
        0] - 1e-12))
    write.csv(res, file.path(out, "pp_intervention_yields.csv"), row.names = FALSE, na = "")
    write.csv(g, file.path(out, "pp_intervention_global_tests.csv"), row.names = FALSE, na = "")
}

precise_core_pathways <- function(root) {
    portable_pdf <- function(...) {
        if (capabilities("cairo")) 
            grDevices::cairo_pdf(...)
        else grDevices::pdf(...)
    }
    suppressPackageStartupMessages(library(jsonlite))
    precise_publication_labels()
    out <- file.path(root, "figures")
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    d <- read.csv(file.path(root, "data/PRECISE-AI.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    d$arm <- d$group
    stopifnot(all(d$ai_ecg_screening_flag %in% 0:1), identical(d$screening_group, ifelse(d$ai_ecg_screening_flag == 1, "High risk", 
        "Low risk")))
    d$risk <- ifelse(d$ai_ecg_screening_flag == 1, "High", "Low")
    d$route <- ifelse(d$selected_matched_cohort == 1, "Trial-arranged ambulatory ECG monitoring", ifelse(!is.na(d[["recorded_holter_count_180d"]]) & d[["recorded_holter_count_180d"]] > 
        0, "Usual-care Holter", "No recorded Holter"))
    d$indication <- d$cied_indication == 1
    stopifnot(all(d$indication[d$cied_implantation == 1]))
    wilson <- function(e, n) {
        if (n == 0) 
            return(c(NA_real_, NA_real_))
        z <- qnorm(0.975)
        p <- e/n
        den <- 1 + z * z/n
        m <- (p + z * z/(2 * n))/den
        h <- z/den * sqrt(p * (1 - p)/n + z * z/(4 * n * n))
        100 * c(max(0, m - h), min(1, m + h))
    }
    node_rows <- list()
    edge_rows <- list()
    addnode <- function(id, parent, arm, risk, route, stage, z, x, y, label) {
        n <- nrow(z)
        e <- sum(z$cied_implantation == 1)
        ci <- wilson(e, n)
        node_rows[[length(node_rows) + 1L]] <<- data.frame(node_id = id, parent_id = parent, arm = arm, risk = risk, route = route, 
            stage = stage, label = label, n = n, implantations = e, implantation_percent = if (n > 0) 
                100 * e/n
            else NA_real_, lower_percent = ci[1], upper_percent = ci[2], x = x, y = y, stringsAsFactors = FALSE)
    }
    route_order <- c("Trial-arranged ambulatory ECG monitoring", "Usual-care Holter", "No recorded Holter")
    for (arm in c("Intervention", "Control")) {
        a <- d[d$arm == arm, , drop = FALSE]
        prefix <- tolower(arm)
        xs <- if (arm == "Intervention") 
            c(23, 110, 200, 285)
        else c(617, 530, 440, 355)
        allid <- paste0(prefix, "_all")
        addnode(allid, "", arm, "All", "All", "all", a, xs[1], 189.5, "All participants")
        for (risk in c("High", "Low")) {
            z <- a[a$risk == risk, , drop = FALSE]
            rid <- paste(prefix, tolower(risk), sep = "_")
            ry <- if (risk == "High") 
                103
            else 276
            addnode(rid, allid, arm, risk, "All", "risk", z, xs[2], ry, screening_label(risk))
            for (j in seq_along(route_order)) {
                route <- route_order[j]
                s <- z[z$route == route, , drop = FALSE]
                nid <- paste(rid, j, sep = "_")
                yy <- if (risk == "High") 
                  c(50, 103, 156)[j]
                else c(223, 276, 329)[j]
                addnode(nid, rid, arm, risk, route, "route", s, xs[3], yy, route)
                addnode(paste0(nid, "_indication"), nid, arm, risk, route, "indication", s[s$indication, , drop = FALSE], xs[4], yy, 
                  "CIED indication")
            }
        }
    }
    nodes <- do.call(rbind, node_rows)
    for (i in seq_len(nrow(nodes))) if (nzchar(nodes$parent_id[i])) {
        child <- nodes[i, ]
        parent <- nodes[match(child$parent_id, nodes$node_id), ]
        edge_rows[[length(edge_rows) + 1L]] <- data.frame(parent_id = parent$node_id, child_id = child$node_id, arm = child$arm, child_n = child$n, 
            parent_n = parent$n, child_percent_of_parent = if (parent$n > 0) 
                100 * child$n/parent$n
            else NA_real_, stringsAsFactors = FALSE)
    }
    edges <- do.call(rbind, edge_rows)
    stopifnot(nrow(nodes) == 30L, nrow(edges) == 28L, sum(nodes$n[nodes$stage == "all"]) == nrow(d), sum(nodes$n[nodes$stage == "route"]) == 
        nrow(d), sum(nodes$n[nodes$stage == "indication"]) == sum(d$indication), sum(nodes$implantations[nodes$stage == "indication"]) == 
        sum(d$cied_implantation == 1), all(nodes$implantations <= nodes$n))
    WIDTH <- 640
    HEIGHT <- 380
    FIGURE_WIDTH_MM <- 180
    FIGURE_HEIGHT_MM <- FIGURE_WIDTH_MM * HEIGHT/WIDTH
    BLUE <- "#3C5488"
    GREY <- "#E5E7EB"
    INK <- "#20252B"
    RADIUS <- 10
    draw_scale <- FIGURE_WIDTH_MM/WIDTH
    put <- function(x, y, s, size = 10, ha = 0.5, bold = FALSE, color = INK) {
        ce <- size * 1.6 * draw_scale/12
        text(x, y, s, cex = ce, adj = c(ha, 0.5), font = if (bold) 
            2
        else 1, col = color)
    }
    curve <- function(x0, y0, x1, y1) {
        t <- seq(0, 1, length.out = 120)
        dx <- x1 - x0
        x <- (1 - t)^3 * x0 + 3 * (1 - t)^2 * t * (x0 + dx * 0.55) + 3 * (1 - t) * t^2 * (x1 - dx * 0.55) + t^3 * x1
        y <- (1 - t)^3 * y0 + 3 * (1 - t)^2 * t * y0 + 3 * (1 - t) * t^2 * y1 + t^3 * y1
        lines(x, y, col = "#8F98A2", lwd = 0.9 * draw_scale)
        arrows(x[115], y[115], x[120], y[120], length = 0.055 * draw_scale, angle = 22, col = "#8F98A2", lwd = 0.9 * draw_scale)
    }
    draw_pie <- function(r) {
        theta <- seq(0, 2 * pi, length.out = 240)
        x <- r$x
        y <- r$y
        if (r$n == 0) {
            polygon(x + RADIUS * cos(theta), y + RADIUS * sin(theta), col = "white", border = "#AEB6C0", lwd = 0.8 * draw_scale, lty = 2)
            put(x, y, "NE", 9, color = "#7A838D")
        }
        else {
            polygon(x + RADIUS * cos(theta), y + RADIUS * sin(theta), col = GREY, border = "white", lwd = 0.5 * draw_scale)
            if (r$implantations > 0) {
                tt <- seq(-pi/2, -pi/2 + 2 * pi * r$implantations/r$n, length.out = max(10, ceiling(240 * r$implantations/r$n)))
                polygon(c(x, x + RADIUS * cos(tt)), c(y, y + RADIUS * sin(tt)), col = BLUE, border = NA)
            }
            lines(x + RADIUS * cos(theta), y + RADIUS * sin(theta), col = "#8F98A2", lwd = 0.55 * draw_scale)
        }
        label <- r$label
        if (r$stage == "risk") {
            label <- if (r$risk == "High") 
                "Flagged by the AI-ECG\nscreening algorithm"
            else "Not flagged by the AI-ECG\nscreening algorithm"
            ce <- 10.5 * 1.6 * draw_scale/12
            lw <- max(strwidth(strsplit(label, "\n", fixed = TRUE)[[1]], cex = ce, font = 2))
            lh <- strheight(label, cex = ce, font = 2)
            ly <- y - RADIUS - 10
            rect(x - lw/2 - 1, ly - lh/2 - 1, x + lw/2 + 1, ly + lh/2 + 1, col = "white", border = NA)
            put(x, ly, label, 10.5, bold = TRUE)
        }
        else {
            if (r$stage == "route" && label == "Trial-arranged ambulatory ECG monitoring")
                label <- "Trial-arranged ambulatory\nECG monitoring"
            put(x, y - RADIUS - if (grepl("\n", label, fixed = TRUE)) 7 else 5,
                label, if (r$stage == "route") 10 else 10.5, bold = TRUE)
        }
        parts <- c(format(r$implantations, big.mark = ",", trim = TRUE), paste0("/", format(r$n, big.mark = ",", trim = TRUE)))
        ww <- strwidth(parts, cex = 16 * draw_scale/12)
        xx <- x - sum(ww)/2
        put(xx, y + 15, parts[1], 10, ha = 0)
        xx <- xx + ww[1]
        put(xx, y + 15, parts[2], 10, ha = 0)
        if (r$n > 0) {
            vals <- c(r$implantation_percent, r$lower_percent, r$upper_percent)
            digits <- if (any(vals > 0 & vals < 0.1)) 
                2
            else 1
            put(x, y + 20, sprintf(paste0("%.", digits, "f%%"), r$implantation_percent), 10, bold = TRUE)
            put(x, y + 25, sprintf(paste0("(%.", digits, "f-%.", digits, "f%%)"), r$lower_percent, r$upper_percent), 9)
        }
        else {
            put(x, y + 20, "NE", 10)
            put(x, y + 25, "(no participants)", 8, color = "#6B7280")
        }
    }
    draw <- function() {
        par(mar = rep(0, 4), xaxs = "i", yaxs = "i", family = "sans", xpd = NA)
        plot.new()
        plot.window(c(0, WIDTH), c(HEIGHT, 0), asp = 1)
        put(160, 7, "Intervention", 17, bold = TRUE)
        put(480, 7, "Control", 17, bold = TRUE)
        segments(320, 31, 320, 358, col = "#D5DAE1", lwd = 0.8 * draw_scale, lty = 3)
        for (i in seq_len(nrow(edges))) {
            e <- edges[i, ]
            p <- nodes[match(e$parent_id, nodes$node_id), ]
            c <- nodes[match(e$child_id, nodes$node_id), ]
            sgn <- sign(c$x - p$x)
            x0 <- p$x + sgn * (RADIUS + 1)
            x1 <- c$x - sgn * (RADIUS + 1)
            curve(x0, p$y, x1, c$y)
            xx <- x0 + (x1 - x0) * 0.63
            yy <- p$y + (c$y - p$y) * 0.69 - 4
            lab <- if (is.finite(e$child_percent_of_parent)) 
                sprintf(if (e$child_percent_of_parent > 0 && e$child_percent_of_parent < 0.1) 
                  "%.2f%%"
                else "%.1f%%", e$child_percent_of_parent)
            else "NE"
            put(xx, yy, lab, 10)
        }
        for (i in seq_len(nrow(nodes))) draw_pie(nodes[i, ])
        rect(230, 366, 235, 371, col = BLUE, border = NA)
        put(238, 368.5, "CIED implantation", 10, ha = 0)
        rect(350, 366, 355, 371, col = GREY, border = "#8F98A2", lwd = 0.4 * draw_scale)
        put(358, 368.5, "No implantation", 10, ha = 0)
    }
    export_figure <- function(ext) {
        target <- file.path(out, paste0("Supplemental_Fig_07.", ext))
        tmp <- paste0(target, ".", Sys.getpid(), ".pending")
        if (ext == "pdf") 
            portable_pdf(tmp, width = FIGURE_WIDTH_MM/25.4, height = FIGURE_HEIGHT_MM/25.4, pointsize = 12, family = "sans", bg = "white")
        device <- dev.cur()
        tryCatch(draw(), finally = {
            if (device %in% dev.list()) 
                dev.off(device)
        })
        stopifnot(file.exists(tmp), file.info(tmp)$size > 1000)
        if (ext == "pdf") {
            bytes <- readBin(tmp, "raw", n = file.info(tmp)$size)
            stopifnot(rawToChar(bytes[1:4]) == "%PDF", grepl("%%EOF", rawToChar(tail(bytes, 40)), fixed = TRUE))
        }
        stopifnot(file.rename(tmp, target))
    }
    for (ext in "pdf") export_figure(ext)
}

precise_core_figures <- function(root) {
    portable_pdf <- function(...) {
        if (capabilities("cairo")) 
            grDevices::cairo_pdf(...)
        else grDevices::pdf(...)
    }
    suppressPackageStartupMessages(library(jsonlite))
    out <- file.path(root, "figures")
    dir.create(out, showWarnings = FALSE, recursive = TRUE)
    read_csv <- function(f) read.csv(file.path(root, f), check.names = FALSE, stringsAsFactors = FALSE)
    E <- read_csv("analysis/effects.csv")
    I <- read_csv("analysis/interactions.csv")
    Y <- read_csv("analysis/pp_intervention_yields.csv")
    YG <- read_csv("analysis/pp_intervention_global_tests.csv")
    D <- read_csv("data/PRECISE-AI.csv")
    stopifnot(all(D$ai_ecg_screening_flag %in% 0:1), identical(D$screening_group, ifelse(D$ai_ecg_screening_flag == 1, "High risk", 
        "Low risk")))
    precise_core_pathways(root)
    BLUE <- "#4B638F"
    GRAY <- "#F1F2F4"
    INK <- "#303135"
    RULE <- "#7B7D80"
    PITCH <- 3.8
    GAP <- 0.8
    strata <- c(I = "Potentially suboptimal rhythm evaluation (Stratum I)", II = "Risk indicators supporting rhythm monitoring (Stratum II)", 
        III = "No prespecified rhythm risk indicators (Stratum III)")
    holter <- c(`I-a` = "Stratum I (without Holter within 1 year)", `I-b` = "Stratum I (with Holter within 1 year)", `II-a` = "Stratum II (without Holter within 1 year)", 
        `II-b` = "Stratum II (with Holter within 1 year)", `III-a` = "Stratum III (without Holter within 1 year)", `III-b` = "Stratum III (with Holter within 1 year)")
    blocks <- list(clinical_stratum = list(title = "Clinical strata", levels = strata), clinical_holter_stratum = list(title = "Strata by recent Holter monitoring", 
        levels = holter), age_group = list(title = "Age", levels = c(`65-79 years` = "65-79 years", `80-90 years` = "80-90 years")), 
        sex_group = list(title = "Sex", levels = c(Female = "Female", Male = "Male")))
    num <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
    pstr <- function(p, method = "") {
        if (length(p) == 0L || is.na(p) || !is.finite(p)) 
            return("NE")
        paste0(if (p < 0.001) 
            "<0.001"
        else sprintf("%.3f", p), if (grepl("Firth", method, ignore.case = TRUE)) 
            "†"
        else "")
    }
    arm <- function(r, k) {
        n <- r[[paste0("n", k)]]
        e <- r[[paste0("event", k)]]
        if (n == 0) 
            "0/0 (NE)"
        else sprintf("%s/%s (%.1f%%)", num(e), num(n), 100 * e/n)
    }
    pick <- function(...) {
        q <- E
        v <- list(...)
        for (k in names(v)) q <- q[q[[k]] == v[[k]], , drop = FALSE]
        stopifnot(nrow(q) == 1)
        q
    }
    interaction <- function(pop, outcome, var) {
        q <- I[I$population == pop & I$outcome == outcome & I$subgroup_var == var, , drop = FALSE]
        stopifnot(nrow(q) == 1)
        q
    }
    canvas <- function(height) {
        par(fig = c(0, 1, 0, 1), new = FALSE)
        par(mar = rep(0, 4), xaxs = "i", yaxs = "i", family = "sans", xpd = NA)
        plot.new()
        plot.window(xlim = c(0, 180), ylim = c(height, 0))
    }
    put <- function(x, y, s, size = 5.8, bold = FALSE, ha = "left", color = INK) {
        if (length(s) == 0L || (is.character(s) && is.na(s))) 
            s <- ""
        a <- switch(ha, left = 0, center = 0.5, right = 1)
        text(x, y, labels = s, cex = size/12, font = if (bold) 
            2
        else 1, adj = c(a, 0.5), col = color)
    }
    line <- function(x1, y1, x2, y2, col = RULE, lwd = 0.65, lty = 1) segments(x1, y1, x2, y2, col = col, lwd = lwd, lty = lty)
    band <- function(y, label, x = 2, w = 176.5, size = 6) {
        rect(x, y - PITCH * 0.43, x + w, y + PITCH * 0.43, col = GRAY, border = NA)
        put(x + 1.4, y, label, size, bold = TRUE)
    }
    arrow <- function(x1, y1, x2, y2, col = BLUE, lwd = 0.7) arrows(x1, y1, x2, y2, length = 0.035, angle = 25, code = 2, col = col, 
        lwd = lwd)
    X <- list(label = 3.4, inter = 60, control = 86, left = 102, right = 129, or = 144.8, p = 159.5, interaction = 171.7)
    tx <- function(v, lo, hi, x = X) x$left + (log(v) - log(lo))/(log(hi) - log(lo)) * (x$right - x$left)
    header <- function(top, label = "Outcome / population", x = X) {
        put(x$label, top + 2, label, 5.9, TRUE)
        put(x$inter, top + 2, "Intervention", 5.7, TRUE, "center")
        put(x$control, top + 2, "Control", 5.7, TRUE, "center")
        put(x$or, top + 2, "Odds ratio (95% CI)", 5.6, TRUE, "center")
        put(x$p, top + 2, expression(bolditalic(P) ~ bold("value")), 5.3, TRUE, "center")
        put(x$interaction, top + 2, expression(bolditalic(P) ~ bold("for interaction")), 5, TRUE, "center")
        line(2, top + 5.3, 178.5, top + 5.3, lwd = 0.75)
    }
    forest_axis <- function(top, bottom, lo = 0.25, hi = 32, x = X) {
        line(tx(1, lo, hi, x), top, tx(1, lo, hi, x), bottom - 3, col = "#AAADB3", lwd = 0.85, lty = 3)
        line(x$left, bottom, x$right, bottom, lwd = 0.75)
        for (t in c(0.25, 1, 4, 32)) {
            xx <- tx(t, lo, hi, x)
            line(xx, bottom, xx, bottom + 1.5)
            put(xx, bottom + 3.3, as.character(t), 5.2, ha = "center")
        }
        put((x$left + x$right)/2, bottom + 7.2, "Intervention / control odds ratio", 5.3, ha = "center")
    }
    estimate <- function(y, r, lo = 0.25, hi = 32, x = X, counts = TRUE, open = FALSE) {
        a <- if (counts) 
            arm(r, "1")
        else ""
        b <- if (counts) 
            arm(r, "0")
        else ""
        put(x$inter, y, a, ha = "center")
        put(x$control, y, b, ha = "center")
        if (is.finite(r$OR)) {
            line(tx(max(r$lower, lo), lo, hi, x), y, tx(min(r$upper, hi), lo, hi, x), y, col = BLUE, lwd = 1.2, lty = if (open) 
                2
            else 1)
            if (r$lower < lo) 
                arrow(x$left + 2.4, y, x$left - 0.6, y, lwd = 1)
            if (r$upper > hi) 
                arrow(x$right - 2.4, y, x$right + 0.6, y, lwd = 1)
            if (r$OR < lo || r$OR > hi) 
                arrow(tx(min(max(r$OR, lo), hi), lo, hi, x), y, if (r$OR < lo) 
                  x$left - 0.8
                else x$right + 0.8, y)
            else points(tx(r$OR, lo, hi, x), y, pch = 22, cex = 0.48, bg = if (open) 
                "white"
            else BLUE, col = BLUE, lwd = 0.65)
            fmt <- function(v) if (v > 0 && v < 0.01) 
                sprintf("%.3f", v)
            else sprintf("%.2f", v)
            ci <- sprintf("%s (%s-%s)", fmt(r$OR), fmt(r$lower), fmt(r$upper))
        }
        else ci <- "NE"
        pp <- pstr(r$p, r$method)
        put(x$or, y, ci, 5.6, ha = "center")
        put(x$p, y, pp, 5.6, ha = "center")
    }
    draw_subgroup <- function(pop, outcome) {
        canvas(92)
        x <- X
        x$inter <- 74
        x$control <- 96
        x$left <- 108
        header(3, "Subgroup", x)
        y <- 11.9
        forest_axis(8.4, 80.6, x = x)
        for (v in names(blocks)) {
            b <- blocks[[v]]
            band(y, b$title)
            it <- interaction(pop, outcome, v)
            put(x$interaction, y, pstr(it$p, it$method), 5.6, ha = "center")
            for (lev in names(b$levels)) {
                y <- y + PITCH
                r <- pick(population = pop, adjustment = "Crude", outcome = outcome, subgroup_var = v, subgroup_level = lev)
                put(x$label + 1.6, y, b$levels[[lev]], 5.8)
                estimate(y, r, x = x)
            }
            y <- y + PITCH + GAP
        }
    }
    draw_yield <- function(fig) {
        canvas(100)
        q <- Y[Y$figure == fig, , drop = FALSE]
        stopifnot(nrow(q) == 14)
        pal <- c("#3C5488", "#1F8A8A", "#646A9D", "#8491B4", "#5E7773", "#507D93")
        xl <- 68
        xr <- 111
        put(3, 4, "Subgroup", 6, TRUE)
        put((xl + xr)/2, 4, "Proportion (%)", 6, TRUE, "center")
        put(140, 4, "n/N; % (95% CI)", 6, TRUE, "center")
        put(172, 4, expression(bold("Global") ~ bolditalic(P) ~ bold("value")), 5.7, TRUE, "center")
        line(2, 7.6, 178, 7.6, col = INK, lwd = 0.55)
        for (t in seq(0, 100, 25)) {
            x <- xl + (xr - xl) * t/100
            line(x, 8.2, x, 83.7, col = "#DDDDDD", lwd = 0.4)
            line(x, 85.3, x, 86.1, col = INK, lwd = 0.45)
            put(x, 89, as.character(t), 5.6, ha = "center")
        }
        line(xl, 85.3, xr, 85.3, col = INK, lwd = 0.45)
        put((xl + xr)/2, 94, "Proportion (%)", 6, ha = "center")
        ys <- list(Overall = 11.8, clinical_stratum = c(21, 24.8, 28.6), clinical_holter_stratum = c(37, 40.8, 44.6, 48.4, 52.2, 56), 
            age_group = c(64.4, 68.2), sex_group = c(76.6, 80.4))
        heads <- c(clinical_stratum = 17.2, clinical_holter_stratum = 33.2, age_group = 60.6, sex_group = 72.8)
        for (v in names(ys)) {
            z <- q[q$subgroup_var == v, , drop = FALSE]
            stopifnot(nrow(z) == length(ys[[v]]))
            if (v != "Overall") {
                band(heads[v], c(clinical_stratum = "Prespecified clinical strata", clinical_holter_stratum = "Strata by recent Holter monitoring", 
                  age_group = "Age group", sex_group = "Recorded sex")[v], w = 176)
                gp <- YG[YG$figure == fig & YG$subgroup_var == v, , drop = FALSE]
                stopifnot(nrow(gp) == 1)
                gp$p_display <- pstr(gp$p_value)
                put(172, heads[v], gp$p_display, ha = "center")
            }
            for (j in seq_len(nrow(z))) {
                r <- z[j, , drop = FALSE]
                y <- ys[[v]][j]
                cl <- if (v == "Overall") 
                  INK
                else pal[j]
                label <- if (v == "clinical_stratum") 
                  strata[r$level]
                else if (v == "clinical_holter_stratum") 
                  holter[r$level]
                else r$label
                put(3, y, label, bold = v == "Overall")
                stopifnot(strwidth(label, cex = 5.8/12, font = if (v == "Overall") 
                  2
                else 1) < xl - 5)
                if (r$total > 0) {
                  lo <- xl + (xr - xl) * r$lower
                  hi <- xl + (xr - xl) * r$upper
                  line(lo, y, hi, y, cl, 1.1)
                  line(lo, y - 0.55, lo, y + 0.55, cl, 0.7)
                  line(hi, y - 0.55, hi, y + 0.55, cl, 0.7)
                  points(xl + (xr - xl) * r$proportion, y, pch = 21, cex = 0.45, bg = cl, col = cl, lwd = 0.3)
                  st <- sprintf("%.1f (%.1f-%.1f)", r$percent, r$lower_percent, r$upper_percent)
                  ss <- c(as.character(r$events), paste0("/", r$total, "; "), st)
                  ww <- strwidth(ss, cex = 5.8/12)
                  xx <- 140 - sum(ww)/2
                  for (k in 1:3) {
                    put(xx, y, ss[k], color = cl)
                    xx <- xx + ww[k]
                  }
                }
                else {
                  put(140, y, "0/0; NE", ha = "center", color = cl)
                }
            }
        }
    }
    figs <- list(
        Supplemental_Fig_03 = list(h = 92, f = function() draw_subgroup("Overall", "cied_implantation")),
        Supplemental_Fig_04 = list(h = 92, f = function() draw_subgroup("High risk", "cied_implantation")),
        Supplemental_Fig_05 = list(h = 92, f = function() draw_subgroup("Overall", "sick_sinus_syndrome")),
        Supplemental_Fig_06 = list(h = 92, f = function() draw_subgroup("High risk", "sick_sinus_syndrome")),
        Supplemental_Fig_08 = list(h = 100, f = function() draw_yield("indication")),
        Supplemental_Fig_09 = list(h = 100, f = function() draw_yield("implantation")))
    for (stem in names(figs)) {
        z <- figs[[stem]]
        for (ext in "pdf") {
            target <- file.path(out, paste0(stem, ".", ext))
            tmp <- paste0(target, ".", Sys.getpid(), ".pending")
            if (ext == "pdf") 
                portable_pdf(tmp, width = 180/25.4, height = z$h/25.4, family = "sans", pointsize = 12)
            dv <- dev.cur()
            tryCatch(z$f(), finally = {
                if (dv %in% dev.list()) 
                  dev.off(dv)
            })
            stopifnot(file.info(tmp)$size > 1000, file.rename(tmp, target))
        }
    }
}

precise_core_risk_window <- function(root) {
    portable_pdf <- function(...) {
        if (capabilities("cairo")) 
            grDevices::cairo_pdf(...)
        else grDevices::pdf(...)
    }
    precise_publication_labels()
    precise_core_helpers()
    d <- prepare_analysis_factors(read.csv(file.path(root, "data", "PRECISE-AI.csv"), check.names = FALSE, stringsAsFactors = FALSE, 
        na.strings = ""))
    outcomes <- c("cied_implantation", precise_secondary_outcomes)
    stopifnot(nrow(d) > 0L, !anyNA(d[c("group", "ai_ecg_screening_flag", "high_risk_ecg_count_1y", outcomes)]), all(d$ai_ecg_screening_flag %in% 
        0:1), all(d$high_risk_ecg_count_1y >= 0), all(d$high_risk_ecg_count_1y == floor(d$high_risk_ecg_count_1y)), all(d$ecg_count_1y >= 
        1))
    levels2 <- c(PROTOCOL_HIGH_LABEL, PROTOCOL_LOW_LABEL)
    d$protocol_risk_group <- factor(ifelse(d$ai_ecg_screening_flag == 1 | d$high_risk_ecg_count_1y > 0, levels2[1], levels2[2]), levels = levels2)
    stopifnot(!anyNA(d$protocol_risk_group), all(d$protocol_risk_group[d$ai_ecg_screening_flag == 1] == levels2[1]), sum(d$protocol_risk_group == 
        levels2[1]) == 1021L, sum(d$protocol_risk_group == levels2[2]) == 10471L, all(d$high_risk_ecg_count_1y[d$ai_ecg_screening_flag == 
        1] > 0))
    counts <- do.call(rbind, lapply(levels2, function(s) do.call(rbind, lapply(c("Intervention", "Control"), function(a) {
        z <- d[d$protocol_risk_group == s & d$group == a, , drop = FALSE]
        data.frame(stratum = s, arm = a, n = nrow(z), stringsAsFactors = FALSE)
    }))))
    stopifnot(all(counts$n > 0), sum(counts$n) == nrow(d))
    effects <- interactions <- list()
    for (v in outcomes) {
        rule <- block_rule(d, v, "protocol_risk_group")
        it <- interaction(d, v, "protocol_risk_group", "Randomized population", use_firth = rule$firth, block_id = paste0("protocol_risk:", 
            v), trigger_levels = rule$trigger_levels)
        it$outcome_label <- outlabs[[v]]
        interactions[[v]] <- it
        for (s in levels2) {
            z <- d[d$protocol_risk_group == s, , drop = FALSE]
            e <- effect(z, v, s, sg = "protocol_risk_group", lv = s, use_firth = rule$firth, block_id = paste0("protocol_risk:", v), 
                trigger_levels = rule$trigger_levels)
            e$interaction_p <- it$p
            e$interaction_method <- it$method
            e$interaction_status <- it$status
            effects[[length(effects) + 1L]] <- e
        }
    }
    E <- do.call(rbind, effects)
    I <- do.call(rbind, interactions)
    stopifnot(nrow(E) == 10L, nrow(I) == 5L, all(E$status %in% c("OK", "No events in either arm")), all(I$status %in% c("OK", "No events in either arm within a subgroup", 
        "Ordinary interaction NE: nonfinite MLE in treatment-by-subgroup outcome cells")))
    for (v in outcomes) {
        q <- E[E$outcome == v, , drop = FALSE]
        stopifnot(sum(q$event1) == sum(d[[v]][d$group == "Intervention"]), sum(q$event0) == sum(d[[v]][d$group == "Control"]), sum(q$n1) == 
            sum(d$group == "Intervention"), sum(q$n0) == sum(d$group == "Control"))
    }
    fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
    fmt_arm <- function(e, n) sprintf("%s/%s (%.1f%%)", fmt_n(e), fmt_n(n), 100 * e/n)
    fmt_est <- function(x) {
        if (!is.finite(x)) 
            return("NE")
        if (x > 0 && x < 0.001) 
            return("<0.001")
        sprintf(if (x > 0 && x < 0.01) 
            "%.3f"
        else "%.2f", x)
    }
    pdisplay <- function(p, method = "") {
        if (is.na(p) || !is.finite(p)) 
            return("NE")
        paste0(if (p < 0.001) 
            "<0.001"
        else sprintf("%.3f", p), if (grepl("Firth", method, ignore.case = TRUE)) 
            "†"
        else "")
    }
    E$intervention_display <- mapply(fmt_arm, E$event1, E$n1)
    E$control_display <- mapply(fmt_arm, E$event0, E$n0)
    E$OR_CI <- vapply(seq_len(nrow(E)), function(j) if (!is.finite(E$OR[j])) 
        "NE"
    else paste0(fmt_est(E$OR[j]), " (", fmt_est(E$lower[j]), "-", fmt_est(E$upper[j]), ")"), character(1))
    E$p_display <- mapply(pdisplay, E$p, E$method)
    E$interaction_p_display <- mapply(pdisplay, E$interaction_p, E$interaction_method)
    I$p_display <- mapply(pdisplay, I$p, I$method)
    BLUE <- "#4B638F"
    GRAY <- "#F1F2F4"
    INK <- "#303135"
    RULE <- "#7B7D80"
    H <- 83
    PITCH <- 3.8
    GAP <- 0.8
    X <- list(label = 3.4, inter = 60, control = 86, left = 102, right = 129, or = 144.8, p = 159.5, interaction = 171.7)
    put <- function(x, y, s, size = 5.8, bold = FALSE, ha = "left", color = INK) text(x, y, labels = s, cex = size/12, font = if (bold) 
        2
    else 1, adj = c(switch(ha, left = 0, center = 0.5, right = 1), 0.5), col = color)
    line <- function(x1, y1, x2, y2, col = RULE, lwd = 0.65, lty = 1) segments(x1, y1, x2, y2, col = col, lwd = lwd, lty = lty)
    tx <- function(v) X$left + (log(v) - log(0.25))/(log(32) - log(0.25)) * (X$right - X$left)
    arrow <- function(x1, y1, x2, y2) arrows(x1, y1, x2, y2, length = 0.035, angle = 25, code = 2, col = BLUE, lwd = 1)
    display_labels <- setNames(c(PROTOCOL_HIGH_LABEL, PROTOCOL_LOW_LABEL), levels2)
    draw <- function() {
        par(mar = rep(0, 4), xaxs = "i", yaxs = "i", family = "sans", xpd = NA)
        plot.new()
        plot.window(c(0, 180), c(H, 0))
        put(X$label, 4.4, "Outcome / protocol defined risk", 5.9, TRUE)
        put(X$inter, 4.4, "Intervention", 5.7, TRUE, "center")
        put(X$control, 4.4, "Control", 5.7, TRUE, "center")
        put(X$or, 4.4, "Odds ratio (95% CI)", 5.6, TRUE, "center")
        put(X$p, 4.4, expression(bolditalic(P) ~ bold("value")), 5.3, TRUE, "center")
        put(X$interaction, 4.4, expression(bolditalic(P) ~ bold("for interaction")), 5, TRUE, "center")
        line(2, 7.7, 178.5, 7.7, lwd = 0.75)
        bottom <- 71.9
        line(tx(1), 7.8, tx(1), bottom - 3, col = "#AAADB3", lwd = 0.85, lty = 3)
        line(X$left, bottom, X$right, bottom, lwd = 0.75)
        for (v in c(0.25, 1, 4, 32)) {
            xx <- tx(v)
            line(xx, bottom, xx, bottom + 1.5)
            put(xx, bottom + 3.3, as.character(v), 5.2, ha = "center")
        }
        put((X$left + X$right)/2, bottom + 7.2, "Intervention / control odds ratio", 5.3, ha = "center")
        y <- 10.8
        for (v in outcomes) {
            prefix <- if (v == "cied_implantation") 
                "Primary outcome: "
            else "Secondary outcome: "
            rect(2, y - PITCH * 0.43, 178.5, y + PITCH * 0.43, col = GRAY, border = NA)
            put(3.4, y, paste0(prefix, outlabs[[v]]), 6, TRUE)
            it <- I[I$outcome == v, , drop = FALSE]
            put(X$interaction, y, it$p_display, 5.6, ha = "center")
            for (s in levels2) {
                y <- y + PITCH
                r <- E[E$outcome == v & E$subgroup_level == s, , drop = FALSE]
                stopifnot(nrow(r) == 1)
                lab <- display_labels[[s]]
                put(5, y, lab, 5.6)
                stopifnot(5 + strwidth(lab, cex = 5.6/12) < X$inter - strwidth(r$intervention_display, cex = 5.8/12)/2 - 1)
                put(X$inter, y, r$intervention_display, ha = "center")
                put(X$control, y, r$control_display, ha = "center")
                if (is.finite(r$OR)) {
                  line(tx(max(r$lower, 0.25)), y, tx(min(r$upper, 32)), y, col = BLUE, lwd = 1.2)
                  if (r$lower < 0.25) 
                    arrow(X$left + 2.4, y, X$left - 0.6, y)
                  if (r$upper > 32) 
                    arrow(X$right - 2.4, y, X$right + 0.6, y)
                  if (r$OR < 0.25 || r$OR > 32) 
                    arrow(tx(min(max(r$OR, 0.25), 32)), y, if (r$OR < 0.25) 
                      X$left - 0.8
                    else X$right + 0.8, y)
                  else points(tx(r$OR), y, pch = 22, cex = 0.48, bg = BLUE, col = BLUE, lwd = 0.65)
                }
                put(X$or, y, r$OR_CI, 5.6, ha = "center")
                put(X$p, y, r$p_display, 5.6, ha = "center")
            }
            y <- y + PITCH + GAP
        }
    }
    for (ext in "pdf") {
        target <- file.path(root, "figures", paste0("Supplemental_Fig_10.", ext))
        tmp <- paste0(target, ".pending")
        if (ext == "pdf") 
            portable_pdf(tmp, width = 180/25.4, height = H/25.4, family = "sans", pointsize = 12)
        draw()
        dev.off()
        stopifnot(file.rename(tmp, target))
    }
}

precise_core_risk_history <- function(root) {
    portable_pdf <- function(...) {
        if (capabilities("cairo")) 
            grDevices::cairo_pdf(...)
        else grDevices::pdf(...)
    }
    precise_publication_labels()
    precise_core_helpers()
    d <- prepare_analysis_factors(read.csv(file.path(root, "data", "PRECISE-AI.csv"), check.names = FALSE, stringsAsFactors = FALSE, 
        na.strings = ""))
    outcomes <- c("cied_implantation", precise_secondary_outcomes)
    stopifnot(nrow(d) > 0L, !anyNA(d[c("group", "ai_ecg_screening_flag", "high_risk_ecg_count_1y", outcomes)]), all(d$ai_ecg_screening_flag %in% 
        0:1), all(d$high_risk_ecg_count_1y >= 0), all(d$high_risk_ecg_count_1y == floor(d$high_risk_ecg_count_1y)), all(d$ecg_count_1y >= 
        1))
    levels3 <- THREE_GROUP_LABELS
    d$screening_history_group <- factor(ifelse(d$ai_ecg_screening_flag == 1, levels3[1], ifelse(d$high_risk_ecg_count_1y > 0, levels3[2], 
        levels3[3])), levels = levels3)
    stopifnot(!anyNA(d$screening_history_group), sum(d$screening_history_group == levels3[1]) == sum(d$ai_ecg_screening_flag == 1))
    counts <- do.call(rbind, lapply(levels3, function(s) do.call(rbind, lapply(c("Intervention", "Control"), function(a) {
        z <- d[d$screening_history_group == s & d$group == a, , drop = FALSE]
        data.frame(stratum = s, arm = a, n = nrow(z), stringsAsFactors = FALSE)
    }))))
    stopifnot(all(counts$n > 0), sum(counts$n) == nrow(d))
    effects <- interactions <- list()
    for (v in outcomes) {
        rule <- block_rule(d, v, "screening_history_group")
        it <- interaction(d, v, "screening_history_group", "Randomized population", use_firth = rule$firth, block_id = paste0("screening_history:", 
            v), trigger_levels = rule$trigger_levels)
        it$outcome_label <- outlabs[[v]]
        interactions[[v]] <- it
        for (s in levels3) {
            z <- d[d$screening_history_group == s, , drop = FALSE]
            e <- effect(z, v, s, sg = "screening_history_group", lv = s, use_firth = rule$firth, block_id = paste0("screening_history:", 
                v), trigger_levels = rule$trigger_levels)
            e$interaction_p <- it$p
            e$interaction_method <- it$method
            e$interaction_status <- it$status
            effects[[length(effects) + 1L]] <- e
        }
    }
    E <- do.call(rbind, effects)
    I <- do.call(rbind, interactions)
    stopifnot(nrow(E) == 15L, nrow(I) == 5L, all(E$status %in% c("OK", "No events in either arm")), all(I$status %in% c("OK", "No events in either arm within a subgroup", 
        "Ordinary interaction NE: nonfinite MLE in treatment-by-subgroup outcome cells")))
    for (v in outcomes) {
        q <- E[E$outcome == v, , drop = FALSE]
        stopifnot(sum(q$event1) == sum(d[[v]][d$group == "Intervention"]), sum(q$event0) == sum(d[[v]][d$group == "Control"]), sum(q$n1) == 
            sum(d$group == "Intervention"), sum(q$n0) == sum(d$group == "Control"))
    }
    fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
    fmt_arm <- function(e, n) sprintf("%s/%s (%.1f%%)", fmt_n(e), fmt_n(n), 100 * e/n)
    fmt_est <- function(x) {
        if (!is.finite(x)) 
            return("NE")
        if (x > 0 && x < 0.001) 
            return("<0.001")
        sprintf(if (x > 0 && x < 0.01) 
            "%.3f"
        else "%.2f", x)
    }
    pdisplay <- function(p, method = "") {
        if (is.na(p) || !is.finite(p)) 
            return("NE")
        paste0(if (p < 0.001) 
            "<0.001"
        else sprintf("%.3f", p), if (grepl("Firth", method, ignore.case = TRUE)) 
            "†"
        else "")
    }
    E$intervention_display <- mapply(fmt_arm, E$event1, E$n1)
    E$control_display <- mapply(fmt_arm, E$event0, E$n0)
    E$OR_CI <- vapply(seq_len(nrow(E)), function(j) if (!is.finite(E$OR[j])) 
        "NE"
    else paste0(fmt_est(E$OR[j]), " (", fmt_est(E$lower[j]), "-", fmt_est(E$upper[j]), ")"), character(1))
    E$p_display <- mapply(pdisplay, E$p, E$method)
    E$interaction_p_display <- mapply(pdisplay, E$interaction_p, E$interaction_method)
    I$p_display <- mapply(pdisplay, I$p, I$method)
    BLUE <- "#4B638F"
    GRAY <- "#F1F2F4"
    INK <- "#303135"
    RULE <- "#7B7D80"
    PITCH <- 3.8
    GAP <- 0.8
    X <- list(label = 3.4, inter = 60, control = 86, left = 102, right = 129, or = 144.8, p = 159.5, interaction = 171.7)
    put <- function(x, y, s, size = 5.8, bold = FALSE, ha = "left", color = INK) text(x, y, labels = s, cex = size/12, font = if (bold) 
        2
    else 1, adj = c(switch(ha, left = 0, center = 0.5, right = 1), 0.5), col = color)
    line <- function(x1, y1, x2, y2, col = RULE, lwd = 0.65, lty = 1) segments(x1, y1, x2, y2, col = col, lwd = lwd, lty = lty)
    tx <- function(v) X$left + (log(v) - log(0.25))/(log(32) - log(0.25)) * (X$right - X$left)
    arrow <- function(x1, y1, x2, y2) arrows(x1, y1, x2, y2, length = 0.035, angle = 25, code = 2, col = BLUE, lwd = 1)
    display_labels <- setNames(THREE_GROUP_LABELS, levels3)
    pdf(NULL, width = 180/25.4, height = 100/25.4, pointsize = 12, family = "Helvetica")
    wrapped_labels <- vapply(display_labels, wrap_publication_label, character(1), width_mm = 44, size = 5.6)
    invisible(dev.off())
    row_height <- setNames(pmax(PITCH, (lengths(strsplit(wrapped_labels, "\n", fixed = TRUE))) * 2.4 + 1), names(display_labels))
    block_height <- PITCH + sum(row_height) + GAP
    H <- 10.8 + length(outcomes) * block_height + 9
    draw <- function() {
        par(mar = rep(0, 4), xaxs = "i", yaxs = "i", family = "sans", xpd = NA)
        plot.new()
        plot.window(c(0, 180), c(H, 0))
        put(X$label, 4.4, "Outcome / screening-status group", 5.9, TRUE)
        put(X$inter, 4.4, "Intervention", 5.7, TRUE, "center")
        put(X$control, 4.4, "Control", 5.7, TRUE, "center")
        put(X$or, 4.4, "Odds ratio (95% CI)", 5.6, TRUE, "center")
        put(X$p, 4.4, expression(bolditalic(P) ~ bold("value")), 5.3, TRUE, "center")
        put(X$interaction, 4.4, expression(bolditalic(P) ~ bold("for interaction")), 5, TRUE, "center")
        line(2, 7.7, 178.5, 7.7, lwd = 0.75)
        bottom <- H - 10.8
        line(tx(1), 7.8, tx(1), bottom - 3, col = "#AAADB3", lwd = 0.85, lty = 3)
        line(X$left, bottom, X$right, bottom, lwd = 0.75)
        for (v in c(0.25, 1, 4, 32)) {
            xx <- tx(v)
            line(xx, bottom, xx, bottom + 1.5)
            put(xx, bottom + 3.3, as.character(v), 5.2, ha = "center")
        }
        put((X$left + X$right)/2, bottom + 7.2, "Intervention / control odds ratio", 5.3, ha = "center")
        y <- 10.8
        for (v in outcomes) {
            prefix <- if (v == "cied_implantation") 
                "Primary outcome: "
            else "Secondary outcome: "
            rect(2, y - PITCH * 0.43, 178.5, y + PITCH * 0.43, col = GRAY, border = NA)
            put(3.4, y, paste0(prefix, outlabs[[v]]), 6, TRUE)
            it <- I[I$outcome == v, , drop = FALSE]
            put(X$interaction, y, it$p_display, 5.6, ha = "center")
            last_height <- PITCH
            for (s in levels3) {
                y <- y + (last_height + row_height[[s]])/2
                last_height <- row_height[[s]]
                r <- E[E$outcome == v & E$subgroup_level == s, , drop = FALSE]
                stopifnot(nrow(r) == 1)
                lab <- display_labels[[s]]
                put(5, y, wrapped_labels[[s]], 5.6)
                stopifnot(5 + max(strwidth(strsplit(wrapped_labels[[s]], "\n", fixed = TRUE)[[1]], cex = 5.6/12)) < X$inter - strwidth(r$intervention_display, 
                  cex = 5.8/12)/2 - 1)
                put(X$inter, y, r$intervention_display, ha = "center")
                put(X$control, y, r$control_display, ha = "center")
                if (is.finite(r$OR)) {
                  line(tx(max(r$lower, 0.25)), y, tx(min(r$upper, 32)), y, col = BLUE, lwd = 1.2)
                  if (r$lower < 0.25) 
                    arrow(X$left + 2.4, y, X$left - 0.6, y)
                  if (r$upper > 32) 
                    arrow(X$right - 2.4, y, X$right + 0.6, y)
                  if (r$OR < 0.25 || r$OR > 32) 
                    arrow(tx(min(max(r$OR, 0.25), 32)), y, if (r$OR < 0.25) 
                      X$left - 0.8
                    else X$right + 0.8, y)
                  else points(tx(r$OR), y, pch = 22, cex = 0.48, bg = BLUE, col = BLUE, lwd = 0.65)
                }
                put(X$or, y, r$OR_CI, 5.6, ha = "center")
                put(X$p, y, r$p_display, 5.6, ha = "center")
            }
            y <- y + (last_height + PITCH)/2 + GAP
        }
    }
    for (ext in "pdf") {
        target <- file.path(root, "figures", paste0("Supplemental_Fig_11.", ext))
        tmp <- paste0(target, ".pending")
        if (ext == "pdf") 
            portable_pdf(tmp, width = 180/25.4, height = H/25.4, family = "sans", pointsize = 12)
        draw()
        dev.off()
        stopifnot(file.rename(tmp, target))
    }
}

precise_core_tables <- function(root) {
    suppressPackageStartupMessages(library(jsonlite))
    out <- file.path(root, "tables")
    dir.create(out, showWarnings = FALSE, recursive = TRUE)
    d <- read.csv(file.path(root, "data", "PRECISE-AI.csv"), check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
    stopifnot(!anyDuplicated(d$row_number), all(d$group %in% c("Control", "Intervention")))
    stopifnot(all(d$ai_ecg_screening_flag %in% 0:1), identical(d$screening_group, ifelse(d$ai_ecg_screening_flag == 1, "High risk", 
        "Low risk")))
    d$group <- factor(d$group, levels = c("Control", "Intervention"))
    w <- read.csv(file.path(root, "analysis", "propensity_weights.csv"), check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", 
        "NA"))
    stopifnot(!anyDuplicated(w$row_number))
    for (v in intersect(c("PS", "IPTW", "FIRTH_WEIGHT"), names(w))) d[[v]] <- w[[v]][match(d$row_number, w$row_number)]
    if (!"IPTW" %in% names(d)) 
        stop("Required column IPTW absent from propensity_weights.csv")
    map <- data.frame(row_index = integer(), variable = character(), level = character(), stringsAsFactors = FALSE)
    addmap <- function(rows, idx, levels = "1") {
        n <- length(rows)
        map <<- rbind(map, data.frame(row_index = rows, variable = precise_baseline_variables[rep(idx, length.out = n)], level = rep(levels, 
            length.out = n), stringsAsFactors = FALSE))
    }
    addmap(c(4, 5), 1, c("1", "0"))
    addmap(6, 2)
    addmap(8, 3)
    addmap(9, 4, "Continuous")
    addmap(c(11, 12), 5, c("0", "1"))
    addmap(14:16, 6, c("I", "II", "III"))
    addmap(c(19, 20), 7, c("Recent-onset", "Long-standing"))
    addmap(21:33, 8:20)
    addmap(35:36, 21:22)
    addmap(38:39, 23, c("Recent-onset", "Long-standing"))
    addmap(41:42, 24, c("Recent-onset", "Long-standing"))
    addmap(43, 25)
    addmap(45:46, 26, c("Recent-onset", "Long-standing"))
    addmap(48:49, 27, c("5.0-9.9%", ">=10.0%"))
    addmap(50:51, 28:29)
    addmap(53:55, 30:32)
    addmap(57:62, 33:38)
    precise_table2_map()
    table2_template <- fromJSON(file.path(root, "reference", "manuscript_tables.json"), simplifyVector = FALSE)[[2]]
    map2 <- table2_display_map(table2_template)
    continuous_variables <- "age_years"
    stopifnot(all(c(map$variable, map2$variable) %in% names(d)), !anyNA(d[unique(c(map$variable, map2$variable))]))
    write.csv(map, file.path(root, "reference", "baseline_row_map.csv"), row.names = FALSE)
    main <- fromJSON(file.path(root, "reference", "manuscript_tables.json"), simplifyVector = FALSE)[1:3]
    supp <- fromJSON(file.path(root, "reference", "supplement_tables.json"), simplifyVector = FALSE)[1:5]
    fmtp <- function(p) if (length(p) == 0 || !is.finite(p)) 
        "NE"
    else if (p < 0.001) 
        "<0.001"
    else sprintf("%.3f", p)
    fmtn <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)
    records <- list()
    key <- function(scope, tid, ri, ci) paste(scope, tid, ri, ci, sep = ":")
    for (scope in c("main", "supp")) for (tid in seq_along(if (scope == "main") main else supp)) {
        tab <- (if (scope == "main") 
            main
        else supp)[[tid]]
        for (rr in tab$rows) for (cc in rr$cells) {
            keep <- cc$column == 0 || (rr$row == 0 && !(scope == "supp" && tid >= 4 && cc$column %in% 1:4)) || (rr$row == 1 && grepl("SMD|P value", 
                cc$text))
            records[[key(scope, tid, rr$row, cc$column)]] <- list(scope = scope, table_index = tid, row_index = rr$row, col_index = cc$column, 
                span = cc$span, text = if (keep) cc$text else "", template_text = cc$text, kind = "structure", variable = "", level = "", 
                arm = "", population = "", adjustment = "", n_total = NA_integer_, n_nonmissing = NA_integer_, n_missing = NA_integer_, 
                n_event = NA_real_, percentage = NA_real_, mean = NA_real_, sd = NA_real_, smd = NA_real_, p = NA_real_, method = "", 
                status = "structure")
        }
    }
    put <- function(scope, tid, ri, ci, val, kind = "numeric") {
        k <- key(scope, tid, ri, ci)
        if (is.null(records[[k]])) 
            stop("Unrecognized source cell ", k)
        z <- records[[k]]
        if (is.list(val)) 
            z[names(val)] <- val
        else z$text <- as.character(val)
        z$kind <- kind
        if (z$status == "structure") 
            z$status <- "observed"
        records[[k]] <<- z
    }
    indicator <- function(x, lv) if (is.character(x) || is.factor(x)) 
        as.numeric(x == lv)
    else as.numeric(x == if (lv == "Yes") 1 else suppressWarnings(as.numeric(lv)))
    summary_cell <- function(z, v, lv, adj = "Crude", digits = 1L, show_n = FALSE, binary_n = FALSE) {
        if (!v %in% names(z)) 
            stop("Missing mapped variable ", v)
        x <- z[[v]]
        n <- nrow(z)
        ok <- !is.na(x)
        nn <- sum(ok)
        ww <- if (adj == "Crude") 
            rep(1, n)
        else z$IPTW
        ans <- list(variable = v, level = lv, adjustment = adj, n_total = n, n_nonmissing = nn, n_missing = n - nn, status = "observed")
        if (n == 0) 
            return(c(ans, list(text = if (v %in% continuous_variables || show_n) "—" else "0 (—)", n_event = if (v %in% continuous_variables || 
                show_n) NA_real_ else 0, percentage = NA_real_)))
        if (adj != "Crude" && any(!is.finite(ww) | ww <= 0)) 
            return(c(ans, list(text = "NE", status = "Invalid propensity estimation; weighted statistic unavailable")))
        if (v %in% continuous_variables || show_n) {
            if (nn == 0) 
                return(c(ans, list(text = "NE", status = "No observed values")))
            mn <- weighted.mean(x[ok], ww[ok])
            sdv <- if (adj == "Crude") 
                sd(x[ok])
            else sqrt(sum(ww[ok] * (x[ok] - mn)^2)/sum(ww[ok]))
            text <- paste0(formatC(mn, format = "f", digits = digits), " ± ", if (is.finite(sdv)) 
                formatC(sdv, format = "f", digits = digits)
            else "—")
            if (show_n && adj == "Crude") 
                text <- paste0(text, " (N = ", nn, ")")
            return(c(ans, list(text = text, mean = mn, sd = sdv)))
        }
        y <- indicator(x, lv)
        ne <- sum(y, na.rm = TRUE)
        pct <- 100 * sum(ww * ifelse(is.na(y), 0, y))/sum(ww)
        text <- if (adj != "Crude") 
            sprintf("%.1f%%", pct)
        else if (binary_n) 
            sprintf("%s/%s (%.1f%%)", fmtn(ne), fmtn(n), pct)
        else sprintf("%s (%.1f)", fmtn(ne), pct)
        c(ans, list(text = text, n_event = ne, percentage = pct))
    }
    smd_value <- function(z, v, lv, adj = "Crude") {
        x <- z[[v]]
        g <- z$group
        w <- if (adj == "Crude") 
            rep(1, nrow(z))
        else z$IPTW
        if (any(!is.finite(w)) || any(table(g) == 0)) 
            return(NA_real_)
        if (!v %in% continuous_variables) 
            x <- as.numeric(!is.na(indicator(x, lv)) & indicator(x, lv) == 1)
        ok <- !is.na(x)
        x <- x[ok]
        g <- g[ok]
        w <- w[ok]
        if (any(table(g) == 0)) 
            return(NA_real_)
        mu <- vapply(levels(g), function(a) weighted.mean(x[g == a], w[g == a]), 0)
        vv <- if (v %in% continuous_variables) 
            vapply(levels(g), function(a) var(x[g == a]), 0)
        else vapply(levels(g), function(a) {
            p <- mean(x[g == a])
            p * (1 - p)
        }, 0)
        den <- sqrt(mean(vv))
        if (!is.finite(den)) 
            return(NA_real_)
        if (den == 0) 
            return(if (abs(diff(mu)) < 1e-12) 0 else Inf)
        unname(abs(diff(mu))/den)
    }
    for (pop in c("Overall", "High risk", "Low risk")) for (arm in c("Intervention", "Control")) {
        z <- d[d$group == arm & (pop == "Overall" | d$screening_group == pop), , drop = FALSE]
        ci <- switch(pop, Overall = 1L, `High risk` = 4L, `Low risk` = 7L) + as.integer(arm == "Control")
        put("main", 1, 1, ci, list(text = sprintf("%s\n(N = %s)", arm, fmtn(nrow(z))), n_total = nrow(z), arm = arm, population = pop), 
            "header")
        for (i in seq_len(nrow(map))) {
            m <- map[i, ]
            if (m$row_index %in% c(4, 5) && pop != "Overall") 
                next
            val <- summary_cell(z, m$variable, m$level)
            val$arm <- arm
            val$population <- pop
            put("main", 1, m$row_index, ci, val)
        }
    }
    pp <- d[d$selected_matched_cohort == 1, , drop = FALSE]
    for (adj in c("Crude", "IPTW")) for (arm in c("Intervention", "Control")) {
        ci <- if (adj == "Crude") 
            1L
        else 5L
        ci <- ci + as.integer(arm == "Control")
        z <- pp[pp$group == arm, , drop = FALSE]
        ess <- if (adj == "Crude") 
            nrow(z)
        else sum(z$IPTW)^2/sum(z$IPTW^2)
        header <- if (adj == "Crude") 
            sprintf("%s\n(N = %s)", arm, fmtn(nrow(z)))
        else sprintf("%s\n(ESS = %.1f)", arm, ess)
        put("main", 2, 1, ci, list(text = header, n_total = nrow(z), arm = arm, adjustment = adj), "header")
        for (i in seq_len(nrow(map2))) {
            m <- map2[i, ]
            ri <- m$row_index
            val <- summary_cell(z, m$variable, m$level, adj)
            val$arm <- arm
            val$population <- "Monitored"
            put("main", 2, ri, ci, val)
            if (arm == "Intervention" && m$show_smd) {
                ss <- smd_value(pp, m$variable, m$level, adj)
                put("main", 2, ri, if (adj == "Crude") 
                  3
                else 7, list(text = if (is.na(ss)) "NE" else if (is.infinite(ss)) "Inf" else sprintf("%.3f", ss), smd = ss, variable = m$variable, 
                  level = m$level, adjustment = adj, method = "Absolute standardized difference; fixed unweighted pooled SD denominator"))
            }
        }
    }
    for (tid in 1:3) {
        base <- if (tid == 1) 
            d
        else d[!is.na(d$holter_within_1y) & d$holter_within_1y == if (tid == 2) 
            0
        else 1, , drop = FALSE]
        for (k in 1:3) for (arm in c("Intervention", "Control")) {
            st <- c("I", "II", "III")[k]
            z <- base[base$clinical_stratum == st & base$group == arm, , drop = FALSE]
            ci <- c(1, 4, 7)[k] + as.integer(arm == "Control")
            put("supp", tid, 1, ci, list(text = sprintf("%s\n(N = %s)", arm, fmtn(nrow(z))), n_total = nrow(z), arm = arm, population = paste("Stratum", 
                st)), "header")
            for (i in seq_len(nrow(map))) {
                m <- map[i, ]
                r <- m$row_index
                if (r %in% c(14, 15, 16)) 
                  next
                ri <- if (r <= 12) 
                  r
                else r - 4L
                val <- summary_cell(z, m$variable, m$level)
                val$arm <- arm
                val$population <- paste("Stratum", st)
                put("supp", tid, ri, ci, val)
            }
        }
    }
    fisher_global <- function(tab) {
        tab <- as.matrix(tab)
        tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
        if (nrow(tab) < 2 || ncol(tab) < 2) 
            return(list(p = NA_real_, method = "NE: no variation or fewer than two nonempty groups"))
        test <- tryCatch(fisher.test(tab, workspace = 2e+07), error = function(e) {
        })
        if (is.null(test)) 
            return(list(p = NA_real_, method = "Exact Fisher computation failed; P value not estimable"))
        list(p = unname(test$p.value), method = "Two-sided exact Fisher-Freeman-Halton test")
    }
    binary_p <- function(subs, v, lv = "Yes") {
        ns <- vapply(subs, function(z) sum(!is.na(z[[v]])), 0L)
        es <- vapply(subs, function(z) sum(indicator(z[[v]], lv), na.rm = TRUE), 0)
        fisher_global(rbind(es, ns - es))
    }
    category_p <- function(subs, vs) {
        cats <- lapply(subs, function(z) {
            a <- z[, vs, drop = FALSE]
            a <- a[complete.cases(a), , drop = FALSE]
            if (!nrow(a)) 
                character()
            else apply(a, 1, paste, collapse = "|")
        })
        lev <- sort(unique(unlist(cats)))
        if (!length(lev)) 
            return(list(p = NA_real_, method = "NE: no observed categories"))
        fisher_global(vapply(cats, function(x) as.numeric(table(factor(x, levels = lev))), numeric(length(lev))))
    }
    age_p <- function(subs) {
        xs <- lapply(subs, function(z) z[["age_years"]][is.finite(z[["age_years"]])])
        xs <- xs[lengths(xs) >= 2]
        if (length(xs) < 2) 
            return(list(p = NA_real_, method = "NE: fewer than two groups with at least two available ages"))
        age_data <- data.frame(age = unlist(xs, use.names = FALSE),
            comparison_group = factor(rep(seq_along(xs), lengths(xs))))
        test <- stats::anova(stats::lm(age ~ comparison_group, data = age_data))
        list(p = unname(test[["Pr(>F)"]][1L]), method = "Ordinary one-way ANOVA among groups with at least two available observations")
    }
    indvars <- c("indication_high_grade_av_block", "indication_complete_av_block", "indication_sick_sinus_syndrome",
        "indication_ventricular_arrhythmia", "other_cied_indication")
    for (tid in 4:5) {
        sel <- if (tid == 4) 
            d[["cied_indication"]] == 1
        else d$cied_implantation == 1
        base <- d[!is.na(sel) & sel, , drop = FALSE]
        conditional_row <- function(ri) if (tid == 5L && ri >= 8L) 
            ri + 3L
        else ri
        if (tid == 5L) 
            stopifnot("cied_type" %in% names(base), !anyNA(base$cied_type), all(base$cied_type %in% c("Pacemaker", "ICD")))
        arms <- c("Intervention", "Control", "Intervention", "Control")
        risks <- c("High risk", "High risk", "Low risk", "Low risk")
        subs <- lapply(1:4, function(k) base[base$group == arms[k] & base$screening_group == risks[k], , drop = FALSE])
        for (ci in 1:4) {
            z <- subs[[ci]]
            put("supp", tid, 0, ci, list(text = sprintf("%s:\n%s\n(N = %s)", ifelse(risks[ci] == "High risk", "Flagged by the AI-ECG screening algorithm", 
                "Not flagged by the AI-ECG screening algorithm"), arms[ci], fmtn(nrow(z))), n_total = nrow(z), arm = arms[ci], population = risks[ci]), 
                "header")
            put("supp", tid, 1, ci, summary_cell(z, "selected_matched_cohort", "Yes"))
            for (k in 1:5) put("supp", tid, k + 2, ci, summary_cell(z, indvars[k], "Yes"))
            if (tid == 5L) {
                put("supp", tid, 9, ci, summary_cell(z, "cied_type", "Pacemaker"))
                put("supp", tid, 10, ci, summary_cell(z, "cied_type", "ICD"))
            }
            for (i in seq_len(nrow(map))) {
                m <- map[i, ]
                r <- m$row_index
                if (!(r %in% c(8, 9, 11, 12)) && r < 18) 
                  next
                ri <- if (r %in% c(8, 9, 11, 12)) 
                  r + 1L
                else r - 3L
                put("supp", tid, conditional_row(ri), ci, summary_cell(z, m$variable, m$level))
            }
        }
        test_put <- function(ri, test, v) {
            put("supp", tid, conditional_row(ri), 5, list(text = fmtp(test$p), p = test$p, method = test$method, variable = v), "p")
        }
        if (tid == 5L) {
            ct <- category_p(subs, "cied_type")
            put("supp", tid, 8, 5, list(text = fmtp(ct$p), p = ct$p, method = ct$method, variable = "cied_type"), "p")
        }
        jobs <- list(`1` = c("selected_matched_cohort", "Yes"), `9` = c("male", "1"))
        for (i in seq_len(nrow(map))) {
            m <- map[i, ]
            r <- m$row_index
            if (r >= 21 && !r %in% c(38, 39, 41, 42, 45, 46, 48, 49)) 
                jobs[[as.character(r - 3L)]] <- c(m$variable, m$level)
        }
        for (ri in names(jobs)) test_put(as.integer(ri), binary_p(subs, jobs[[ri]][1], jobs[[ri]][2]), jobs[[ri]][1])
        test_put(10, age_p(subs), "age_years")
        test_put(11, binary_p(subs, "age_ge80", "1"), "age_ge80")
        catjobs <- list(`15` = "af_recency", `34` = "rbbb_recency", `37` = "lbbb_recency", `41` = "pvc_ecg_recency", `44` = "holter_pvc_burden_category")
        for (ri in names(catjobs)) test_put(as.integer(ri), category_p(subs, catjobs[[ri]]), catjobs[[ri]])
        for (k in 1:5) test_put(k + 2, binary_p(subs, indvars[k]), indvars[k])
    }
    precise_table3()
    allcells <- unname(records)
    flat <- fromJSON(toJSON(allcells, auto_unbox = TRUE, na = "null", null = "null", digits = 16))
    for (scope in c("main", "supp")) for (tid in seq_along(if (scope == "main") main else supp)) {
        tt <- (if (scope == "main") 
            main
        else supp)[[tid]]
        nc <- max(vapply(tt$rows[[1]]$cells, function(x) x$column + x$span, 0))
        nr <- length(tt$rows)
        grid <- matrix("", nr, nc)
        rr <- flat[flat$scope == scope & flat$table_index == tid, ]
        grid[cbind(rr$row_index + 1, rr$col_index + 1)] <- rr$text
        write.table(grid, file.path(out, sprintf("%s_Table_%d.csv", if (scope == "main") 
            "Main"
        else "Supplementary", if (scope == "main") 
            tid
        else tid + 2L)), sep = ",", row.names = FALSE, col.names = FALSE, na = "", qmethod = "double")
    }
    stopifnot(all(flat$n_event[is.finite(flat$n_event)] >= 0), all(flat$n_event[is.finite(flat$n_event)] <= flat$n_total[is.finite(flat$n_event)]))
}

precise_core_table8 <- function(root) {
    precise_publication_labels()
    out <- file.path(root, "tables")
    dir.create(out, showWarnings = FALSE, recursive = TRUE)
    d <- read.csv(file.path(root, "data/PRECISE-AI.csv"), check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
    stopifnot(!anyDuplicated(d$row_number), all(d$group %in% c("Intervention", "Control")), all(d$ai_ecg_screening_flag %in% 0:1), 
        !anyNA(d[c("ai_ecg_screening_flag", "high_risk_ecg_count_1y")]))
    strata <- THREE_GROUP_LABELS
    d$screening_history_group <- factor(ifelse(d$ai_ecg_screening_flag == 1, strata[1], ifelse(d$high_risk_ecg_count_1y > 0, strata[2], 
        strata[3])), levels = strata)
    map <- data.frame(row_index = integer(), variable = character(), level = character(), stringsAsFactors = FALSE)
    addmap <- function(rows, idx, levels = "1") {
        n <- length(rows)
        map <<- rbind(map, data.frame(row_index = rows, variable = precise_baseline_variables[rep(idx, length.out = n)], level = rep(levels, 
            length.out = n), stringsAsFactors = FALSE))
    }
    addmap(c(4, 5), 1, c("1", "0"))
    addmap(6, 2)
    addmap(8, 3)
    addmap(9, 4, "Continuous")
    addmap(c(11, 12), 5, c("0", "1"))
    addmap(14:16, 6, c("I", "II", "III"))
    addmap(c(19, 20), 7, c("Recent-onset", "Long-standing"))
    addmap(21:33, 8:20)
    addmap(35:36, 21:22)
    addmap(38:39, 23, c("Recent-onset", "Long-standing"))
    addmap(41:42, 24, c("Recent-onset", "Long-standing"))
    addmap(43, 25)
    addmap(45:46, 26, c("Recent-onset", "Long-standing"))
    addmap(48:49, 27, c("5.0-9.9%", ">=10.0%"))
    addmap(50:51, 28:29)
    addmap(53:55, 30:32)
    addmap(57:62, 33:38)
    reference_map <- read.csv(file.path(root, "reference/baseline_row_map.csv"), stringsAsFactors = FALSE)
    stopifnot(isTRUE(all.equal(map, reference_map, check.attributes = FALSE)), nrow(map) == 47L, !anyNA(d[unique(map$variable)]))
    template <- jsonlite::fromJSON(file.path(root, "reference/manuscript_tables.json"), simplifyVector = FALSE)[[1]]
    stopifnot(length(template$rows) == 63L)
    labels <- vapply(template$rows, function(r) {
        cc <- Filter(function(z) z$column == 0, r$cells)
        if (length(cc)) 
            cc[[1]]$text
        else ""
    }, character(1))
    column_strata <- rep(strata, each = 2)
    column_arms <- rep(c("Intervention", "Control"), 3)
    headers <- THREE_GROUP_LABELS
    grid <- matrix("", nrow = 63, ncol = 7)
    grid[, 1] <- labels
    grid[1, c(2, 4, 6)] <- headers
    fmtn <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
    for (ci in 1:6) {
        z <- d[d$screening_history_group == column_strata[ci] & d$group == column_arms[ci], , drop = FALSE]
        n <- nrow(z)
        grid[2, ci + 1] <- sprintf("%s\n(N = %s)", column_arms[ci], fmtn(n))
        for (i in seq_len(nrow(map))) {
            m <- map[i, ]
            x <- z[[m$variable]]
            nn <- sum(!is.na(x))
            cell <- list(row_index = m$row_index, column = ci, label = labels[m$row_index + 1], variable = m$variable, level = m$level, 
                stratum = column_strata[ci], arm = column_arms[ci], n_total = n, n_nonmissing = nn, n_missing = n - nn, n_event = NA_real_, 
                percentage = NA_real_, mean = NA_real_, sd = NA_real_)
            if (m$level == "Continuous") {
                cell$mean <- mean(x)
                cell$sd <- sd(x)
                cell$text <- sprintf("%.1f ± %.1f", cell$mean, cell$sd)
            }
            else {
                indicator <- if (is.character(x) || is.factor(x)) 
                  x == m$level
                else x == as.numeric(m$level)
                cell$n_event <- sum(indicator)
                cell$percentage <- 100 * cell$n_event/n
                cell$text <- sprintf("%s (%.1f)", fmtn(cell$n_event), cell$percentage)
            }
            grid[m$row_index + 1, ci + 1] <- cell$text
        }
    }
    physical_grid <- matrix("", nrow = 63, ncol = 9)
    physical_grid[, c(1, 2, 3, 5, 6, 8, 9)] <- grid
    file_base <- file.path(out, "Supplemental_Table_08_Risk_History_Baseline")
    write.table(physical_grid, paste0(file_base, ".csv"), sep = ",", row.names = FALSE, col.names = FALSE, na = "", qmethod = "double")
}

precise_table_templates <- list(main = list(list(labels = c("", "", "Key profile", "Stratification by AI-ECG", "Flagged by the AI-ECG screening algorithm, n (%)", 
    "Not flagged by the AI-ECG screening algorithm, n (%)", "  Holter monitoring within 1 year, n (%)", "Demographics and prespecified clinical condition", 
    "  Male, n (%)", "  Age (years, mean ± SD)", "Age group", "65-79 years, n (%)", "80-90 years, n (%)", "  Prespecified clinical condition subgroupϵ", 
    "Potentially suboptimal rhythm evaluation, n (%)", "Risk indicators supporting rhythm monitoring, n (%)", "No prespecified rhythm risk indicators, n (%)", 
    "Clinical history", "  Atrial fibrillation/flutter", "    Recent-onset (≤1 years), n (%)", "Long-standing (>1 years), n (%)", 
    "Hypertension, n (%)", "Diabetes mellitus, n (%)", "Heart failure, n (%)", "Chronic kidney disease, n (%)", "Chronic obstructive pulmonary disease, n (%)", 
    "Peripheral arterial disease, n (%)", "Valvular heart disease, n (%)", "Cerebrovascular disease, n (%)", "Obstructive sleep apnea, n (%)", 
    "Thyroid disease, n (%)", "Prior myocardial infarction, n (%)", "History of coronary artery bypass grafting, n (%)", "History of electrophysiology study or catheter ablation, n (%)", 
    "Past 12-lead ECG and Holter profile", "Mobitz type I AV block in 12-lead ECG or Holter, n (%)", "Borderline low heart rate in 12-lead ECG or Holter, n (%)", 
    "Complete right bundle branch block in 12-lead ECG", "Recent-onset (≤1 years), n (%)", "Long-standing (>1 years), n (%)", "Complete left bundle branch block in 12-lead ECG", 
    "Recent-onset (≤1 years), n (%)", "Long-standing (>1 years), n (%)", "Left fascicular block in 12-lead ECG, n (%)", "Premature ventricular complex in 12-lead ECG", 
    "Recent-onset (≤1 years), n (%)", "Long-standing (>1 years), n (%)", "Premature ventricular complex burden in Holter", "Premature ventricular complex burden 5.0-9.9%, n (%)", 
    "Premature ventricular complex burden ≥10.0%, n (%)", "Non-sustained ventricular tachycardia in Holter, n (%)", "Borderline pause or asymptomatic pause in Holter, n (%)", 
    "PR interval group", "Marked PR prolongation (≥250) in 12-lead ECG, n (%)", "Mild PR prolongation (200-249) in 12-lead ECG, n (%)", 
    "Documeted PR prolongation (≥200) in Holter, n (%)", "Medication history", "Beta blocker, n (%)", "Non-dihydropyridine calcium channel blockers, n (%)", 
    "Class IA antiarrhythmia drug, n (%)", "Class IC antiarrhythmia drug, n (%)", "Class III antiarrhythmia drug, n (%)", "Digoxin, n (%)"), 
    ncol = 9, headers = list(c("", "All patients", "", "", "Flagged by the AI-ECG screening algorithm", "", "", "Not flagged by the AI-ECG screening algorithm", 
        ""), c("", "", "", "", "", "", "", "", ""))), list(labels = c("", "", "Key profile", "Stratification by AI-ECG", "Flagged by the AI-ECG screening algorithm, n (%)", 
    "Not flagged by the AI-ECG screening algorithm, n (%)", "Demographics", "  Male, n (%)", "  Age (years, mean ± SD)", "Age group", 
    "65-79 years, n (%)", "80-90 years, n (%)", "Clinical history", "  Atrial fibrillation/flutter, n (%)", "Hypertension, n (%)", 
    "Diabetes mellitus, n (%)", "Heart failure, n (%)", "Chronic kidney disease, n (%)", "Chronic obstructive pulmonary disease, n (%)", 
    "Peripheral arterial disease, n (%)", "Valvular heart disease, n (%)", "Cerebrovascular disease, n (%)", "Obstructive sleep apnea, n (%)", 
    "Thyroid disease, n (%)", "Prior myocardial infarction, n (%)", "History of coronary artery bypass grafting, n (%)", "History of electrophysiology study or catheter ablation, n (%)", 
    "Past 12-lead ECG and Holter profile", "Mobitz type I AV block in 12-lead ECG or Holter, n (%)", "Borderline low heart rate in 12-lead ECG or Holter, n (%)", 
    "Complete right bundle branch block in 12-lead ECG", "Recent-onset (≤1 years), n (%)", "Long-standing (>1 years), n (%)", "Complete left bundle branch block in 12-lead ECG", 
    "Recent-onset (≤1 years), n (%)", "Long-standing (>1 years), n (%)", "Left fascicular block in 12-lead ECG, n (%)", "Premature ventricular complex in 12-lead ECG", 
    "Recent-onset (≤1 years), n (%)", "Long-standing (>1 years), n (%)", "Premature ventricular complex burden in Holter", "Premature ventricular complex burden 5.0-9.9%, n (%)", 
    "Premature ventricular complex burden ≥10.0%, n (%)", "Non-sustained ventricular tachycardia in Holter, n (%)", "Borderline pause or asymptomatic pause in Holter, n (%)", 
    "PR interval group", "Marked PR prolongation (≥250) in 12-lead ECG, n (%)", "Mild PR prolongation (200-249) in 12-lead ECG, n (%)", 
    "Documeted PR prolongation (≥200) in Holter, n (%)", "Medication history", "Beta blocker, n (%)", "Non-dihydropyridine calcium channel blockers, n (%)", 
    "Class IA antiarrhythmia drug, n (%)", "Class IC antiarrhythmia drug, n (%)", "Class III antiarrhythmia drug, n (%)", "Digoxin, n (%)"), 
    ncol = 8, headers = list(c("", "Before IPTW", "", "", "", "After IPTW", "", ""), c("", "", "", "SMD", "", "", "", "SMD"))), list(labels = c("", 
    "", "New-onset atrial fibrillation, n/N (%)", "Average heart rate (bpm), mean ± SD (N)", "Maximal heart rate (bpm), mean ± SD (N)", 
    "Minimal heart rate (bpm), mean ± SD (N)", "Longest RR interval (second), mean ± SD (N)", "Premature atrial complex burden (%), mean ± SD (N)", 
    "Premature ventricular complex burden (%), mean ± SD (N)", "New-onset non-sustained ventricular tachycardia, n/N (%)", "New-onset bundle branch block, n/N (%)", 
    "New-onset PR prolongation, n/N (%)", "New-onset Mobitz type I AV block, n/N (%)"), ncol = 8, headers = list(c("", "Before IPTW", 
    "", "", "", "After IPTW", "", ""), c("", "", "", "P value", "", "", "", "P value")))), supp = list(list(labels = c("", "", "Key profile", 
    "  Stratification by AI-ECG", "    Flagged by the AI-ECG screening algorithm, n (%)", "    Not flagged by the AI-ECG screening algorithm, n (%)", 
    "  Holter monitoring within 1 year, n (%)", "Demographics and prespecified clinical condition", "  Male, n (%)", "  Age (years, mean ± SD)", 
    "Age group", "65-79 years, n (%)", "80-90 years, n (%)", "Clinical history", "  Atrial fibrillation/flutter", "    Recent-onset (≤1 years), n (%)", 
    "    Long-standing (>1 years), n (%)", "  Hypertension, n (%)", "  Diabetes mellitus, n (%)", "  Heart failure, n (%)", "  Chronic kidney disease, n (%)", 
    "  Chronic obstructive pulmonary disease, n (%)", "  Peripheral arterial disease, n (%)", "  Valvular heart disease, n (%)", "  Cerebrovascular disease, n (%)", 
    "  Obstructive sleep apnea, n (%)", "  Thyroid disease, n (%)", "  Prior myocardial infarction, n (%)", "  History of coronary artery bypass grafting, n (%)", 
    "  History of electrophysiology study or catheter ablation, n (%)", "Past 12-lead ECG and Holter profile", "  Mobitz type I AV block in 12-lead ECG or Holter, n (%)", 
    "  Borderline low heart rate in 12-lead ECG or Holter, n (%)", "  Complete right bundle branch block in 12-lead ECG", "    Recent-onset (≤1 years), n (%)", 
    "    Long-standing (>1 years), n (%)", "  Complete left bundle branch block in 12-lead ECG", "    Recent-onset (≤1 years), n (%)", 
    "    Long-standing (>1 years), n (%)", "  Left fascicular block in 12-lead ECG, n (%)", "  Premature ventricular complex in 12-lead ECG", 
    "    Recent-onset (≤1 years), n (%)", "    Long-standing (>1 years), n (%)", "  Premature ventricular complex burden in Holter", 
    "    burden 5.0-9.9%, n (%)", "    burden ≥10.0%, n (%)", "  Non-sustained ventricular tachycardia in Holter, n (%)", "  Borderline pause or asymptomatic pause in Holter, n (%)", 
    "  PR interval group", "    Marked PR prolongation (≥250) in 12-lead ECG, n (%)", "    Mild PR prolongation (200-249) in 12-lead ECG, n (%)", 
    "    Documeted PR prolongation (≥200) in Holter, n (%)", "Medication history", "  Beta blocker, n (%)", "  Non-dihydropyridine calcium channel blockers, n (%)", 
    "  Class IA antiarrhythmia drug, n (%)", "  Class IC antiarrhythmia drug, n (%)", "  Class III antiarrhythmia drug, n (%)", "  Digoxin, n (%)"), 
    ncol = 9, headers = list(c("", "Potentially suboptimal rhythm evaluation", "", "", "Risk indicators supporting rhythm monitoring", 
        "", "", "No prespecified rhythm risk indicators", ""), c("", "", "", "", "", "", "", "", ""))), list(labels = c("", "", "Key profile", 
    "  Stratification by AI-ECG", "    Flagged by the AI-ECG screening algorithm, n (%)", "    Not flagged by the AI-ECG screening algorithm, n (%)", 
    "  Holter monitoring within 1 year, n (%)", "Demographics and prespecified clinical condition", "  Male, n (%)", "  Age (years, mean ± SD)", 
    "Age group", "65-79 years, n (%)", "80-90 years, n (%)", "Clinical history", "  Atrial fibrillation/flutter", "    Recent-onset (≤1 years), n (%)", 
    "    Long-standing (>1 years), n (%)", "  Hypertension, n (%)", "  Diabetes mellitus, n (%)", "  Heart failure, n (%)", "  Chronic kidney disease, n (%)", 
    "  Chronic obstructive pulmonary disease, n (%)", "  Peripheral arterial disease, n (%)", "  Valvular heart disease, n (%)", "  Cerebrovascular disease, n (%)", 
    "  Obstructive sleep apnea, n (%)", "  Thyroid disease, n (%)", "  Prior myocardial infarction, n (%)", "  History of coronary artery bypass grafting, n (%)", 
    "  History of electrophysiology study or catheter ablation, n (%)", "Past 12-lead ECG and Holter profile", "  Mobitz type I AV block in 12-lead ECG or Holter, n (%)", 
    "  Borderline low heart rate in 12-lead ECG or Holter, n (%)", "  Complete right bundle branch block in 12-lead ECG", "    Recent-onset (≤1 years), n (%)", 
    "    Long-standing (>1 years), n (%)", "  Complete left bundle branch block in 12-lead ECG", "    Recent-onset (≤1 years), n (%)", 
    "    Long-standing (>1 years), n (%)", "  Left fascicular block in 12-lead ECG, n (%)", "  Premature ventricular complex in 12-lead ECG", 
    "    Recent-onset (≤1 years), n (%)", "    Long-standing (>1 years), n (%)", "  Premature ventricular complex burden in Holter", 
    "    burden 5.0-9.9%, n (%)", "    burden ≥10.0%, n (%)", "  Non-sustained ventricular tachycardia in Holter, n (%)", "  Borderline pause or asymptomatic pause in Holter, n (%)", 
    "  PR interval group", "    Marked PR prolongation (≥250) in 12-lead ECG, n (%)", "    Mild PR prolongation (200-249) in 12-lead ECG, n (%)", 
    "    Documeted PR prolongation (≥200) in Holter, n (%)", "Medication history", "  Beta blocker, n (%)", "  Non-dihydropyridine calcium channel blockers, n (%)", 
    "  Class IA antiarrhythmia drug, n (%)", "  Class IC antiarrhythmia drug, n (%)", "  Class III antiarrhythmia drug, n (%)", "  Digoxin, n (%)"), 
    ncol = 9, headers = list(c("", "Potentially suboptimal rhythm evaluation", "", "", "Risk indicators supporting rhythm monitoring", 
        "", "", "No prespecified rhythm risk indicators", ""), c("", "", "", "", "", "", "", "", ""))), list(labels = c("", "", "Key profile", 
    "  Stratification by AI-ECG", "    Flagged by the AI-ECG screening algorithm, n (%)", "    Not flagged by the AI-ECG screening algorithm, n (%)", 
    "  Holter monitoring within 1 year, n (%)", "Demographics and prespecified clinical condition", "  Male, n (%)", "  Age (years, mean ± SD)", 
    "Age group", "65-79 years, n (%)", "80-90 years, n (%)", "Clinical history", "  Atrial fibrillation/flutter", "    Recent-onset (≤1 years), n (%)", 
    "    Long-standing (>1 years), n (%)", "  Hypertension, n (%)", "  Diabetes mellitus, n (%)", "  Heart failure, n (%)", "  Chronic kidney disease, n (%)", 
    "  Chronic obstructive pulmonary disease, n (%)", "  Peripheral arterial disease, n (%)", "  Valvular heart disease, n (%)", "  Cerebrovascular disease, n (%)", 
    "  Obstructive sleep apnea, n (%)", "  Thyroid disease, n (%)", "  Prior myocardial infarction, n (%)", "  History of coronary artery bypass grafting, n (%)", 
    "  History of electrophysiology study or catheter ablation, n (%)", "Past 12-lead ECG and Holter profile", "  Mobitz type I AV block in 12-lead ECG or Holter, n (%)", 
    "  Borderline low heart rate in 12-lead ECG or Holter, n (%)", "  Complete right bundle branch block in 12-lead ECG", "    Recent-onset (≤1 years), n (%)", 
    "    Long-standing (>1 years), n (%)", "  Complete left bundle branch block in 12-lead ECG", "    Recent-onset (≤1 years), n (%)", 
    "    Long-standing (>1 years), n (%)", "  Left fascicular block in 12-lead ECG, n (%)", "  Premature ventricular complex in 12-lead ECG", 
    "    Recent-onset (≤1 years), n (%)", "    Long-standing (>1 years), n (%)", "  Premature ventricular complex burden in Holter", 
    "    burden 5.0-9.9%, n (%)", "    burden ≥10.0%, n (%)", "  Non-sustained ventricular tachycardia in Holter, n (%)", "  Borderline pause or asymptomatic pause in Holter, n (%)", 
    "  PR interval group", "    Marked PR prolongation (≥250) in 12-lead ECG, n (%)", "    Mild PR prolongation (200-249) in 12-lead ECG, n (%)", 
    "    Documeted PR prolongation (≥200) in Holter, n (%)", "Medication history", "  Beta blocker, n (%)", "  Non-dihydropyridine calcium channel blockers, n (%)", 
    "  Class IA antiarrhythmia drug, n (%)", "  Class IC antiarrhythmia drug, n (%)", "  Class III antiarrhythmia drug, n (%)", "  Digoxin, n (%)"), 
    ncol = 9, headers = list(c("", "Potentially suboptimal rhythm evaluation", "", "", "Risk indicators supporting rhythm monitoring", 
        "", "", "No prespecified rhythm risk indicators", ""), c("", "", "", "", "", "", "", "", ""))), list(labels = c("", "Protocol Holter, n (%)", 
    "CIED indication", "  High-grade atrioventricular block, n (%)", "  Complete atrioventricular block, n (%)", "  Sick sinus syndrome, n (%)", 
    "  Ventricular arrhythmia, n (%)", "  Other, n (%)", "Demographics and prespecified clinical condition", "  Male, n (%)", "  Age (years, mean ± SD)", 
    "Age group", "65-79 years, n (%)", "80-90 years, n (%)", "Clinical history", "  Atrial fibrillation/flutter", "    Recent-onset (≤1 years), n (%)", 
    "    Long-standing (>1 years), n (%)", "  Hypertension, n (%)", "  Diabetes mellitus, n (%)", "  Heart failure, n (%)", "  Chronic kidney disease, n (%)", 
    "  Chronic obstructive pulmonary disease, n (%)", "  Peripheral arterial disease, n (%)", "  Valvular heart disease, n (%)", "  Cerebrovascular disease, n (%)", 
    "  Obstructive sleep apnea, n (%)", "  Thyroid disease, n (%)", "  Prior myocardial infarction, n (%)", "  History of coronary artery bypass grafting, n (%)", 
    "  History of electrophysiology study or catheter ablation, n (%)", "Past 12-lead ECG and Holter profile", "  Mobitz type I AV block in 12-lead ECG or Holter, n (%)", 
    "  Borderline low heart rate in 12-lead ECG or Holter, n (%)", "  Complete right bundle branch block in 12-lead ECG", "    Recent-onset (≤1 years), n (%)", 
    "    Long-standing (>1 years), n (%)", "  Complete left bundle branch block in 12-lead ECG", "    Recent-onset (≤1 years), n (%)", 
    "    Long-standing (>1 years), n (%)", "  Left fascicular block in 12-lead ECG, n (%)", "  Premature ventricular complex in 12-lead ECG", 
    "    Recent-onset (≤1 years), n (%)", "    Long-standing (>1 years), n (%)", "  Premature ventricular complex burden in Holter", 
    "    burden 5.0-9.9%, n (%)", "    burden ≥10.0%, n (%)", "  Non-sustained ventricular tachycardia in Holter, n (%)", "  Borderline pause or asymptomatic pause in Holter, n (%)", 
    "  PR interval group", "    Marked PR prolongation (≥250) in 12-lead ECG, n (%)", "    Mild PR prolongation (200-249) in 12-lead ECG, n (%)", 
    "    Documeted PR prolongation (≥200) in Holter, n (%)", "Medication history", "  Beta blocker, n (%)", "  Non-dihydropyridine calcium channel blockers, n (%)", 
    "  Class IA antiarrhythmia drug, n (%)", "  Class IC antiarrhythmia drug, n (%)", "  Class III antiarrhythmia drug, n (%)", "  Digoxin, n (%)"), 
    ncol = 6, headers = list(c("", "Flagged by the AI-ECG screening algorithm:\nIntervention", "Flagged by the AI-ECG screening algorithm:\nControl", 
        "Not flagged by the AI-ECG screening algorithm:\nIntervention", "Not flagged by the AI-ECG screening algorithm:\nControl", 
        "P value"), c("Protocol Holter, n (%)", "", "", "", "", ""))), list(labels = c("", "Protocol Holter, n (%)", "CIED indication", 
    "  High-grade atrioventricular block, n (%)", "  Complete atrioventricular block, n (%)", "  Sick sinus syndrome, n (%)", "  Ventricular arrhythmia, n (%)", 
    "  Other, n (%)†", "CIED type", "  Pacemaker, n (%)", "  ICD, n (%)", "Demographics and prespecified clinical condition", 
    "  Male, n (%)", "  Age (years, mean ± SD)", "Age group", "65-79 years, n (%)", "80-90 years, n (%)", "Clinical history", "  Atrial fibrillation/flutter", 
    "    Recent-onset (≤1 years), n (%)", "    Long-standing (>1 years), n (%)", "  Hypertension, n (%)", "  Diabetes mellitus, n (%)", 
    "  Heart failure, n (%)", "  Chronic kidney disease, n (%)", "  Chronic obstructive pulmonary disease, n (%)", "  Peripheral arterial disease, n (%)", 
    "  Valvular heart disease, n (%)", "  Cerebrovascular disease, n (%)", "  Obstructive sleep apnea, n (%)", "  Thyroid disease, n (%)", 
    "  Prior myocardial infarction, n (%)", "  History of coronary artery bypass grafting, n (%)", "  History of electrophysiology study or catheter ablation, n (%)", 
    "Past 12-lead ECG and Holter profile", "  Mobitz type I AV block in 12-lead ECG or Holter, n (%)", "  Borderline low heart rate in 12-lead ECG or Holter, n (%)", 
    "  Complete right bundle branch block in 12-lead ECG", "    Recent-onset (≤1 years), n (%)", "    Long-standing (>1 years), n (%)", 
    "  Complete left bundle branch block in 12-lead ECG", "    Recent-onset (≤1 years), n (%)", "    Long-standing (>1 years), n (%)", 
    "  Left fascicular block in 12-lead ECG, n (%)", "  Premature ventricular complex in 12-lead ECG", "    Recent-onset (≤1 years), n (%)", 
    "    Long-standing (>1 years), n (%)", "  Premature ventricular complex burden in Holter", "    burden 5.0-9.9%, n (%)", "    burden ≥10.0%, n (%)", 
    "  Non-sustained ventricular tachycardia in Holter, n (%)", "  Borderline pause or asymptomatic pause in Holter, n (%)", "  PR interval group", 
    "    Marked PR prolongation (≥250) in 12-lead ECG, n (%)", "    Mild PR prolongation (200-249) in 12-lead ECG, n (%)", "    Documeted PR prolongation (≥200) in Holter, n (%)", 
    "Medication history", "  Beta blocker, n (%)", "  Non-dihydropyridine calcium channel blockers, n (%)", "  Class IA antiarrhythmia drug, n (%)", 
    "  Class IC antiarrhythmia drug, n (%)", "  Class III antiarrhythmia drug, n (%)", "  Digoxin, n (%)"), ncol = 6, headers = list(c("", 
    "Flagged by the AI-ECG screening algorithm:\nIntervention", "Flagged by the AI-ECG screening algorithm:\nControl", 
    "Not flagged by the AI-ECG screening algorithm:\nIntervention", "Not flagged by the AI-ECG screening algorithm:\nControl", 
    "P value"), c("Protocol Holter, n (%)", "", "", "", "", "")))))

prepare_analysis_stage <- function(d, stage) {
    for (module in c("core", "safety", "recent", "timing")) {
        p <- file.path(stage, module)
        for (folder in c("data", "analysis", "tables", "figures", "reference")) dir.create(file.path(p, folder), recursive = TRUE, 
            showWarnings = FALSE)
        write.csv(d, file.path(p, "data", "PRECISE-AI.csv"), row.names = FALSE, na = "")
        if (module %in% c("core", "safety"))
            writeLines(precise_propensity_variables, file.path(p, "data/ps_variables.txt"))
    }
    make_template <- function(t) {
        list(rows = lapply(seq_along(t$labels), function(i) {
            values <- rep("", t$ncol)
            values[1] <- t$labels[i]
            if (i <= 2) values <- t$headers[[i]]
            list(row = i - 1L, cells = lapply(seq_len(t$ncol), function(j) list(column = j - 1L, span = 1L, text = values[j])))
        }))
    }
    for (scope in c("main", "supp")) jsonlite::write_json(lapply(precise_table_templates[[scope]], make_template), file.path(stage, 
        "core/reference", if (scope == "main") 
            "manuscript_tables.json"
        else "supplement_tables.json"), auto_unbox = TRUE)
    endpoints <- data.frame(outcome = c("cied_implantation", precise_secondary_outcomes, "monitor_adverse_event", "cied_adverse_event"), 
        label = c("CIED implantation", "High-grade atrioventricular block", "Complete atrioventricular block", "Sick sinus syndrome", 
            "Ventricular arrhythmia", "Adverse events related to continuous ECG monitoring", "CIED-related complications"))
    write.csv(endpoints, file.path(stage, "safety/data/Fig_03_endpoint_labels.csv"), row.names = FALSE)
    invisible(stage)
}

run_analysis_engines <- function(d, stage) {
    prepare_analysis_stage(d, stage)
    core <- file.path(stage, "core")
    safety <- file.path(stage, "safety")
    recent <- file.path(stage, "recent")
    timing <- file.path(stage, "timing")
    precise_core_analysis(core)
    precise_core_diagnostics(core)
    precise_core_yields(core)
    precise_core_figures(core)
    precise_core_risk_window(core)
    precise_core_risk_history(core)
    precise_core_tables(core)
    precise_core_table8(core)
    precise_figure2(safety)
    precise_figure3_analysis(safety)
    precise_figure3_plot(safety)
    tab1 <- as.matrix(read.csv(file.path(core, "tables/Main_Table_1.csv"), header = FALSE, check.names = FALSE, stringsAsFactors = FALSE, 
        na.strings = character()))
    tab1[is.na(tab1)] <- ""
    template <- list(rows = lapply(seq_len(nrow(tab1)), function(i) list(cells = as.list(unname(tab1[i, ])))))
    jsonlite::write_json(template, file.path(recent, "reference/table1_template.json"), auto_unbox = TRUE)
    precise_recent_analysis(recent)
    precise_recent_plot(recent)
    precise_recent_table(recent)
    precise_primary_timing(timing)
    precise_primary_timing_plot(timing)
    precise_sss_timing(timing)
    precise_sss_timing_plot(timing)
    precise_time_distribution(timing)
    invisible(stage)
}

