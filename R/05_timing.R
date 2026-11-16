precise_exact_test <- function (envir = parent.frame()) 
{
    evalq({
        exact_rank_sum <- function(x, y) {
            stopifnot(length(x) > 0L, length(y) > 0L, all(is.finite(c(x, y))))
            z <- c(x, y)
            n <- length(z)
            nx <- length(x)
            ny <- length(y)
            ranks2 <- as.integer(round(2 * rank(z, ties.method = "average")))
            if (nx <= ny) {
                m <- nx
                selected <- seq_len(nx)
            }
            else {
                m <- ny
                selected <- nx + seq_len(ny)
            }
            observed <- sum(ranks2[selected])
            expected <- m * (n + 1)
            maxsum <- sum(sort(ranks2, decreasing = TRUE)[seq_len(m)])
            dp <- matrix(0, nrow = m + 1L, ncol = maxsum + 1L)
            dp[1L, 1L] <- 1
            used <- 0L
            for (r in ranks2) {
                for (k in seq.int(min(m, used + 1L), 1L)) {
                  ix <- which(dp[k, seq_len(maxsum - r + 1L)] > 0)
                  if (length(ix)) 
                    dp[k + 1L, ix + r] <- dp[k + 1L, ix + r] + dp[k, ix]
                }
                used <- used + 1L
            }
            ways <- dp[m + 1L, ]
            total <- choose(n, m)
            stopifnot(abs(sum(ways) - total) <= 1e-12 * total)
            sums2 <- seq.int(0, maxsum)
            tail <- abs(sums2 - expected) >= abs(observed - expected)
            p <- sum(ways[tail])/total
            list(p.value = p)
        }
    }, envir)
} 

open_figure_device <- function(path, width_mm, height_mm) {
    width_in <- width_mm/25.4
    height_in <- height_mm/25.4
    if (capabilities("cairo")) {
        grDevices::cairo_pdf(path, width = width_in, height = height_in, family = "sans", pointsize = 12)
    }
    else {
        grDevices::pdf(path, width = width_in, height = height_in, family = "Helvetica", pointsize = 12, useDingbats = FALSE)
    }
}

precise_primary_timing <- function (root) 
{
    suppressPackageStartupMessages({
        library(survival)
    })
    dir.create(file.path(root, "analysis"), showWarnings = FALSE, recursive = TRUE)
    precise_core_helpers()
    precise_screening_labels()
    d <- read.csv(file.path(root, "data/PRECISE-AI.csv"), check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
    stopifnot(!anyNA(d$group), !anyNA(d$screening_group), all(d$group %in% c("Intervention", "Control")), all(d$screening_group %in% c("High risk", "Low risk")), all(d$ai_ecg_screening_flag %in% 
        0:1), all(d$screening_group == ifelse(d$ai_ecg_screening_flag == 1, "High risk", "Low risk")), all(is.na(d$cied_implantation_days) | (is.finite(d$cied_implantation_days) & 
        d$cied_implantation_days >= 0)))
    d <- prepare_analysis_factors(d)
    d$TIME <- ifelse(is.na(d$cied_implantation_days), 180, pmin(d$cied_implantation_days, 180))
    d$EVENT <- as.integer(!is.na(d$cied_implantation_days) & d$cied_implantation_days <= 180)
    cutoffs <- seq(30, 180, by = 30)
    effects <- interactions <- list()
    for (t in cutoffs) {
        d$EVENT_BY_TIME <- as.integer(!is.na(d$cied_implantation_days) & d$cied_implantation_days <= t)
        br0 <- block_rule(d, "EVENT_BY_TIME")
        br1 <- block_rule(d, "EVENT_BY_TIME", "screening_group")
        i <- interaction(d, "EVENT_BY_TIME", "screening_group", "ITT", use_firth = br1$firth, block_id = paste0("day", t, "_risk"), trigger_levels = br1$trigger_levels)
        i$day <- t
        interactions[[length(interactions) + 1L]] <- i
        for (pop in c("Overall", "High risk", "Low risk")) {
            z <- if (pop == "Overall") 
                d
            else d[d$screening_group == pop, , drop = FALSE]
            br <- if (pop == "Overall") 
                br0
            else br1
            a <- effect(z, "EVENT_BY_TIME", "ITT", sg = if (pop == "Overall") 
                "Overall"
            else "screening_group", lv = if (pop == "Overall") 
                "All"
            else pop, use_firth = br$firth, block_id = paste0("day", t, if (pop == "Overall") 
                "_overall"
            else "_risk"), trigger_levels = br$trigger_levels)
            a$day <- t
            a$population_label <- screening_label(pop)
            a$outcome_label <- paste("CIED implantation by day", t)
            a$interaction_p <- i$p
            a$interaction_method <- i$method
            effects[[length(effects) + 1L]] <- a
        }
    }
    effects <- do.call(rbind, effects)
    interactions <- do.call(rbind, interactions)
    effects <- effects[c("day", "population_label", setdiff(names(effects), c("day", "population_label")))]
    interactions <- interactions[c("day", setdiff(names(interactions), "day"))]
    write.csv(effects, file.path(root, "analysis/Supplemental_Fig_13_effects.csv"), row.names = FALSE, na = "")
    write.csv(interactions, file.path(root, "analysis/Supplemental_Fig_13_interactions.csv"), row.names = FALSE, na = "")
    km <- risk <- list()
    for (pop in c("Overall", "High risk", "Low risk")) {
        zz <- if (pop == "Overall") 
            d
        else d[d$screening_group == pop, , drop = FALSE]
        pl <- screening_label(pop)
        for (arm in c("Intervention", "Control")) {
            z <- zz[zz$group == arm, , drop = FALSE]
            f <- survfit(Surv(TIME, EVENT) ~ 1, data = z, conf.int = 0.95, conf.type = "log-log")
            n <- nrow(z)
            baseline <- data.frame(population = pl, arm = arm, time = 0, survival = 1, lower_survival = 1, upper_survival = 1, cumulative_incidence = 0, 
                lower_incidence = 0, upper_incidence = 0, n_risk = n, n_event = 0, n_censor = 0)
            coord <- data.frame(population = pl, arm = arm, time = f$time, survival = f$surv, lower_survival = f$lower, upper_survival = f$upper, 
                cumulative_incidence = 1 - f$surv, lower_incidence = 1 - f$upper, upper_incidence = 1 - f$lower, n_risk = f$n.risk, n_event = f$n.event, 
                n_censor = f$n.censor)
            stopifnot(all(is.finite(coord$survival)), all(is.finite(coord$lower_survival)), all(is.finite(coord$upper_survival)))
            if (tail(coord$time, 1L) < 180) {
                last <- tail(coord, 1L)
                last$time <- 180
                last$n_risk <- sum(z$TIME >= 180)
                last$n_event <- 0
                last$n_censor <- sum(z$TIME == 180 & z$EVENT == 0)
                coord <- rbind(coord, last)
            }
            km[[length(km) + 1L]] <- rbind(baseline, coord)
            for (t in c(0, cutoffs)) {
                risk[[length(risk) + 1L]] <- data.frame(population = pl, arm = arm, day = t, n_risk = sum(z$TIME >= t), cumulative_events = sum(z$EVENT == 
                  1 & z$TIME <= t), cumulative_censors = sum(z$EVENT == 0 & z$TIME <= t), n_total = n)
            }
        }
    }
    km <- do.call(rbind, km)
    risk <- do.call(rbind, risk)
    write.csv(km, file.path(root, "analysis/Supplemental_Fig_13_km_coordinates.csv"), row.names = FALSE, na = "")
    write.csv(risk, file.path(root, "analysis/Supplemental_Fig_13_km_at_risk.csv"), row.names = FALSE, na = "")
    stopifnot(nrow(effects) == 18L, nrow(interactions) == 6L, all(diff(effects$event1[effects$population_label == "Overall"]) >= 0), all(diff(effects$event0[effects$population_label == 
        "Overall"]) >= 0), all(effects$status %in% c("OK", "No events in either arm")), 
        all(interactions$status %in% c("OK", "No events across subgroup block", "No events in either arm within a subgroup", "Ordinary interaction NE: nonfinite MLE in treatment-by-subgroup outcome cells")))
    for (pop in unique(effects$population_label)) {
        a <- effects[effects$population_label == pop, ]
        stopifnot(length(unique(a$n1)) == 1L, length(unique(a$n0)) == 1L, all(diff(a$event1) >= 0), all(diff(a$event0) >= 0))
    }
} 

precise_primary_timing_plot <- function (root) 
{
    precise_screening_labels()
    E <- read.csv(file.path(root, "analysis/Supplemental_Fig_13_effects.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    IT <- read.csv(file.path(root, "analysis/Supplemental_Fig_13_interactions.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    K <- read.csv(file.path(root, "analysis/Supplemental_Fig_13_km_coordinates.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    NR <- read.csv(file.path(root, "analysis/Supplemental_Fig_13_km_at_risk.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    E$population_label <- screening_label(E$population_label)
    K$population <- screening_label(K$population)
    NR$population <- screening_label(NR$population)
    out <- file.path(root, "figures")
    dir.create(out, showWarnings = FALSE, recursive = TRUE)
    H <- 194
    W <- 180
    PITCH <- 3.8
    GAP <- 0.8
    BLUE <- "#4B638F"
    RED <- "#C75D47"
    INK <- "#303135"
    RULE <- "#7B7D80"
    GRAY <- "#F1F2F4"
    COL <- c(Intervention = RED, Control = BLUE)
    num <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
    put <- function(x, y, s, size = 5.8, bold = FALSE, ha = "left", color = INK, srt = 0) {
        text(x, y, labels = s, cex = size/12, font = if (bold) 
            2
        else 1, adj = c(switch(ha, left = 0, center = 0.5, right = 1), 0.5), col = color, srt = srt)
    }
    line <- function(x1, y1, x2, y2, col = RULE, lwd = 0.65, lty = 1) segments(x1, y1, x2, y2, col = col, lwd = lwd, lty = lty)
    arrow <- function(x1, y1, x2, y2, col = BLUE, lwd = 1) arrows(x1, y1, x2, y2, length = 0.035, angle = 25, code = 2, col = col, lwd = lwd)
    pstr <- function(p, method = "") if (!is.finite(p)) 
        "NE"
    else paste0(if (p < 0.001) 
        "<0.001"
    else sprintf("%.3f", p), if (grepl("Firth", method)) 
        "†"
    else "")
    armstr <- function(r, k) sprintf("%s/%s (%.1f%%)", num(r[[paste0("event", k)]]), num(r[[paste0("n", k)]]), 100 * r[[paste0("event", k)]]/r[[paste0("n", 
        k)]])
    stepxy <- function(x, y) list(x = c(x[1], rep(x[-1], each = 2)), y = c(rep(y[-length(y)], each = 2), tail(y, 1)))
    X <- list(label = 3.4, inter = 60, control = 86, left = 102, right = 129, or = 144.6, p = 159.5, interaction = 171.7)
    tx <- function(v) X$left + (log(v) - log(0.25))/(log(32) - log(0.25)) * (X$right - X$left)
    plot_a <- function() {
        put(1.5, 3.4, "a", 9, TRUE)
        line(119, 3.4, 123.5, 3.4, RED, 1.5)
        put(125, 3.4, "Intervention", 6.1)
        line(149, 3.4, 153.5, 3.4, BLUE, 1.5, lty = 2)
        put(155, 3.4, "Control", 6.1)
        pops <- screening_label(c("Overall", "High risk", "Low risk"))
        yticks <- list(c(0, 0.2, 0.4, 0.6), c(0, 2.5, 5, 7.5, 10), c(0, 0.1, 0.2, 0.3, 0.4))
        put(1.6, 32.5, "CIED implantation (%)", 6, ha = "center", srt = 90)
        for (j in seq_along(pops)) {
            pop <- pops[j]
            start <- 2 + 60 * (j - 1)
            xl <- start + 8
            xr <- start + 55
            top <- 15.5
            bot <- 48.5
            yy <- yticks[[j]]
            ymax <- max(yy)
            fx <- function(t) xl + t/180 * (xr - xl)
            fy <- function(p) bot - p/ymax * (bot - top)
            title <- if (j == 1) 
                "Overall"
            else if (j == 2) 
                "Flagged by the AI-ECG\nscreening algorithm"
            else "Not flagged by the AI-ECG\nscreening algorithm"
            put((xl + xr)/2, 9, title, 6.5, TRUE, "center")
            stopifnot(max(100 * K$upper_incidence[K$population == pop], na.rm = TRUE) <= ymax)
            for (t in yy[-1]) line(xl, fy(t), xr, fy(t), col = "#E9EAEC", lwd = 0.45)
            line(fx(90), top, fx(90), bot, col = "#A3A6AB", lwd = 0.65, lty = 3)
            for (a in c("Intervention", "Control")) {
                z <- K[K$population == pop & K$arm == a, , drop = FALSE]
                z <- z[order(z$time), ]
                lower <- stepxy(z$time, 100 * z$lower_incidence)
                upper <- stepxy(z$time, 100 * z$upper_incidence)
                polygon(c(fx(lower$x), rev(fx(upper$x))), c(fy(lower$y), rev(fy(upper$y))), col = adjustcolor(COL[a], alpha.f = 0.1), border = NA)
            }
            for (a in c("Intervention", "Control")) {
                z <- K[K$population == pop & K$arm == a, , drop = FALSE]
                z <- z[order(z$time), ]
                s <- stepxy(z$time, 100 * z$cumulative_incidence)
                lines(fx(s$x), fy(s$y), col = COL[a], lwd = 1.25, lty = if (a == "Control") 
                  2
                else 1)
            }
            line(xl, top, xl, bot, lwd = 0.7)
            line(xl, bot, xr, bot, lwd = 0.7)
            for (t in yy) {
                line(xl - 0.6, fy(t), xl, fy(t))
                lab <- if (j == 2) 
                  format(t, trim = TRUE)
                else sprintf("%.1f", t)
                put(xl - 1.1, fy(t), lab, 5.5, ha = "right")
            }
            for (t in seq(0, 180, 30)) {
                line(fx(t), bot, fx(t), bot + 0.7)
                put(fx(t), 51.4, t, 5.4, ha = "center")
            }
            put(xl, 60.7, "No. at risk", 5.6)
            for (a in c("Intervention", "Control")) {
                y <- if (a == "Intervention") 
                  64.9
                else 68.3
                put(start + 3.8, y, if (a == "Intervention") 
                  "I"
                else "C", 5.7, TRUE, ha = "center", color = COL[a])
                for (t in seq(0, 180, 30)) {
                  z <- NR[NR$population == pop & NR$arm == a & NR$day == t, , drop = FALSE]
                  stopifnot(nrow(z) == 1)
                  put(fx(t), y, z$n_risk, 5.2, ha = "center", color = COL[a])
                }
            }
        }
        put(90, 56.1, "Days from randomization", 6, ha = "center")
    }
    plot_b <- function() {
        put(1.5, 76.8, "b", 9, TRUE)
        put(9, 76.8, "Follow-up / population", 5.9, TRUE)
        put(X$inter, 76.8, "Intervention", 5.7, TRUE, "center")
        put(X$control, 76.8, "Control", 5.7, TRUE, "center")
        put(X$or, 76.8, "Odds ratio (95% CI)", 5.6, TRUE, "center")
        put(X$p, 76.8, expression(bolditalic(P) ~ bold("value")), 5.3, TRUE, "center")
        put(X$interaction, 76.8, expression(bolditalic(P) ~ bold("for interaction")), 5, TRUE, "center")
        line(2, 80.1, 178.5, 80.1, lwd = 0.75)
        line(tx(1), 80.5, tx(1), 179.9, col = "#AAADB3", lwd = 0.85, lty = 3)
        y <- 83.2
        for (day in seq(30, 180, 30)) {
            rect(2, y - PITCH * 0.43, 178.5, y + PITCH * 0.43, col = if (day == 90) 
                "#E5EAF1"
            else GRAY, border = NA)
            put(3.4, y, if (day == 90) 
                "90 days"
            else paste(day, "days"), 6, TRUE)
            it <- IT[IT$day == day, , drop = FALSE]
            stopifnot(nrow(it) == 1)
            put(X$interaction, y, pstr(it$p, it$method), 5.6, ha = "center")
            for (pop in screening_label(c("Overall", "High risk", "Low risk"))) {
                y <- y + PITCH
                r <- E[E$day == day & E$population_label == pop, , drop = FALSE]
                stopifnot(nrow(r) == 1)
                put(X$label + 1.6, y, pop, 5.6)
                a <- armstr(r, "1")
                b <- armstr(r, "0")
                stopifnot(X$label + 1.6 + strwidth(pop, cex = 5.6/12) < X$inter - strwidth(a, cex = 5.8/12)/2 - 0.5)
                put(X$inter, y, a, ha = "center")
                put(X$control, y, b, ha = "center")
                if (is.finite(r$OR)) {
                  line(tx(max(r$lower, 0.25)), y, tx(min(r$upper, 32)), y, col = BLUE, lwd = 1.2)
                  if (r$lower < 0.25) 
                    arrow(X$left + 2.4, y, X$left - 0.6, y)
                  if (r$upper > 32) 
                    arrow(X$right - 2.4, y, X$right + 0.6, y)
                  if (r$OR >= 0.25 && r$OR <= 32) 
                    points(tx(r$OR), y, pch = 22, cex = 0.48, bg = BLUE, col = BLUE, lwd = 0.65)
                  else arrow(tx(min(max(r$OR, 0.25), 32)), y, if (r$OR < 0.25) 
                    X$left - 0.8
                  else X$right + 0.8, y)
                  ci <- sprintf("%.2f (%.2f–%.2f)", r$OR, r$lower, r$upper)
                }
                else ci <- "NE"
                p <- pstr(r$p, r$method)
                put(X$or, y, ci, 5.5, ha = "center")
                put(X$p, y, p, 5.6, ha = "center")
            }
            y <- y + PITCH + GAP
        }
        bottom <- 181.9
        line(X$left, bottom, X$right, bottom)
        for (t in c(0.25, 1, 4, 32)) {
            line(tx(t), bottom, tx(t), bottom + 1.5)
            put(tx(t), bottom + 3.3, t, 5.2, ha = "center")
        }
        put((X$left + X$right)/2, bottom + 7.2, "Intervention / control odds ratio", 5.3, ha = "center")
    }
    draw <- function() {
        par(mar = rep(0, 4), xaxs = "i", yaxs = "i", family = "sans", xpd = NA)
        plot.new()
        plot.window(xlim = c(0, W), ylim = c(H, 0))
        plot_a()
        plot_b()
    }
    open_figure_device(file.path(out, "Supplemental_Fig_13.pdf"), W, H)
    draw()
    dev.off()
} 

precise_sss_timing <- function (root) 
{
    suppressPackageStartupMessages({
        library(survival)
    })
    dir.create(file.path(root, "analysis"), showWarnings = FALSE, recursive = TRUE)
    precise_core_helpers()
    precise_screening_labels()
    d <- read.csv(file.path(root, "data/PRECISE-AI.csv"), check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
    stopifnot(!anyNA(d$group), !anyNA(d$screening_group), all(d$group %in% c("Intervention", "Control")), all(d$screening_group %in% c("High risk", "Low risk")), all(d$ai_ecg_screening_flag %in% 
        0:1), all(d$screening_group == ifelse(d$ai_ecg_screening_flag == 1, "High risk", "Low risk")), all(is.na(d$cied_indication_days) | (is.finite(d$cied_indication_days) & 
        d$cied_indication_days >= 0)))
    d <- prepare_analysis_factors(d)
    stopifnot(!anyNA(d$sick_sinus_syndrome), all(d$sick_sinus_syndrome %in% 0:1), !any(d$sick_sinus_syndrome == 1 & is.na(d$cied_indication_days)))
    d$ANALYSIS_EVENT_DAYS <- ifelse(d$sick_sinus_syndrome == 1, d$cied_indication_days, NA_real_)
    d$TIME <- ifelse(is.na(d$ANALYSIS_EVENT_DAYS), 180, pmin(d$ANALYSIS_EVENT_DAYS, 180))
    d$EVENT <- as.integer(!is.na(d$ANALYSIS_EVENT_DAYS) & d$ANALYSIS_EVENT_DAYS <= 180)
    cutoffs <- seq(30, 180, by = 30)
    effects <- interactions <- list()
    for (t in cutoffs) {
        d$EVENT_BY_TIME <- as.integer(!is.na(d$ANALYSIS_EVENT_DAYS) & d$ANALYSIS_EVENT_DAYS <= t)
        br0 <- block_rule(d, "EVENT_BY_TIME")
        br1 <- block_rule(d, "EVENT_BY_TIME", "screening_group")
        i <- interaction(d, "EVENT_BY_TIME", "screening_group", "ITT", use_firth = br1$firth, block_id = paste0("day", t, "_risk"), trigger_levels = br1$trigger_levels)
        i$day <- t
        interactions[[length(interactions) + 1L]] <- i
        for (pop in c("Overall", "High risk", "Low risk")) {
            z <- if (pop == "Overall") 
                d
            else d[d$screening_group == pop, , drop = FALSE]
            br <- if (pop == "Overall") 
                br0
            else br1
            a <- effect(z, "EVENT_BY_TIME", "ITT", sg = if (pop == "Overall") 
                "Overall"
            else "screening_group", lv = if (pop == "Overall") 
                "All"
            else pop, use_firth = br$firth, block_id = paste0("day", t, if (pop == "Overall") 
                "_overall"
            else "_risk"), trigger_levels = br$trigger_levels)
            a$day <- t
            a$population_label <- screening_label(pop)
            a$outcome_label <- paste("Sick sinus syndrome by day", t)
            a$interaction_p <- i$p
            a$interaction_method <- i$method
            effects[[length(effects) + 1L]] <- a
        }
    }
    effects <- do.call(rbind, effects)
    interactions <- do.call(rbind, interactions)
    effects <- effects[c("day", "population_label", setdiff(names(effects), c("day", "population_label")))]
    interactions <- interactions[c("day", setdiff(names(interactions), "day"))]
    write.csv(effects, file.path(root, "analysis/Supplemental_Fig_14_effects.csv"), row.names = FALSE, na = "")
    write.csv(interactions, file.path(root, "analysis/Supplemental_Fig_14_interactions.csv"), row.names = FALSE, na = "")
    km <- risk <- list()
    for (pop in c("Overall", "High risk", "Low risk")) {
        zz <- if (pop == "Overall") 
            d
        else d[d$screening_group == pop, , drop = FALSE]
        pl <- screening_label(pop)
        for (arm in c("Intervention", "Control")) {
            z <- zz[zz$group == arm, , drop = FALSE]
            f <- survfit(Surv(TIME, EVENT) ~ 1, data = z, conf.int = 0.95, conf.type = "log-log")
            n <- nrow(z)
            baseline <- data.frame(population = pl, arm = arm, time = 0, survival = 1, lower_survival = 1, upper_survival = 1, cumulative_incidence = 0, 
                lower_incidence = 0, upper_incidence = 0, n_risk = n, n_event = 0, n_censor = 0)
            coord <- data.frame(population = pl, arm = arm, time = f$time, survival = f$surv, lower_survival = f$lower, upper_survival = f$upper, 
                cumulative_incidence = 1 - f$surv, lower_incidence = 1 - f$upper, upper_incidence = 1 - f$lower, n_risk = f$n.risk, n_event = f$n.event, 
                n_censor = f$n.censor)
            stopifnot(all(is.finite(coord$survival)), all(is.finite(coord$lower_survival)), all(is.finite(coord$upper_survival)))
            if (tail(coord$time, 1L) < 180) {
                last <- tail(coord, 1L)
                last$time <- 180
                last$n_risk <- sum(z$TIME >= 180)
                last$n_event <- 0
                last$n_censor <- sum(z$TIME == 180 & z$EVENT == 0)
                coord <- rbind(coord, last)
            }
            km[[length(km) + 1L]] <- rbind(baseline, coord)
            for (t in c(0, cutoffs)) {
                risk[[length(risk) + 1L]] <- data.frame(population = pl, arm = arm, day = t, n_risk = sum(z$TIME >= t), cumulative_events = sum(z$EVENT == 
                  1 & z$TIME <= t), cumulative_censors = sum(z$EVENT == 0 & z$TIME <= t), n_total = n)
            }
        }
    }
    km <- do.call(rbind, km)
    risk <- do.call(rbind, risk)
    write.csv(km, file.path(root, "analysis/Supplemental_Fig_14_km_coordinates.csv"), row.names = FALSE, na = "")
    write.csv(risk, file.path(root, "analysis/Supplemental_Fig_14_km_at_risk.csv"), row.names = FALSE, na = "")
    stopifnot(nrow(effects) == 18L, nrow(interactions) == 6L, all(diff(effects$event1[effects$population_label == "Overall"]) >= 0), all(diff(effects$event0[effects$population_label == 
        "Overall"]) >= 0), all(effects$status %in% c("OK", "No events in either arm")), 
        all(interactions$status %in% c("OK", "No events across subgroup block", "No events in either arm within a subgroup", "Ordinary interaction NE: nonfinite MLE in treatment-by-subgroup outcome cells")))
    for (pop in unique(effects$population_label)) {
        a <- effects[effects$population_label == pop, ]
        stopifnot(length(unique(a$n1)) == 1L, length(unique(a$n0)) == 1L, all(diff(a$event1) >= 0), all(diff(a$event0) >= 0))
    }
} 

precise_sss_timing_plot <- function (root) 
{
    precise_screening_labels()
    E <- read.csv(file.path(root, "analysis/Supplemental_Fig_14_effects.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    IT <- read.csv(file.path(root, "analysis/Supplemental_Fig_14_interactions.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    K <- read.csv(file.path(root, "analysis/Supplemental_Fig_14_km_coordinates.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    NR <- read.csv(file.path(root, "analysis/Supplemental_Fig_14_km_at_risk.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    E$population_label <- screening_label(E$population_label)
    K$population <- screening_label(K$population)
    NR$population <- screening_label(NR$population)
    out <- file.path(root, "figures")
    dir.create(out, showWarnings = FALSE, recursive = TRUE)
    H <- 194
    W <- 180
    PITCH <- 3.8
    GAP <- 0.8
    BLUE <- "#4B638F"
    RED <- "#C75D47"
    INK <- "#303135"
    RULE <- "#7B7D80"
    GRAY <- "#F1F2F4"
    COL <- c(Intervention = RED, Control = BLUE)
    num <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
    put <- function(x, y, s, size = 5.8, bold = FALSE, ha = "left", color = INK, srt = 0) {
        text(x, y, labels = s, cex = size/12, font = if (bold) 
            2
        else 1, adj = c(switch(ha, left = 0, center = 0.5, right = 1), 0.5), col = color, srt = srt)
    }
    line <- function(x1, y1, x2, y2, col = RULE, lwd = 0.65, lty = 1) segments(x1, y1, x2, y2, col = col, lwd = lwd, lty = lty)
    arrow <- function(x1, y1, x2, y2, col = BLUE, lwd = 1) arrows(x1, y1, x2, y2, length = 0.035, angle = 25, code = 2, col = col, lwd = lwd)
    pstr <- function(p, method = "") if (!is.finite(p)) 
        "NE"
    else paste0(if (p < 0.001) 
        "<0.001"
    else sprintf("%.3f", p), if (grepl("Firth", method)) 
        "†"
    else "")
    armstr <- function(r, k) sprintf("%s/%s (%.1f%%)", num(r[[paste0("event", k)]]), num(r[[paste0("n", k)]]), 100 * r[[paste0("event", k)]]/r[[paste0("n", 
        k)]])
    stepxy <- function(x, y) list(x = c(x[1], rep(x[-1], each = 2)), y = c(rep(y[-length(y)], each = 2), tail(y, 1)))
    X <- list(label = 3.4, inter = 60, control = 86, left = 102, right = 129, or = 144.6, p = 159.5, interaction = 171.7)
    tx <- function(v) X$left + (log(v) - log(0.25))/(log(32) - log(0.25)) * (X$right - X$left)
    plot_a <- function() {
        put(1.5, 3.4, "a", 9, TRUE)
        line(119, 3.4, 123.5, 3.4, RED, 1.5)
        put(125, 3.4, "Intervention", 6.1)
        line(149, 3.4, 153.5, 3.4, BLUE, 1.5, lty = 2)
        put(155, 3.4, "Control", 6.1)
        pops <- screening_label(c("Overall", "High risk", "Low risk"))
        yticks <- list(c(0, 0.2, 0.4, 0.6), c(0, 2.5, 5, 7.5, 10), c(0, 0.1, 0.2, 0.3))
        put(1.6, 32.5, "Sick sinus syndrome (%)", 6, ha = "center", srt = 90)
        for (j in seq_along(pops)) {
            pop <- pops[j]
            start <- 2 + 60 * (j - 1)
            xl <- start + 8
            xr <- start + 55
            top <- 15.5
            bot <- 48.5
            yy <- yticks[[j]]
            ymax <- max(yy)
            fx <- function(t) xl + t/180 * (xr - xl)
            fy <- function(p) bot - p/ymax * (bot - top)
            title <- if (j == 1) 
                "Overall"
            else if (j == 2) 
                "Flagged by the AI-ECG\nscreening algorithm"
            else "Not flagged by the AI-ECG\nscreening algorithm"
            put((xl + xr)/2, 9, title, 6.5, TRUE, "center")
            stopifnot(max(100 * K$upper_incidence[K$population == pop], na.rm = TRUE) <= ymax)
            for (t in yy[-1]) line(xl, fy(t), xr, fy(t), col = "#E9EAEC", lwd = 0.45)
            line(fx(90), top, fx(90), bot, col = "#A3A6AB", lwd = 0.65, lty = 3)
            for (a in c("Intervention", "Control")) {
                z <- K[K$population == pop & K$arm == a, , drop = FALSE]
                z <- z[order(z$time), ]
                lower <- stepxy(z$time, 100 * z$lower_incidence)
                upper <- stepxy(z$time, 100 * z$upper_incidence)
                polygon(c(fx(lower$x), rev(fx(upper$x))), c(fy(lower$y), rev(fy(upper$y))), col = adjustcolor(COL[a], alpha.f = 0.1), border = NA)
            }
            for (a in c("Intervention", "Control")) {
                z <- K[K$population == pop & K$arm == a, , drop = FALSE]
                z <- z[order(z$time), ]
                s <- stepxy(z$time, 100 * z$cumulative_incidence)
                lines(fx(s$x), fy(s$y), col = COL[a], lwd = 1.25, lty = if (a == "Control") 
                  2
                else 1)
            }
            line(xl, top, xl, bot, lwd = 0.7)
            line(xl, bot, xr, bot, lwd = 0.7)
            for (t in yy) {
                line(xl - 0.6, fy(t), xl, fy(t))
                lab <- if (j == 2) 
                  format(t, trim = TRUE)
                else sprintf("%.1f", t)
                put(xl - 1.1, fy(t), lab, 5.5, ha = "right")
            }
            for (t in seq(0, 180, 30)) {
                line(fx(t), bot, fx(t), bot + 0.7)
                put(fx(t), 51.4, t, 5.4, ha = "center")
            }
            put(xl, 60.7, "No. at risk", 5.6)
            for (a in c("Intervention", "Control")) {
                y <- if (a == "Intervention") 
                  64.9
                else 68.3
                put(start + 3.8, y, if (a == "Intervention") 
                  "I"
                else "C", 5.7, TRUE, ha = "center", color = COL[a])
                for (t in seq(0, 180, 30)) {
                  z <- NR[NR$population == pop & NR$arm == a & NR$day == t, , drop = FALSE]
                  stopifnot(nrow(z) == 1)
                  put(fx(t), y, z$n_risk, 5.2, ha = "center", color = COL[a])
                }
            }
        }
        put(90, 56.1, "Days from randomization", 6, ha = "center")
    }
    plot_b <- function() {
        put(1.5, 76.8, "b", 9, TRUE)
        put(9, 76.8, "Follow-up / population", 5.9, TRUE)
        put(X$inter, 76.8, "Intervention", 5.7, TRUE, "center")
        put(X$control, 76.8, "Control", 5.7, TRUE, "center")
        put(X$or, 76.8, "Odds ratio (95% CI)", 5.6, TRUE, "center")
        put(X$p, 76.8, expression(bolditalic(P) ~ bold("value")), 5.3, TRUE, "center")
        put(X$interaction, 76.8, expression(bolditalic(P) ~ bold("for interaction")), 5, TRUE, "center")
        line(2, 80.1, 178.5, 80.1, lwd = 0.75)
        line(tx(1), 80.5, tx(1), 179.9, col = "#AAADB3", lwd = 0.85, lty = 3)
        y <- 83.2
        for (day in seq(30, 180, 30)) {
            rect(2, y - PITCH * 0.43, 178.5, y + PITCH * 0.43, col = if (day == 90) 
                "#E5EAF1"
            else GRAY, border = NA)
            put(3.4, y, if (day == 90) 
                "90 days"
            else paste(day, "days"), 6, TRUE)
            it <- IT[IT$day == day, , drop = FALSE]
            stopifnot(nrow(it) == 1)
            put(X$interaction, y, pstr(it$p, it$method), 5.6, ha = "center")
            for (pop in screening_label(c("Overall", "High risk", "Low risk"))) {
                y <- y + PITCH
                r <- E[E$day == day & E$population_label == pop, , drop = FALSE]
                stopifnot(nrow(r) == 1)
                put(X$label + 1.6, y, pop, 5.6)
                a <- armstr(r, "1")
                b <- armstr(r, "0")
                stopifnot(X$label + 1.6 + strwidth(pop, cex = 5.6/12) < X$inter - strwidth(a, cex = 5.8/12)/2 - 0.5)
                put(X$inter, y, a, ha = "center")
                put(X$control, y, b, ha = "center")
                if (is.finite(r$OR)) {
                  line(tx(max(r$lower, 0.25)), y, tx(min(r$upper, 32)), y, col = BLUE, lwd = 1.2)
                  if (r$lower < 0.25) 
                    arrow(X$left + 2.4, y, X$left - 0.6, y)
                  if (r$upper > 32) 
                    arrow(X$right - 2.4, y, X$right + 0.6, y)
                  if (r$OR >= 0.25 && r$OR <= 32) 
                    points(tx(r$OR), y, pch = 22, cex = 0.48, bg = BLUE, col = BLUE, lwd = 0.65)
                  else arrow(tx(min(max(r$OR, 0.25), 32)), y, if (r$OR < 0.25) 
                    X$left - 0.8
                  else X$right + 0.8, y)
                  ci <- sprintf("%.2f (%.2f–%.2f)", r$OR, r$lower, r$upper)
                }
                else ci <- "NE"
                p <- pstr(r$p, r$method)
                put(X$or, y, ci, 5.5, ha = "center")
                put(X$p, y, p, 5.6, ha = "center")
            }
            y <- y + PITCH + GAP
        }
        bottom <- 181.9
        line(X$left, bottom, X$right, bottom)
        for (t in c(0.25, 1, 4, 32)) {
            line(tx(t), bottom, tx(t), bottom + 1.5)
            put(tx(t), bottom + 3.3, t, 5.2, ha = "center")
        }
        put((X$left + X$right)/2, bottom + 7.2, "Intervention / control odds ratio", 5.3, ha = "center")
    }
    draw <- function() {
        par(mar = rep(0, 4), xaxs = "i", yaxs = "i", family = "sans", xpd = NA)
        plot.new()
        plot.window(xlim = c(0, W), ylim = c(H, 0))
        plot_a()
        plot_b()
    }
    open_figure_device(file.path(out, "Supplemental_Fig_14.pdf"), W, H)
    draw()
    dev.off()
} 

precise_time_distribution <- function (root) 
{
    dir.create(file.path(root, "figures"), recursive = TRUE, showWarnings = FALSE)
    precise_exact_test()
    d <- read.csv(file.path(root, "data/PRECISE-AI.csv"), check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
    required <- c("group", "cied_implantation", "cied_indication", "cied_implantation_days", "cied_indication_days")
    stopifnot(all(required %in% names(d)), !anyNA(d$group), all(d$group %in% c("Intervention", "Control")), !anyNA(d$cied_implantation), all(d$cied_implantation %in% 
        0:1))
    arms <- c("Intervention", "Control")
    specs <- list(a = list(title = "With CIED indication", subset = d$cied_indication == 1, values = d$cied_indication_days), 
        b = list(title = "With CIED implantation", subset = d$cied_implantation == 1, values = d$cied_indication_days), 
        c = list(title = "Without CIED implantation", subset = d$cied_implantation == 0, values = d$cied_indication_days), 
        d = list(title = "Indication to implantation", subset = d$cied_implantation == 1, values = d$cied_implantation_days - d$cied_indication_days))
    S <- tests <- list()
    samples <- list()
    for (panel in names(specs)) {
        s <- specs[[panel]]
        samples[[panel]] <- list()
        for (arm in arms) {
            z <- s$values[s$subset & d$group == arm]
            x <- z[!is.na(z)]
            stopifnot(length(x) > 0, all(is.finite(x)), all(x >= 0 & x <= 180))
            samples[[panel]][[arm]] <- x
            q <- quantile(x, c(0, 0.25, 0.5, 0.75, 1), names = FALSE, type = 7)
            iqr <- q[4] - q[2]
            lo <- min(x[x >= q[2] - 1.5 * iqr])
            hi <- max(x[x <= q[4] + 1.5 * iqr])
            S[[length(S) + 1L]] <- data.frame(panel = panel, arm = arm, observed_n = length(x),
                q1 = q[2], median = q[3], q3 = q[4], whisker_lower = lo, whisker_upper = hi)
        }
        r <- exact_rank_sum(samples[[panel]][["Intervention"]], samples[[panel]][["Control"]])
        tests[[panel]] <- data.frame(panel = panel, p_raw = r$p.value)
    }
    S <- do.call(rbind, S)
    T <- do.call(rbind, tests)
    fmt <- function(x) sub("\\.?0+$", "", sprintf("%.2f", x))
    S$plot_annotation <- paste0("Median (IQR): ", fmt(S$median), " (", fmt(S$q1), "-", fmt(S$q3), ") d")
    for (arm in arms) stopifnot(identical(sort(samples$a[[arm]]), sort(c(samples$b[[arm]], samples$c[[arm]]))))
    stopifnot(all(!is.na(d$cied_implantation_days[d$cied_implantation == 1])), all(is.na(d$cied_implantation_days[d$cied_implantation == 0])))
    swarm_offsets <- function(y, diameter = 1.25) {
        stopifnot(all(is.finite(y)), diameter > 0)
        ord <- order(y, seq_along(y))
        sy <- y[ord]
        sx <- numeric(length(y))
        for (i in seq_along(sy)) {
            prev <- if (i > 1L) 
                seq_len(i - 1L)
            else integer()
            near <- prev[abs(sy[prev] - sy[i]) < diameter]
            if (!length(near)) {
                sx[i] <- 0
                next
            }
            delta <- sqrt(pmax(0, diameter^2 - (sy[near] - sy[i])^2))
            candidates <- unique(c(0, sx[near] - delta, sx[near] + delta))
            valid <- vapply(candidates, function(a) all((a - sx[near])^2 + (sy[i] - sy[near])^2 >= diameter^2 - 1e-08), logical(1))
            candidates <- candidates[valid]
            stopifnot(length(candidates) > 0)
            preferred_side <- if (sum(sx[prev]) > 0) 
                -1
            else 1
            sx[i] <- candidates[order(abs(candidates), -preferred_side * candidates)][1L]
        }
        ans <- numeric(length(y))
        ans[ord] <- sx
        if (length(y) > 1L) {
            dist2 <- outer(ans, ans, "-")^2 + outer(y, y, "-")^2
            diag(dist2) <- Inf
            stopifnot(min(dist2) >= diameter^2 - 1e-08)
        }
        ans
    }
    W <- 180
    H <- 149
    RED <- "#C75D47"
    BLUE <- "#4B638F"
    INK <- "#303135"
    RULE <- "#7B7D80"
    COL <- c(Intervention = RED, Control = BLUE)
    put <- function(x, y, s, size = 6.4, bold = FALSE, ha = "left", color = INK, srt = 0) {
        text(x, y, labels = s, cex = size/12, font = if (bold) 
            2
        else 1, adj = c(switch(ha, left = 0, center = 0.5, right = 1), 0.5), col = color, srt = srt)
    }
    line <- function(x1, y1, x2, y2, col = RULE, lwd = 0.65) segments(x1, y1, x2, y2, col = col, lwd = lwd)
    coords <- list()
    for (i in seq_along(specs)) {
        panel <- names(specs)[i]
        y0 <- 5 + 72 * ((i - 1L)%/%2L)
        for (j in seq_along(arms)) {
            x <- sort(samples[[panel]][[arms[j]]])
            py <- y0 + 62 - x/180 * 37
            dx <- swarm_offsets(py)
            stopifnot(max(abs(dx)) < 12)
            coords[[length(coords) + 1L]] <- data.frame(panel = panel, arm = arms[j], horizontal_offset_mm = dx, vertical_coordinate_mm = py)
        }
    }
    C <- do.call(rbind, coords)
    draw <- function() {
        par(mar = rep(0, 4), xaxs = "i", yaxs = "i", family = "sans", xpd = NA)
        plot.new()
        plot.window(xlim = c(0, W), ylim = c(H, 0))
        for (i in seq_along(specs)) {
            panel <- names(specs)[i]
            s <- specs[[panel]]
            x0 <- 2 + 90 * ((i - 1L)%%2L)
            y0 <- 5 + 72 * ((i - 1L)%/%2L)
            xl <- x0 + 11
            xr <- x0 + 84
            top <- y0 + 25
            bot <- y0 + 62
            fy <- function(v) bot - v/180 * (bot - top)
            put(x0, y0, panel, 9, TRUE)
            put(xl, y0, s$title, 7.5, TRUE)
            ss <- S[S$panel == panel, ]
            tt <- T[T$panel == panel, ]
            for (j in seq_along(arms)) {
                arm <- arms[j]
                cx <- x0 + c(29, 67)[j]
                z <- ss[ss$arm == arm, ]
                put(cx, y0 + 5.4, paste0(arm, ", n = ", z$observed_n), 6.3, ha = "center", color = COL[arm])
                stopifnot(strwidth(z$plot_annotation, cex = 6.1/12, units = "user") < 37)
                put(cx, y0 + 10.3, z$plot_annotation, 6.1, ha = "center", color = COL[arm])
            }
            ptext <- bquote(italic(P) == .(sprintf("%.3f", tt$p_raw)))
            put((xl + xr)/2, y0 + 17, ptext, 6.6, ha = "center")
            for (v in seq(0, 180, 30)) {
                if (v > 0) 
                  line(xl, fy(v), xr, fy(v), col = "#ECEDEE", lwd = 0.45)
                line(xl - 0.7, fy(v), xl, fy(v))
                put(xl - 1.5, fy(v), v, 6, ha = "right")
            }
            line(xl, top, xl, bot)
            line(xl, bot, xr, bot)
            for (j in seq_along(arms)) {
                arm <- arms[j]
                cx <- x0 + c(29, 67)[j]
                z <- ss[ss$arm == arm, ]
                line(cx, fy(z$whisker_lower), cx, fy(z$whisker_upper), col = COL[arm], lwd = 0.8)
                line(cx - 2, fy(z$whisker_lower), cx + 2, fy(z$whisker_lower), col = COL[arm], lwd = 0.8)
                line(cx - 2, fy(z$whisker_upper), cx + 2, fy(z$whisker_upper), col = COL[arm], lwd = 0.8)
                rect(cx - 3.6, fy(z$q3), cx + 3.6, fy(z$q1), col = adjustcolor(COL[arm], alpha.f = 0.12), border = COL[arm], lwd = 0.8)
                zz <- C[C$panel == panel & C$arm == arm, ]
                symbols(cx + zz$horizontal_offset_mm, zz$vertical_coordinate_mm, circles = rep(0.55, nrow(zz)), inches = FALSE, add = TRUE, fg = "white", 
                  bg = adjustcolor(COL[arm], alpha.f = 0.9), lwd = 0.3)
                line(cx - 3.6, fy(z$median), cx + 3.6, fy(z$median), col = INK, lwd = 1.1)
                put(cx, bot + 4.7, arm, 6.6, ha = "center", color = COL[arm])
            }
            yl <- if (panel == "d") 
                "Indication to implantation (days)"
            else "Time to CIED indication (days)"
            put(x0 + 1.3, (top + bot)/2, yl, 6.3, ha = "center", srt = 90)
        }
    }
    open_figure_device(file.path(root, "figures", "Supplemental_Fig_15.pdf"), W, H)
    draw()
    dev.off()
} 

