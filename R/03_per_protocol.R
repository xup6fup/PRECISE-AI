precise_figure2 <- function (root) 
{
    precise_screening_labels()
    portable_pdf <- function(...) {
        if (capabilities("cairo")) 
            grDevices::cairo_pdf(...)
        else grDevices::pdf(...)
    }
    for (z in c("analysis", "figures")) dir.create(file.path(root, z), recursive = TRUE, showWarnings = FALSE)
    core <- new.env(parent = environment())
    precise_core_helpers(core)
    D <- read.csv(file.path(root, "data/PRECISE-AI.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    stopifnot(nrow(D) == 11492L, all(D$screening_group == ifelse(D$ai_ecg_screening_flag == 1, "High risk", "Low risk")))
    D <- core$prepare_analysis_factors(D)
    vars <- c("cied_implantation", precise_secondary_outcomes, "cied_adverse_event")
    core$outlabs <- c(core$outlabs, cied_adverse_event = "CIED-related complications")
    el <- il <- list()
    for (v in vars) {
        r0 <- core$block_rule(D, v)
        id0 <- paste("ITT", v, "Overall", sep = ":")
        el[[length(el) + 1L]] <- core$effect(D, v, "Overall", use_firth = r0$firth, block_id = id0, trigger_levels = r0$trigger_levels)
        r1 <- core$block_rule(D, v, "screening_group")
        id1 <- paste("ITT", v, "screening_group", sep = ":")
        it <- core$interaction(D, v, "screening_group", "Overall", use_firth = r1$firth, block_id = id1, trigger_levels = r1$trigger_levels)
        il[[length(il) + 1L]] <- it
        for (lev in levels(D$screening_group)) {
            ef <- core$effect(D[D$screening_group == lev, , drop = FALSE], v, lev, sg = "screening_group", lv = lev, use_firth = r1$firth, block_id = id1, trigger_levels = r1$trigger_levels)
            el[[length(el) + 1L]] <- ef
        }
    }
    E <- do.call(rbind, el)
    I <- do.call(rbind, il)
    stopifnot(nrow(E) == 18L, nrow(I) == 6L, all(E$status %in% c("OK", "No events in either arm")), all(I$status %in% c("OK", "No events across subgroup block", "No events in either arm within a subgroup")))
    out <- file.path(root, "figures")
    BLUE <- "#4B638F"
    GRAY <- "#F1F2F4"
    INK <- "#303135"
    RULE <- "#7B7D80"
    PITCH <- 3.8
    GAP <- 0.8
    outcomes <- c(cied_implantation = "CIED implantation", high_grade_av_block = "High-grade atrioventricular block", complete_av_block = "Complete atrioventricular block", 
        sick_sinus_syndrome = "Sick sinus syndrome", ventricular_arrhythmia = "Ventricular arrhythmia")
    safety_outcomes <- c(cied_adverse_event = "CIED-related complications")
    figure2_outcomes <- c(outcomes, safety_outcomes)
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
    arrow <- function(x1, y1, x2, y2, col = BLUE, lwd = 0.7) arrows(x1, y1, x2, y2, length = 0.035, angle = 25, code = 2, col = col, lwd = lwd)
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
            fmt <- function(v) if (v > 0 && v < 0.001) 
                sprintf("%.4f", v)
            else if (v > 0 && v < 0.01) 
                sprintf("%.3f", v)
            else sprintf("%.2f", v)
            ci <- sprintf("%s (%s-%s)", fmt(r$OR), fmt(r$lower), fmt(r$upper))
        }
        else ci <- "NE"
        pp <- pstr(r$p, r$method)
        put(x$or, y, ci, 5.6, ha = "center")
        put(x$p, y, pp, 5.6, ha = "center")
    }
    draw_figure2 <- function() {
        canvas(118)
        header(2.4)
        y <- 10.8
        forest_axis(7.8, 106.9)
        for (j in seq_along(figure2_outcomes)) {
            oc <- names(figure2_outcomes)[j]
            prefix <- if (oc == "cied_implantation") 
                "Primary outcome: "
            else if (oc %in% names(safety_outcomes)) 
                "Safety outcome: "
            else "Secondary outcome: "
            ttl <- paste0(prefix, figure2_outcomes[j])
            band(y, ttl)
            it <- interaction("Overall", oc, "screening_group")
            put(X$interaction, y, pstr(it$p, it$method), 5.6, ha = "center")
            for (lev in c("All", "High risk", "Low risk")) {
                y <- y + PITCH
                r <- pick(population = if (lev == "All") 
                  "Overall"
                else lev, adjustment = "Crude", outcome = oc, subgroup_var = if (lev == "All") 
                  "Overall"
                else "screening_group", subgroup_level = lev)
                label <- screening_label(lev)
                put(X$label + 1.6, y, label, 5.6)
                stopifnot(X$label + 1.6 + strwidth(label, cex = 5.6/12) < X$inter - strwidth(arm(r, "1"), cex = 5.8/12)/2 - 0.5)
                estimate(y, r)
            }
            y <- y + PITCH + GAP
        }
    }
    f <- file.path(out, "Fig_02.pdf")
    portable_pdf(f, width = 180/25.4, height = 118/25.4, family = "sans", pointsize = 12)
    draw_figure2()
    dev.off()
} 

precise_figure3_analysis <- function (root) 
{
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
        X <- fullx[, keep, drop = FALSE]
        fitted <- safe(glm.fit(x = X, y = pp$treatment, family = binomial(), control = glm.control(maxit = 1000, epsilon = 1e-12)))
        ps_fit <- fitted$value
        ps_warn <- fitted$warnings
        ps_ok <- !is.null(ps_fit) && isTRUE(ps_fit$converged) && !isTRUE(ps_fit$boundary) && all(is.finite(ps_fit$coefficients)) && all(is.finite(ps_fit$fitted.values)) && 
            all(ps_fit$fitted.values > 0 & ps_fit$fitted.values < 1) && !grepl("numerically 0 or 1|did not converge|separat", ps_warn, ignore.case = TRUE)
        if (ps_ok) {
            pp$PS <- as.numeric(ps_fit$fitted.values)
            pp$IPTW <- ifelse(pp$treatment == 1, 1/pp$PS, 1/(1 - pp$PS))
            pp$FIRTH_WEIGHT <- nrow(pp) * pp$IPTW/sum(pp$IPTW)
            ps_ok <- all(is.finite(pp$IPTW) & pp$IPTW > 0) & all(is.finite(pp$FIRTH_WEIGHT) & pp$FIRTH_WEIGHT > 0)
        }
        if (!ps_ok) 
            pp$PS <- pp$IPTW <- pp$FIRTH_WEIGHT <- NA_real_
    }
    mapping <- read.csv(file.path(root, "data/Fig_03_endpoint_labels.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    outlabs[mapping$outcome] <- mapping$label
    endpoints <- c(outs, precise_safety_outcomes)
    stopifnot(ps_ok, nrow(pp) > 0L)
    make_effects <- function(z) {
        rows <- list()
        for (v in endpoints) {
            rule <- block_rule(z, v)
            for (adj in c("Crude", "IPTW")) {
                a <- effect(z, v, "Monitored", adj, use_firth = rule$firth, block_id = paste("PP", v, "Overall", sep = ":"), trigger_levels = rule$trigger_levels)
                rows[[length(rows) + 1L]] <- a
            }
        }
        ef <- do.call(rbind, rows)
        ef
    }
    effects <- make_effects(pp)
    write_csv(effects, "Fig_03_effects")
} 

precise_figure3_plot <- function (root) 
{
    portable_pdf <- function(...) {
        if (capabilities("cairo")) 
            grDevices::cairo_pdf(...)
        else grDevices::pdf(...)
    }
    out <- file.path(root, "figures")
    dir.create(out, showWarnings = FALSE, recursive = TRUE)
    E <- read.csv(file.path(root, "analysis/Fig_03_effects.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    mapping <- read.csv(file.path(root, "data/Fig_03_endpoint_labels.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    figure3_outcomes <- setNames(mapping$label, mapping$outcome)
    BLUE <- "#4B638F"
    GRAY <- "#F1F2F4"
    INK <- "#303135"
    RULE <- "#7B7D80"
    PITCH <- 3.8
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
    arrow <- function(x1, y1, x2, y2, col = BLUE, lwd = 0.7) arrows(x1, y1, x2, y2, length = 0.035, angle = 25, code = 2, col = col, lwd = lwd)
    tx <- function(v, lo, hi, x) x$left + (log(v) - log(lo))/(log(hi) - log(lo)) * (x$right - x$left)
    estimate <- function(y, r, lo, hi, x, counts = TRUE, open = FALSE) {
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
    draw_figure3 <- function() {
        canvas(122.2)
        x <- list(label = 3, inter = 72, control = 95, left = 112, right = 138, or = 154, p = 175)
        for (k in c("label", "inter", "control", "or", "p")) put(x[[k]], 4, if (k == "p") 
            expression(bolditalic(P) ~ bold("value"))
        else c(label = "Outcome", inter = "Intervention", control = "Control", or = "Odds ratio (95% CI)")[k], 6, TRUE, if (k == "label") 
            "left"
        else "center")
        line(2, 8, 178.5, 8)
        titles <- c(14.8, 31.6 + c(0, 12.2, 24.4, 36.6), 84.2, 96.4)
        ay <- tail(titles, 1) + 12.2
        line(tx(1, 0.1, 100, x), 8, tx(1, 0.1, 100, x), ay, col = "#AAADB3", lty = 3)
        band(11, "Primary outcome", size = 6.5)
        band(27.8, "Secondary outcomes", size = 6.5)
        band(80.4, "Safety outcomes", size = 6.5)
        for (j in seq_along(figure3_outcomes)) {
            oc <- names(figure3_outcomes)[j]
            y <- titles[j]
            put(5, y, figure3_outcomes[j], 6.4, TRUE)
            for (adj in c("Crude", "IPTW")) {
                weighted <- adj == "IPTW"
                yy <- y + if (weighted) 
                  2 * PITCH
                else PITCH
                r <- pick(population = "Monitored", adjustment = adj, outcome = oc, subgroup_var = "Overall", subgroup_level = "All")
                put(if (weighted) 
                  9
                else 6.3, yy, if (weighted) 
                  "IPTW adjusted"
                else "Per-protocol", 6.2)
                estimate(yy, r, 0.1, 100, x, counts = !weighted, open = weighted)
            }
        }
        line(x$left, ay, x$right, ay)
        for (v in c(0.1, 1, 10, 100)) {
            xx <- tx(v, 0.1, 100, x)
            line(xx, ay, xx, ay + 1.4)
            put(xx, ay + 3.5, as.character(v), 5.6, ha = "center")
        }
        put((x$left + x$right)/2, ay + 8, "Intervention / control odds ratio", 5.6, ha = "center")
    }
    target <- file.path(out, "Fig_03.pdf")
    portable_pdf(target, width = 180/25.4, height = 122.2/25.4, family = "sans", pointsize = 12)
    draw_figure3()
    dev.off()
} 

