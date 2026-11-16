precise_recent_analysis <- function (root) 
{
    precise_core_helpers()
    dir.create(file.path(root, "analysis"), showWarnings = FALSE, recursive = TRUE)
    d <- prepare_analysis_factors(read.csv(file.path(root, "data", "PRECISE-AI.csv"), check.names = FALSE, stringsAsFactors = FALSE, na.strings = ""))
    stopifnot(nrow(d) == 11492L, !anyDuplicated(d$row_number), all(!is.na(d$ai_ecg_screening_flag)), identical(as.character(d$screening_group), ifelse(d$ai_ecg_screening_flag == 
        1, "High risk", "Low risk")), all(vapply(d[outs], function(x) all(!is.na(x) & x %in% 0:1), logical(1))))
    recent <- is.finite(d$last_ecg_days_before_index) & d$last_ecg_days_before_index >= 0 & d$last_ecg_days_before_index <= 91.5
    z <- d[recent, , drop = FALSE]
    write_csv <- function(x, nm) write.csv(x, file.path(root, "analysis", paste0(nm, ".csv")), row.names = FALSE, na = "")
    E <- I <- list()
    for (v in outs) {
        r0 <- block_rule(z, v)
        rr <- block_rule(z, v, "screening_group")
        b0 <- paste("RecentECG", v, "Overall", sep = ":")
        br <- paste("RecentECG", v, "screening_group", sep = ":")
        it <- interaction(z, v, "screening_group", "Recent ECG cohort", use_firth = rr$firth, block_id = br, trigger_levels = rr$trigger_levels)
        I[[v]] <- it
        for (lv in c("All", "High risk", "Low risk")) {
            zz <- if (lv == "All") 
                z
            else z[z$screening_group == lv, , drop = FALSE]
            rule <- if (lv == "All") 
                r0
            else rr
            e <- effect(zz, v, if (lv == "All") 
                "Overall"
            else lv, sg = if (lv == "All") 
                "Overall"
            else "screening_group", lv = lv, use_firth = rule$firth, block_id = if (lv == "All") 
                b0
            else br, trigger_levels = rule$trigger_levels)
            E[[length(E) + 1L]] <- e
        }
    }
    E <- do.call(rbind, E)
    I <- do.call(rbind, I)
    write_csv(E, "effects")
    write_csv(I, "interactions")
    stopifnot(nrow(E) == 15L, nrow(I) == 5L, all(E$event1 <= E$n1), all(E$event0 <= E$n0),
        !any(grepl("fail|invalid", E$status, ignore.case = TRUE)), !any(grepl("fail|invalid", I$status, ignore.case = TRUE)))
} 

precise_recent_plot <- function (root) 
{
    precise_screening_labels()
    E <- read.csv(file.path(root, "analysis", "effects.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    I <- read.csv(file.path(root, "analysis", "interactions.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    out <- file.path(root, "figures")
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    BLUE <- "#4B638F"
    GRAY <- "#F1F2F4"
    INK <- "#303135"
    RULE <- "#7B7D80"
    PITCH <- 3.8
    GAP <- 0.8
    num <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
    pstr <- function(p, method = "") if (!is.finite(p)) 
        "NE"
    else paste0(if (p < 0.001) 
        "<0.001"
    else sprintf("%.3f", p), if (grepl("Firth", method)) 
        "†"
    else "")
    arm <- function(r, k) {
        n <- r[[paste0("n", k)]]
        e <- r[[paste0("event", k)]]
        if (n == 0) 
            "0/0 (NE)"
        else sprintf("%s/%s (%.1f%%)", num(e), num(n), 100 * e/n)
    }
    outcomes <- c(cied_implantation = "CIED implantation", high_grade_av_block = "High-grade atrioventricular block", complete_av_block = "Complete atrioventricular block", 
        sick_sinus_syndrome = "Sick sinus syndrome", ventricular_arrhythmia = "Ventricular arrhythmia")
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
    H <- 102
    draw <- function() {
        canvas(H)
        header(2.4)
        forest_axis(7.8, 90.9)
        y <- 10.8
        for (oc in names(outcomes)) {
            band(y, paste0(if (oc == "cied_implantation") 
                "Primary outcome: "
            else "Secondary outcome: ", outcomes[[oc]]))
            it <- I[I$outcome == oc, , drop = FALSE]
            stopifnot(nrow(it) == 1)
            put(X$interaction, y, pstr(it$p, it$method), 5.6, ha = "center")
            for (lv in c("All", "High risk", "Low risk")) {
                y <- y + PITCH
                r <- E[E$outcome == oc & E$subgroup_level == lv, , drop = FALSE]
                stopifnot(nrow(r) == 1)
                label <- screening_label(lv)
                put(X$label + 1.6, y, label, 5.6)
                stopifnot(X$label + 1.6 + strwidth(label, cex = 5.6/12) < X$inter - strwidth(arm(r, "1"), cex = 5.8/12)/2 - 0.5)
                estimate(y, r)
            }
            y <- y + PITCH + GAP
        }
    }
    target <- file.path(out, "Supplemental_Fig_12.pdf")
    if (capabilities("cairo"))
        cairo_pdf(target, width = 180/25.4, height = H/25.4, family = "sans", pointsize = 12)
    else pdf(target, width = 180/25.4, height = H/25.4, family = "sans", pointsize = 12, useDingbats = FALSE)
    draw()
    dev.off()
} 

precise_recent_table <- function (root) 
{
    suppressPackageStartupMessages(library(jsonlite))
    input <- file.path(root, "data/PRECISE-AI.csv")
    for (z in c("analysis", "reference")) dir.create(file.path(root, z), recursive = TRUE, showWarnings = FALSE)
    x <- read.csv(input, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
    variables <- precise_baseline_variables
    stopifnot(nrow(x) == 11492L, all(c(variables, "last_ecg_days_before_index") %in% names(x)))
    stopifnot(all(x$group %in% c("Intervention", "Control")), !anyNA(x[variables]), all(x$ai_ecg_screening_flag %in% 0:1))
    if (!"row_number" %in% names(x)) 
        x$row_number <- seq_len(nrow(x))
    stopifnot(!anyDuplicated(x$row_number), identical(x$row_number, seq_len(nrow(x))))
    recent <- is.finite(x$last_ecg_days_before_index) & x$last_ecg_days_before_index >= 0 & x$last_ecg_days_before_index <= 91.5
    d <- x[recent, , drop = FALSE]
    d$screening_group <- ifelse(d$ai_ecg_screening_flag == 1, "High risk", "Low risk")
    template <- fromJSON(file.path(root, "reference/table1_template.json"), simplifyVector = FALSE)
    stopifnot(length(template$rows) == 63L, all(vapply(template$rows, function(r) length(r$cells), integer(1)) == 9L))
    labels <- vapply(template$rows, function(r) r$cells[[1]], character(1))
    template_grid <- do.call(rbind, lapply(template$rows, function(r) unlist(r$cells)))
    grid <- matrix("", nrow = 63L, ncol = 9L)
    grid[, 1] <- labels
    grid[1, ] <- template_grid[1, ]
    map <- data.frame(row_index = integer(), variable = character(), level = character(), stringsAsFactors = FALSE)
    addmap <- function(rows, idx, levels = "1") {
        n <- length(rows)
        map <<- rbind(map, data.frame(row_index = rows, variable = variables[rep(idx, length.out = n)], level = rep(levels, 
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
    stopifnot(nrow(map) == 47L)
    fmtn <- function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE)
    summary_cell <- function(z, m) {
        v <- z[[m$variable]]
        n <- nrow(z)
        ans <- list()
        if (m$level == "Continuous") {
            ans$mean <- mean(v)
            ans$sd <- sd(v)
            ans$text <- sprintf("%.1f ± %.1f", ans$mean, ans$sd)
        }
        else {
            ans$n_event <- sum(if (is.character(v) || is.factor(v)) v == m$level else v == as.numeric(m$level))
            ans$percentage <- 100 * ans$n_event/n
            ans$text <- sprintf("%s (%.1f)", fmtn(ans$n_event), ans$percentage)
        }
        ans
    }
    for (pop in c("Overall", "High risk", "Low risk")) for (arm in c("Intervention", "Control")) {
        selected <- d$group == arm & (pop == "Overall" | d$screening_group == pop)
        z <- d[selected, , drop = FALSE]
        ci <- switch(pop, Overall = 1L, `High risk` = 4L, `Low risk` = 7L) + as.integer(arm == "Control")
        grid[2, ci + 1L] <- sprintf("%s\n(N = %s)", arm, fmtn(nrow(z)))
        for (i in seq_len(nrow(map))) {
            m <- map[i, ]
            if (m$row_index %in% c(4, 5) && pop != "Overall") 
                next
            a <- summary_cell(z, m)
            grid[m$row_index + 1L, ci + 1L] <- a$text
        }
    }
    base <- file.path(root, "analysis/Supplemental_Table_09_Recent_ECG_Baseline")
    write.table(grid, paste0(base, ".csv"), sep = ",", row.names = FALSE, col.names = FALSE, na = "", qmethod = "double", fileEncoding = "UTF-8")
} 

