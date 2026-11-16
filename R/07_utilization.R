calculate_utilization <- function(d) {
  fields <- c('group', 'ai_ecg_screening_flag', 'high_risk_ecg_count_1y',
              'ecg_count_180d', 'recorded_holter_count_180d', 'selected_matched_cohort')
  stopifnot(all(fields %in% names(d)), !anyNA(d[fields]),
            all(d$group %in% c('Intervention', 'Control')),
            all(d$ai_ecg_screening_flag %in% 0:1), all(d$selected_matched_cohort %in% 0:1),
            all(d$ecg_count_180d >= 0), all(d$recorded_holter_count_180d >= 0),
            all(d$ecg_count_180d == floor(d$ecg_count_180d)),
            all(d$recorded_holter_count_180d == floor(d$recorded_holter_count_180d)))
  strata <- ifelse(d$ai_ecg_screening_flag == 1, 1L,
                   ifelse(d$high_risk_ecg_count_1y > 0, 2L, 3L))
  labels <- c('Overall', 'Flagged by the AI-ECG screening algorithm',
    'Not flagged by the AI-ECG screening algorithm (Protocol defined AI-ECG high risk)',
    'Not flagged by the AI-ECG screening algorithm (Protocol defined AI-ECG low risk)')
  metrics <- c(ECG = 'ecg_count_180d', Routine_Holter = 'recorded_holter_count_180d',
               Protocol_Holter = 'selected_matched_cohort')
  res <- do.call(rbind, lapply(names(metrics), function(m) {
    do.call(rbind, lapply(0:3, function(s) {
      selected <- if (s == 0) rep(TRUE, nrow(d)) else strata == s
      a <- d[[metrics[[m]]]][selected & d$group == 'Intervention']
      b <- d[[metrics[[m]]]][selected & d$group == 'Control']
      stopifnot(length(a) > 1, length(b) > 1)
      variance_a <- stats::var(a)
      variance_b <- stats::var(b)
      p_mean <- if (variance_a == 0 && variance_b == 0) {
        if (mean(a) == mean(b)) 1 else NA_real_
      } else stats::t.test(a, b, alternative = 'two.sided', var.equal = FALSE)$p.value
      tab <- rbind(c(sum(a > 0), sum(a == 0)), c(sum(b > 0), sum(b == 0)))
      p_proportion <- if (any(colSums(tab) == 0)) 1 else
        stats::fisher.test(tab, alternative = 'two.sided')$p.value
      data.frame(Metric = m, Stratum = s, Group = labels[s + 1L],
        Intervention_N = length(a), Control_N = length(b),
        Intervention_Total = sum(a), Control_Total = sum(b),
        Intervention_Mean = mean(a), Intervention_SD = stats::sd(a),
        Control_Mean = mean(b), Control_SD = stats::sd(b),
        P_Mean = p_mean, Intervention_With_Exam = sum(a > 0),
        Control_With_Exam = sum(b > 0),
        Intervention_Percent = 100 * mean(a > 0),
        Control_Percent = 100 * mean(b > 0), P_Proportion = p_proportion,
        stringsAsFactors = FALSE)
    }))
  }))
  rownames(res) <- NULL
  res
}

draw_precise_utilization <- function(d, path) {
  res <- calculate_utilization(d)
  metrics <- c('ECG', 'Routine_Holter', 'Protocol_Holter')
  titles <- c('12-lead ECG', 'Usual-care Holter', 'Trial-arranged ambulatory\nECG monitoring')
  cols <- c(Intervention = '#485E8F', Control = '#E85543')
  ink <- '#363435'
  height <- 114
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(path, width = 180 / 25.4, height = height / 25.4,
                      family = 'Arial', bg = 'white')
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  tx <- function(x) grid::unit(x, 'mm')
  ty <- function(y) grid::unit(height - y, 'mm')
  txt <- function(label, x, y, size = 6.0, just = 'centre', bold = FALSE, col = ink) {
    grid::grid.text(label, x = tx(x), y = ty(y), just = just,
      gp = grid::gpar(fontfamily = 'Arial', fontsize = size,
                     fontface = if (bold) 'bold' else 'plain', col = col,
                     lineheight = 1.02))
  }
  line <- function(x, y, col = ink, lwd = .45) {
    grid::grid.lines(x = tx(x), y = ty(y), gp = grid::gpar(col = col, lwd = lwd))
  }
  ptext <- function(p) {
    ifelse(is.na(p), 'NE', ifelse(p < .001, '<0.001', sprintf('%.3f', p)))
  }
  dots <- function(a, b, x, y, upper) {
    xx <- x + c(a, b) / upper * 13
    line(xx, c(y - .65, y + .65), col = '#B6B6B6', lwd = .55)
    grid::grid.points(x = tx(xx), y = ty(c(y - .65, y + .65)), pch = c(16, 17),
      size = tx(c(1.35, 1.5)), gp = grid::gpar(col = unname(cols)))
  }
  axis <- function(x, y, upper, percent = FALSE) {
    line(c(x, x + 13), rep(y, 2), col = '#858585')
    line(rep(x, 2), c(y, y + .6), col = '#858585')
    line(rep(x + 13, 2), c(y, y + .6), col = '#858585')
    txt('0', x, y - 1.7, 6.0)
    txt(paste0(format(upper, trim = TRUE), if (percent) '%' else ''),
        x + 13, y - 1.7, 6.0)
  }
  txt('Examinations within 180 days', 3, 4, 7.7, just = 'left', bold = TRUE)
  txt('Mean (SD) per participant', 77, 4, 6.3, bold = TRUE)
  txt('Participants with \u22651 examination, n/N (%)', 145, 4, 6.3, bold = TRUE)
  txt('Group', 3, 9, 6.0, just = 'left', bold = TRUE)
  txt('Intervention', 67, 9, 6.0, bold = TRUE, col = cols['Intervention'])
  txt('Control', 85, 9, 6.0, bold = TRUE, col = cols['Control'])
  txt(expression(italic(P)), 99, 9, 6.0)
  txt('Intervention', 135, 9, 6.0, bold = TRUE, col = cols['Intervention'])
  txt('Control', 158, 9, 6.0, bold = TRUE, col = cols['Control'])
  txt(expression(italic(P)), 175, 9, 6.0)
  line(c(3, 179), rep(11, 2))
  row_labels <- c('Overall', 'Flagged', 'Not flagged, protocol high risk',
                  'Not flagged, protocol low risk')
  heads <- c(15, 44, 77)
  for (j in seq_along(metrics)) {
    z <- res[res$Metric == metrics[j], ]
    mean_upper <- max(pretty(c(0, z$Intervention_Mean, z$Control_Mean), n = 3))
    pct_upper <- max(pretty(c(0, z$Intervention_Percent, z$Control_Percent), n = 3))
    if (mean_upper <= 0) mean_upper <- 1
    if (pct_upper <= 0) pct_upper <- 1
    txt(paste0(letters[j], '  ', titles[j]), 3, heads[j], 6.5, just = 'left', bold = TRUE)
    axis(44, heads[j] + 1.2, mean_upper)
    axis(111, heads[j] + 1.2, pct_upper, TRUE)
    for (k in seq_len(nrow(z))) {
      r <- z[k, ]
      y <- heads[j] + 6 + (k - 1) * 6
      if (k %% 2 == 0) {
        grid::grid.rect(x = tx(3), y = ty(y), width = tx(176), height = tx(5.6),
          just = c('left', 'centre'), gp = grid::gpar(fill = '#F7F7F7', col = NA))
      }
      txt(row_labels[k], 3, y, 6.0, just = 'left')
      dots(r$Intervention_Mean, r$Control_Mean, 44, y, mean_upper)
      txt(sprintf('%.2f (%.2f)', r$Intervention_Mean, r$Intervention_SD),
          67, y, 6.0, col = cols['Intervention'])
      txt(sprintf('%.2f (%.2f)', r$Control_Mean, r$Control_SD),
          85, y, 6.0, col = cols['Control'])
      txt(ptext(r$P_Mean), 99, y, 6.0)
      dots(r$Intervention_Percent, r$Control_Percent, 111, y, pct_upper)
      txt(sprintf('%d/%d (%.1f)', r$Intervention_With_Exam,
                  r$Intervention_N, r$Intervention_Percent),
          135, y, 6.0, col = cols['Intervention'])
      txt(sprintf('%d/%d (%.1f)', r$Control_With_Exam,
                  r$Control_N, r$Control_Percent),
          158, y, 6.0, col = cols['Control'])
      txt(ptext(r$P_Proportion), 175, y, 6.0)
    }
    line(c(3, 179), rep(heads[j] + 27, 2), col = '#BFBFBF', lwd = .35)
  }
  txt('Flagged = Flagged by the AI-ECG screening algorithm; not flagged = Not flagged by the AI-ECG screening algorithm.',
      3, 107, 6.0, just = 'left')
  txt('Protocol high risk / low risk = Protocol defined AI-ECG high risk / Protocol defined AI-ECG low risk.',
      3, 111, 6.0, just = 'left')
  invisible(res)
}
