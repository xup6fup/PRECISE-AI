draw_precise_ecg <- function(d, path) {
  arms <- c('Intervention', 'Control')
  cols <- c(Intervention = '#485E8F', Control = '#E85543')
  ink <- '#363435'
  stopifnot(nrow(d) > 0L, all(d$group %in% arms),
            !anyNA(d$ai_ecg_screening_flag), all(d$ai_ecg_screening_flag %in% 0:1),
            !anyNA(d$high_risk_ecg_count_1y))
  ecg_lists <- lapply(d$ecg_days_before_index, function(s) {
    if (is.na(s) || !nzchar(s)) stop('Missing ECG list.')
    as.numeric(strsplit(s, ';', fixed = TRUE)[[1]])
  })
  stopifnot(identical(as.numeric(lengths(ecg_lists)), as.numeric(d$ecg_count_1y)))
  all_lags <- unlist(ecg_lists, use.names = FALSE)
  stopifnot(all(is.finite(all_lags)), all(all_lags == floor(all_lags)),
            all(all_lags >= 0 & all_lags <= 365.5))
  panel_titles <- c(a = 'Most recent ECG', b = 'All available ECGs',
    c = 'Flagged by the AI-ECG\nscreening algorithm',
    d = 'Not flagged by the AI-ECG screening algorithm\nwith a high-risk ECG in the preceding year')
  xlabels <- c(a = 'Days before randomization: most recent ECG',
               b = 'Days before randomization: all available ECGs',
               c = 'Days before randomization: most recent high-risk ECG',
               d = 'Days before randomization: most recent high-risk ECG')
  breaks <- c(0, seq(16, 376, 15))
  raw_values <- list(); summary_rows <- list(); count_rows <- list()
  for (p in names(panel_titles)) for (arm in arms) {
    mask <- d$group == arm
    if (p == 'c') mask <- mask & d$ai_ecg_screening_flag == 1
    if (p == 'd') mask <- mask & d$ai_ecg_screening_flag == 0 & d$high_risk_ecg_count_1y > 0
    vals <- switch(p, a = d$last_ecg_days_before_index[mask],
                   b = unlist(ecg_lists[mask], use.names = FALSE),
                   c = d$last_high_risk_ecg_days_before_index[mask],
                   d = d$last_high_risk_ecg_days_before_index[mask])
    stopifnot(length(vals) > 0L, all(is.finite(vals)),
              all(vals == floor(vals)), all(vals >= 0 & vals <= 365.5))
    summary_rows[[paste(p, arm)]] <- data.frame(Panel = p, Arm = arm, N = length(vals))
    counts <- hist(vals, breaks = breaks, plot = FALSE, right = FALSE)$counts
    count_rows[[paste(p, arm)]] <- data.frame(Panel = p, Arm = arm,
      Bin_lower_inclusive = head(breaks, -1), Bin_upper_inclusive = tail(breaks, -1) - 1,
      Count = counts, Percent = 100 * counts / length(vals))
    raw_values[[paste(p, arm)]] <- vals
  }
  summaries <- do.call(rbind, summary_rows); rownames(summaries) <- NULL
  counts <- do.call(rbind, count_rows); rownames(counts) <- NULL
  inset_rows <- list()
  for (p in c('a', 'b', 'c')) for (arm in arms) {
    vals <- raw_values[[paste(p, arm)]]
    lo <- seq(0, 14, 2)
    nn <- vapply(lo, function(x) sum(vals >= x & vals < x + 2), integer(1))
    inset_rows[[paste(p, arm)]] <- data.frame(Panel = p, Arm = arm,
      Bin_lower_inclusive = lo, Bin_upper_inclusive = lo + 1,
      Count = nn, Percent = nn / length(vals) * 100)
    stopifnot(sum(nn) == counts$Count[counts$Panel == p & counts$Arm == arm &
                                      counts$Bin_lower_inclusive == 0])
  }
  inset_counts <- do.call(rbind, inset_rows)
  rownames(inset_counts) <- NULL
  draw_panel <- function(p, region, page_mm) {
    grid::pushViewport(grid::viewport(x = region[1], y = region[3],
      width = region[2] - region[1], height = region[4] - region[3],
      just = c('left', 'bottom')))
    on.exit(grid::popViewport())
    pw <- page_mm[1] * diff(region[1:2]); ph <- page_mm[2] * diff(region[3:4])
    left <- 12.0
    right <- 2.8
    bottom <- 19
    top <- 16
    x0 <- left / pw; x1 <- 1 - right / pw
    y0 <- bottom / ph; y1 <- 1 - top / ph
    size_axis <- 6.4
    size_title <- 7.0
    title <- panel_titles[p]
    grid::grid.text(title, x = (x0 + x1) / 2, y = 1 - 3 / ph, just = 'top',
      gp = grid::gpar(fontfamily = 'Arial', fontsize = size_title, fontface = 'bold',
                col = ink, lineheight = 1.12))
    grid::grid.text(p, x = 1.1 / pw, y = 1 - 2.7 / ph,
      just = c('left', 'top'), gp = grid::gpar(fontfamily = 'Arial', fontsize = 10,
                                      fontface = 'bold', col = ink))
    q <- counts[counts$Panel == p, ]
    ymax <- max(q$Count) * 1.23
    ticks <- pretty(c(0, ymax), n = 4)
    ticks <- ticks[ticks >= 0 & ticks <= ymax]
    if (p %in% c('c', 'd')) {
      ymax <- 100
      ticks <- seq(0, 100, 20)
    }
    tx <- function(x) x0 + x / 378 * (x1 - x0)
    ty <- function(y) y0 + y / ymax * (y1 - y0)
    has_inset <- p %in% c('a', 'b', 'c')
    if (has_inset) {
      inset_bottom <- if (p == 'b') .50 else .345
      ix0 <- x0 + .35 * (x1 - x0); ix1 <- x0 + .97 * (x1 - x0)
      iy0 <- y0 + inset_bottom * (y1 - y0); iy1 <- y0 + .86 * (y1 - y0)
      first_max <- max(q$Count[q$Bin_lower_inclusive == 0])
      high <- min(first_max + ymax * .13, ymax * .98)
      anchor <- first_max * .90
      grid::grid.rect(x = tx(0), y = y0, width = tx(16) - tx(0), height = ty(high) - y0,
        just = c('left', 'bottom'), gp = grid::gpar(fill = '#F0F1F3', col = NA))
      grid::grid.polygon(x = c(tx(16), tx(16), ix0, ix0),
        y = c(ty(high), ty(anchor), iy0, iy1), gp = grid::gpar(fill = '#F0F1F3', col = NA))
      grid::grid.segments(x0 = rep(tx(16), 2), y0 = c(ty(high), ty(anchor)),
        x1 = rep(ix0, 2), y1 = c(iy1, iy0), gp = grid::gpar(col = '#A9AFB8', lwd = .55))
    }
    grid::grid.lines(x = c(x0, x0, x1), y = c(y1, y0, y0),
               gp = grid::gpar(col = ink, lwd = .75))
    for (y in ticks) {
      grid::grid.lines(x = c(x0 - 1.1 / pw, x0), y = rep(ty(y), 2), gp = grid::gpar(col = ink, lwd = .7))
      grid::grid.text(format(y, big.mark = ',', trim = TRUE, scientific = FALSE),
        x = x0 - 1.7 / pw, y = ty(y), just = 'right',
        gp = grid::gpar(fontfamily = 'Arial', fontsize = size_axis, col = ink))
    }
    xticks <- seq(0, 360, 60)
    for (x in xticks) {
      grid::grid.lines(x = rep(tx(x), 2), y = c(y0, y0 - 1.1 / ph), gp = grid::gpar(col = ink, lwd = .7))
      grid::grid.text(x, x = tx(x), y = y0 - 2.3 / ph, just = 'top',
                gp = grid::gpar(fontfamily = 'Arial', fontsize = size_axis, col = ink))
    }
    for (arm in arms) {
      z <- q[q$Arm == arm, ]
      lo <- z$Bin_lower_inclusive; hi <- z$Bin_upper_inclusive + 1
      mid <- (lo + hi) / 2
      l <- if (arm == 'Intervention') lo + .65 else mid
      r <- if (arm == 'Intervention') mid else hi - .65
      grid::grid.rect(x = tx(l), y = y0, width = tx(r) - tx(l), height = ty(z$Count) - y0,
        just = c('left', 'bottom'), gp = grid::gpar(fill = cols[arm], col = NA))
      nz <- z$Count > 0
      grid::grid.text(sprintf('%.1f%%', z$Percent[nz]), x = tx((l[nz] + r[nz]) / 2),
        y = ty(z$Count[nz]) + 1 / ph, rot = 90, just = c('left', 'centre'),
        gp = grid::gpar(fontfamily = 'Arial', fontsize = 4.5,
                  col = cols[arm]))
    }
    if (has_inset) {
      z <- inset_counts[inset_counts$Panel == p, ]
      zymax <- max(z$Count, 1) * 1.29
      zticks <- pretty(c(0, zymax), n = 3)
      zticks <- zticks[zticks >= 0 & zticks <= zymax]
      zx <- function(x) ix0 + (x + .2) / 16.4 * (ix1 - ix0)
      zy <- function(y) iy0 + y / zymax * (iy1 - iy0)
      grid::grid.text('0-15 days', x = (ix0 + ix1) / 2,
        y = y0 + .935 * (y1 - y0), gp = grid::gpar(fontfamily = 'Arial',
          fontsize = 6.2, col = ink))
      grid::grid.lines(x = c(ix0, ix0, ix1), y = c(iy1, iy0, iy0),
        gp = grid::gpar(col = ink, lwd = .6))
      for (arm in arms) {
        zz <- z[z$Arm == arm, ]
        ll <- zz$Bin_lower_inclusive + if (arm == 'Intervention') .18 else 1
        rr <- ll + .82
        grid::grid.rect(x = zx(ll), y = iy0, width = zx(rr) - zx(ll),
          height = zy(zz$Count) - iy0, just = c('left', 'bottom'),
          gp = grid::gpar(fill = cols[arm], col = NA))
        grid::grid.text(sprintf('%.1f%%', zz$Percent), x = zx((ll + rr) / 2),
          y = zy(zz$Count) + .025 * (iy1 - iy0), rot = 90,
          just = c('left', 'centre'), gp = grid::gpar(fontfamily = 'Arial',
            fontsize = 4.5, col = cols[arm]))
      }
      for (x in c(0, 5, 10, 15)) {
        grid::grid.lines(x = rep(zx(x), 2), y = c(iy0, iy0 - .7 / ph), gp = grid::gpar(col = ink, lwd = .6))
        grid::grid.text(x, x = zx(x), y = iy0 - .035 * (y1 - y0), just = 'top',
          gp = grid::gpar(fontfamily = 'Arial', fontsize = 5.3, col = ink))
      }
      for (y in zticks) {
        grid::grid.lines(x = c(ix0 - .7 / pw, ix0), y = rep(zy(y), 2), gp = grid::gpar(col = ink, lwd = .6))
        grid::grid.text(y, x = ix0 - .012 * (x1 - x0), y = zy(y), just = 'right',
          gp = grid::gpar(fontfamily = 'Arial', fontsize = 5.3, col = ink))
      }
    }
    grid::grid.text(if (p == 'b') 'ECGs (n)' else 'Participants (n)',
      x = 2.0 / pw, y = (y0 + y1) / 2, rot = 90,
      gp = grid::gpar(fontfamily = 'Arial', fontsize = size_axis + .5, col = ink))
    grid::grid.text(xlabels[p], x = (x0 + x1) / 2, y = y0 - 7.1 / ph, just = 'top',
      gp = grid::gpar(fontfamily = 'Arial', fontsize = 6.1, col = ink))
    s <- summaries[summaries$Panel == p, ]
    legend_size <- 6.0
    legend_y <- 3.4 / ph
    legend_x <- c(.15, .61)
    for (i in seq_along(arms)) {
      n <- s$N[s$Arm == arms[i]]
      label <- paste0(arms[i], ' (n = ', format(n, big.mark = ',', trim = TRUE), ')')
      grid::grid.rect(x = legend_x[i], y = legend_y, width = grid::unit(1.6, 'mm'),
        height = grid::unit(1.6, 'mm'), gp = grid::gpar(fill = cols[i], col = NA))
      grid::grid.text(label, x = legend_x[i] + 1.8 / pw,
        y = legend_y, just = 'left', gp = grid::gpar(fontfamily = 'Arial', fontsize = legend_size, col = ink))
    }
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(path, width = 180 / 25.4, height = 160 / 25.4,
                      family = 'Arial', bg = 'white')
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  regions <- list(a = c(0, .5, .5, 1), b = c(.5, 1, .5, 1),
                  c = c(0, .5, 0, .5), d = c(.5, 1, 0, .5))
  for (p in names(regions)) draw_panel(p, regions[[p]], c(180, 160))
  invisible(path)
}
