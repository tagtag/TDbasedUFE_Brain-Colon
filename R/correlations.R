compute_correlations <- function(scores, out) {
  stopifnot(nrow(scores) == 75L, !anyDuplicated(scores$sample_id))
  scores$age <- factor(scores$age, levels = ages)
  design <- model.matrix(~ age, scores)
  qr_age <- qr(design)
  n <- nrow(scores)
  df <- n - qr_age$rank - 1L # 75 - 5 - 1 = 69, NOT the ordinary Pearson df 73.
  ordinary <- adjusted <- list()
  k <- 0L
  for (layer in names(layer_files)) for (b in behavior_columns) {
    x <- scores[[layer]]; y <- scores[[b]]
    stopifnot(all(is.finite(x)), all(is.finite(y)))
    k <- k + 1L
    ct <- cor.test(x, y, method = "pearson", alternative = "two.sided")
    ordinary[[k]] <- data.frame(molecular_score = layer, behavior = b, n = n,
      r = unname(ct$estimate), p_value = ct$p.value)
    r <- cor(qr.resid(qr_age, x), qr.resid(qr_age, y))
    statistic <- r * sqrt(df / (1-r^2))
    adjusted[[k]] <- data.frame(molecular_score = layer, behavior = b, n = n,
      age_adjustment = "categorical_5_groups", df = df, partial_r = r,
      p_value = 2 * pt(-abs(statistic), df = df))
  }
  ordinary <- do.call(rbind, ordinary); adjusted <- do.call(rbind, adjusted)
  ordinary$BH_adjusted_p <- p.adjust(ordinary$p_value, "BH")
  adjusted$BH_adjusted_p <- p.adjust(adjusted$p_value, "BH")
  ordinary <- ordinary[order(ordinary$p_value), ]
  adjusted <- adjusted[order(adjusted$p_value), ]
  csv_write(ordinary, file.path(out, "correlations_all_36.csv"))
  csv_write(adjusted, file.path(out, "age_adjusted_correlations_all_36.csv"))
  csv_write(ordinary[ordinary$BH_adjusted_p < 0.05, ], file.path(out, "Table3_correlations.csv"))
  pdf(file.path(out, "cognition_scatterplots.pdf"), width = 10, height = 5)
  par(mfrow = c(1,2))
  cols <- c("#0072B2", "#009E73", "#E69F00", "#D55E00", "#CC79A7")
  plot(scores$colon_methylation, scores$MedianRankCombinedScore,
    col = cols[scores$age], pch = 19, xlab = "Colonic methylation component score",
    ylab = "Overall cognition score", main = "Across age groups")
  legend("topleft", legend = ages, col = cols, pch = 19, bty = "n")
  plot(qr.resid(qr_age, scores$colon_methylation),
    qr.resid(qr_age, scores$MedianRankCombinedScore), col = cols[scores$age], pch = 19,
    xlab = "Age-adjusted molecular score", ylab = "Age-adjusted cognition score",
    main = "After categorical age adjustment")
  dev.off()
  invisible(list(ordinary = ordinary, adjusted = adjusted))
}
