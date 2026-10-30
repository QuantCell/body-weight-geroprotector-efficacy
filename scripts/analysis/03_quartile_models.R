# Matched weight strata, survival increases and compound comparisons.


fit_matched_quartiles <- function() {
  suppressPackageStartupMessages(library(survival))

  args <- commandArgs(trailingOnly = TRUE)
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
  package <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
  input <- if (length(args) >= 1) args[1] else file.path(package, "data", "derived", "predosing_comparison_mice.csv")
  outdir <- if (length(args) >= 2) args[2] else file.path(package, "data", "derived")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  d <- read.csv(input, stringsAsFactors = FALSE)
  d$Tx <- as.integer(d$Tx)
  d$dead <- as.integer(d$dead)
  d$matched_quartile <- as.integer(d$matched_quartile)
  d$q_factor <- factor(d$matched_quartile, levels = 1:4, labels = paste0("Q", 1:4))
  d$experiment <- factor(d$experiment)
  d$stratum <- factor(d$stratum)
  d$mouse <- factor(d$mouse)

  weighted_median <- function(z) {
    fit <- survfit(
      Surv(entry, lifespan_days, dead) ~ 1,
      data = z,
      weights = standardization_weight,
      conf.type = "none"
    )
    as.numeric(summary(fit)$table["median"])
  }

  quartile_rows <- lapply(1:4, function(q) {
    z <- d[d$matched_quartile == q, ]
    control <- z[z$Tx == 0, ]
    treated <- z[z$Tx == 1, ]
    mc <- weighted_median(control)
    mt <- weighted_median(treated)
    fit <- coxph(
      Surv(entry, lifespan_days, dead) ~ Tx + strata(stratum) + cluster(mouse),
      data = z,
      ties = "efron"
    )
    s <- summary(fit)$coefficients["Tx", ]
    data.frame(
      quartile = q,
      comparison_rows = nrow(z),
      unique_mice = length(unique(z$mouse)),
      mean_control_weight = weighted.mean(control$body_weight, control$standardization_weight),
      standardized_control_median = mc,
      standardized_treated_median = mt,
      gain_days = mt - mc,
      gain_pct = 100 * (mt - mc) / mc,
      treatment_HR = exp(s["coef"]),
      treatment_lo = exp(s["coef"] - 1.96 * s["robust se"]),
      treatment_hi = exp(s["coef"] + 1.96 * s["robust se"]),
      treatment_p = s["Pr(>|z|)"]
    )
  })
  quartiles <- do.call(rbind, quartile_rows)

  # Shared ordinal treatment-by-quartile trend while allowing each randomized
  # experiment to retain its own treatment main effect. Quartiles were assigned within
  # experiment x site x treatment, so every comparison contributes across the gradient.
  d$q_centered <- d$matched_quartile - 1
  d$tx_q <- d$Tx * d$q_centered
  trend_fit <- coxph(
    Surv(entry, lifespan_days, dead) ~ 0 + experiment:Tx + q_factor + tx_q +
      strata(stratum) + cluster(mouse),
    data = d,
    ties = "efron",
    singular.ok = TRUE
  )
  trend_s <- summary(trend_fit)$coefficients["tx_q", ]
  trend <- data.frame(
    analysis = "within-experiment-site matched quartile trend",
    experiments = length(unique(d$experiment)),
    comparison_rows = nrow(d),
    unique_mice = length(unique(d$mouse)),
    HR_per_quartile = exp(trend_s["coef"]),
    lo = exp(trend_s["coef"] - 1.96 * trend_s["robust se"]),
    hi = exp(trend_s["coef"] + 1.96 * trend_s["robust se"]),
    p = trend_s["Pr(>|z|)"]
  )

  write.csv(quartiles, file.path(outdir, "matched_quartile_standardized_survival.csv"), row.names = FALSE)
  write.csv(trend, file.path(outdir, "matched_quartile_trend.csv"), row.names = FALSE)
  print(quartiles, row.names = FALSE, digits = 5)
  print(trend, row.names = FALSE, digits = 5)

  # Export the same weighted Kaplan-Meier estimates used for the reported medians.
  curve_rows <- list()
  for (q in 1:4) {
    for (tx in 0:1) {
      z <- d[d$matched_quartile == q & d$Tx == tx, ]
      fit <- survfit(Surv(entry, lifespan_days, dead) ~ 1, data = z,
                     weights = standardization_weight, conf.type = "none")
      curve_rows[[length(curve_rows) + 1L]] <- data.frame(
        quartile = q, Tx = tx, entry = min(z$entry),
        time = fit$time, survival = fit$surv)
    }
  }
  write.csv(do.call(rbind, curve_rows),
            file.path(outdir, "matched_quartile_survival_curves.csv"), row.names = FALSE)

  # Rapamycin monotherapy is a compound subset of the pooled display.
  raw_group <- sub("^[^|]+\\|", "", d$experiment)
  rapamycin <- d[grepl("^[Rr]apa(_|$)", raw_group), ]
  stopifnot(length(unique(rapamycin$experiment)) == 6L)
  rapamycin_curves <- list()
  for (q in c(1L, 4L)) {
    for (tx in 0:1) {
      z <- rapamycin[rapamycin$matched_quartile == q & rapamycin$Tx == tx, ]
      fit <- survfit(Surv(entry, lifespan_days, dead) ~ 1, data = z,
                     weights = standardization_weight, conf.type = "none")
      rapamycin_curves[[length(rapamycin_curves) + 1L]] <- data.frame(
        quartile = q, Tx = tx, entry = min(z$entry),
        time = fit$time, survival = fit$surv)
    }
  }
  write.csv(do.call(rbind, rapamycin_curves),
            file.path(outdir, "rapamycin_quartile_survival_curves.csv"), row.names = FALSE)

}


fit_quartile_robustness <- function() {
  suppressPackageStartupMessages(library(survival))
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  root <- normalizePath(file.path(dirname(script), "../.."), winslash = "/")
  out <- file.path(root, "data/derived")
  d <- read.csv(file.path(out, "predosing_comparison_mice.csv"))
  # The common age horizon is fixed independently of treatment contrasts.
  tau <- 1000
  reps <- 2000L
  set.seed(20261005)
  stopifnot(length(unique(d$entry)) == 1L, all(d$lifespan_days > d$entry))
  entry <- d$entry[1]
  common <- d
  cuts <- lapply(split(d[d$Tx == 0, ], d$stratum[d$Tx == 0]), function(z)
    as.numeric(quantile(z$body_weight, c(.25, .5, .75), type = 7)))
  common$matched_quartile <- vapply(seq_len(nrow(common)), function(i)
    1L + sum(common$body_weight[i] > cuts[[common$stratum[i]]]), integer(1))
  # Keep a constant set of comparison cells in all four quarters, in both groups.
  counts <- table(common$stratum, common$Tx, common$matched_quartile)
  eligible <- dimnames(counts)[[1]][apply(counts > 0, 1, all)]
  common <- common[common$stratum %in% eligible, ]
  stopifnot(length(eligible) > 0L)
  write.csv(do.call(rbind, lapply(names(cuts), function(s) data.frame(
    stratum = s, q25_g = cuts[[s]][1], q50_g = cuts[[s]][2],
    q75_g = cuts[[s]][3], retained = s %in% eligible))),
    file.path(out, "common_weight_thresholds.csv"), row.names = FALSE)

  prepare <- function(z) {
    z$cell <- as.integer(factor(interaction(z$stratum, z$Tx, z$matched_quartile)))
    z$mouse_index <- match(z$mouse, unique(d$mouse))
    z
  }
  datasets <- list(separate_ranks = prepare(d), common_control_thresholds = prepare(common))
  weights <- function(z, multiplier) {
    raw <- multiplier[z$mouse_index]
    # Numeric cell order preserves the row-to-cell mapping.
    total <- rowsum(raw, z$cell, reorder = TRUE)[, 1]
    raw / total[match(z$cell, as.integer(names(total)))]
  }
  rmst <- function(z, w) {
    t <- sort(unique(z$lifespan_days))
    idx <- match(z$lifespan_days, t)
    at_time <- as.numeric(rowsum(w, idx, reorder = TRUE))
    events <- as.numeric(rowsum(w * z$dead, idx, reorder = TRUE))
    risk <- rev(cumsum(rev(at_time)))
    after <- cumprod(pmax(0, 1 - events / risk))
    keep <- t < tau
    sum(diff(c(entry, t[keep], tau)) * c(1, after[keep]))
  }
  estimate <- function(z, multiplier) {
    w <- weights(z, multiplier)
    unlist(lapply(1:4, function(q) sapply(0:1, function(tx) {
      i <- z$matched_quartile == q & z$Tx == tx
      rmst(z[i, ], w[i])
    })))
  }
  point <- lapply(datasets, estimate, rep(1, length(unique(d$mouse))))
  boots <- lapply(datasets, function(z) matrix(NA_real_, reps, 8))
  for (b in seq_len(reps)) {
    # One positive multiplier per unique mouse, reused in every comparison row.
    multiplier <- rexp(length(unique(d$mouse)))
    for (name in names(datasets)) boots[[name]][b, ] <- estimate(datasets[[name]], multiplier)
  }
  ci <- function(x) as.numeric(quantile(x, c(.025, .975)))
  quarter_rows <- list(); contrast_rows <- list(); trend_rows <- list()
  for (name in names(datasets)) {
    z <- datasets[[name]]; p <- point[[name]]; b <- boots[[name]]
    w <- weights(z, rep(1, length(unique(d$mouse))))
    for (q in 1:4) {
      cidx <- 2*q-1; tidx <- 2*q
      gain_ci <- ci(b[, tidx] - b[, cidx])
      i <- z$matched_quartile == q & z$Tx == 0
      quarter_rows[[length(quarter_rows)+1]] <- data.frame(
        grouping = name, quartile = q, horizon_age_days = tau, entry_age_days = entry,
        experiments = length(unique(z$experiment)), sites_by_experiment = length(unique(z$stratum)),
        unique_mice = length(unique(z$mouse)),
        quarter_mice = length(unique(z$mouse[z$matched_quartile == q])),
        mean_control_weight_g = weighted.mean(z$body_weight[i], w[i]),
        control_rmst_days = p[cidx], treated_rmst_days = p[tidx],
        gain_days = p[tidx]-p[cidx], gain_lo = gain_ci[1], gain_hi = gain_ci[2])
    }
    delta <- (p[8]-p[7])-(p[2]-p[1])
    delta_ci <- ci((b[,8]-b[,7])-(b[,2]-b[,1]))
    contrast_rows[[length(contrast_rows)+1]] <- data.frame(
      grouping = name, horizon_age_days = tau, bootstrap_replicates = reps,
      heaviest_minus_lightest_gain_days = delta, lo = delta_ci[1], hi = delta_ci[2])
    if (name == "common_control_thresholds") {
      z$experiment <- factor(z$experiment); z$stratum <- factor(z$stratum)
      z$q_factor <- factor(z$matched_quartile)
      z$tx_q <- z$Tx * (z$matched_quartile-1)
      fit <- coxph(Surv(entry, lifespan_days, dead) ~ 0 + experiment:Tx + q_factor +
        tx_q + strata(stratum) + cluster(mouse), data = z, ties = "efron")
      s <- summary(fit)$coefficients["tx_q", ]
      trend_rows[[1]] <- data.frame(analysis = "control-derived cut-offs in both groups", grouping = name, experiments = length(unique(z$experiment)),
        sites_by_experiment = length(unique(z$stratum)), unique_mice = length(unique(z$mouse)),
        comparison_rows = nrow(z), HR_per_quartile = exp(s["coef"]),
        lo = exp(s["coef"]-1.96*s["robust se"]),
        hi = exp(s["coef"]+1.96*s["robust se"]), p = s["Pr(>|z|)"])
    }
    # Check the integration against survival's weighted Kaplan-Meier implementation.
    i <- z$matched_quartile == 1 & z$Tx == 0
    fit <- survfit(Surv(lifespan_days-entry, dead) ~ 1, data = z[i, ], weights = w[i])
    reference <- summary(fit, rmean = tau-entry)$table["rmean"]
    stopifnot(isTRUE(all.equal(unname(reference), p[1], tolerance = 1e-8)))
  }
  write.csv(do.call(rbind, quarter_rows), file.path(out, "quartile_rmst.csv"), row.names = FALSE)
  write.csv(do.call(rbind, contrast_rows), file.path(out, "quartile_rmst_contrasts.csv"), row.names = FALSE)
  write.csv(do.call(rbind, trend_rows), file.path(out, "common_threshold_quartile_trend.csv"), row.names = FALSE)
  print(do.call(rbind, quarter_rows), row.names = FALSE)
  print(do.call(rbind, contrast_rows), row.names = FALSE)
  print(do.call(rbind, trend_rows), row.names = FALSE)

}


summarize_quartile_sensitivities <- function() {
  suppressPackageStartupMessages(library(survival))
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
  package <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
  outdir <- file.path(package, "data", "derived")

  km_median <- function(z) {
    fit <- survfit(Surv(entry, lifespan_days, dead) ~ 1, data = z, conf.type = "none")
    as.numeric(summary(fit)$table["median"])
  }

  empirical_bins <- function(v, k) {
    cut(v, breaks = quantile(v, probs = 0:k / k, type = 7),
        include.lowest = TRUE, labels = FALSE)
  }

  # Conventional strata use weight cut-offs within site, pooled across experiments and groups.
  summarize_bins <- function(d, analysis, k) {
    d$bin <- ave(d$body_weight, d$site, FUN = function(v) empirical_bins(v, k))
    rows <- lapply(seq_len(k), function(q) {
      z <- d[d$bin == q, ]; c <- z[z$Tx == 0, ]; t <- z[z$Tx == 1, ]
      mc <- km_median(c); mt <- km_median(t)
      data.frame(analysis = analysis, bins = k, bin = q,
                 mean_control_weight = mean(c$body_weight),
                 control_median = mc, treated_median = mt,
                 gain_pct = 100 * (mt - mc) / mc)
    })
    do.call(rbind, rows)
  }

  # Match the primary analysis: ranks within experiment x site x treatment group,
  # with equal cell contributions in each stratum and shared controls retained.
  standardized_bins <- function(d, analysis, k) {
    cell <- interaction(d$stratum, d$Tx, drop = TRUE)
    stopifnot(min(table(cell)) >= k)
    d$bin <- ave(d$body_weight, cell, FUN = function(v)
      pmin(pmax(ceiling(rank(v, ties.method = "first") / length(v) * k), 1), k))
    d$w <- 1 / ave(d$body_weight, cell, d$bin, FUN = length)
    rows <- lapply(seq_len(k), function(q) {
      z <- d[d$bin == q, ]; c <- z[z$Tx == 0, ]; t <- z[z$Tx == 1, ]
      mc <- as.numeric(summary(survfit(Surv(entry, lifespan_days, dead) ~ 1,
                       data = c, weights = w, conf.type = "none"))$table["median"])
      mt <- as.numeric(summary(survfit(Surv(entry, lifespan_days, dead) ~ 1,
                       data = t, weights = w, conf.type = "none"))$table["median"])
      data.frame(analysis = analysis, bins = k, bin = q,
                 mean_control_weight = weighted.mean(c$body_weight, c$w),
                 control_median = mc, treated_median = mt,
                 gain_pct = 100 * (mt - mc) / mc)
    })
    structure(do.call(rbind, rows), bins = d$bin, weights = d$w)
  }

  comparisons <- read.csv(file.path(outdir, "predosing_comparison_mice.csv"),
                          stringsAsFactors = FALSE)
  lower <- comparisons[!duplicated(comparisons$mouse), ]
  nonlower <- read.csv(file.path(outdir, "nonlowering_comparison_mice.csv"),
                       stringsAsFactors = FALSE)

  quarters <- standardized_bins(comparisons, "weight-lowering standardized quarters", 4)
  stopifnot(identical(as.integer(attr(quarters, "bins")), as.integer(comparisons$matched_quartile)))
  stopifnot(isTRUE(all.equal(attr(quarters, "weights"), comparisons$standardization_weight)))
  primary <- read.csv(file.path(outdir, "matched_quartile_standardized_survival.csv"))
  stopifnot(isTRUE(all.equal(quarters$control_median, primary$standardized_control_median)))
  stopifnot(isTRUE(all.equal(quarters$treated_median, primary$standardized_treated_median)))
  stopifnot(isTRUE(all.equal(quarters$mean_control_weight, primary$mean_control_weight)))
  rows <- list(
    summarize_bins(lower, "weight-lowering conventional quartiles", 4),
    summarize_bins(nonlower, "non-lowering conventional quartiles", 4),
    standardized_bins(comparisons, "weight-lowering standardized thirds", 3),
    quarters,
    standardized_bins(comparisons, "weight-lowering standardized fifths", 5)
  )
  rows <- lapply(rows, function(x) { attr(x, "bins") <- NULL; attr(x, "weights") <- NULL; x })
  out <- do.call(rbind, rows)
  write.csv(out, file.path(outdir, "descriptive_sensitivities.csv"), row.names = FALSE)
  print(out, row.names = FALSE, digits = 6)

}


summarize_compound_quartiles <- function() {
  # Descriptive quartiles for the pre-treatment weight-lowering single-compound set.
  suppressPackageStartupMessages(library(survival))

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")
  root <- normalizePath(file.path(dirname(script), "..", ".."), winslash = "/")
  input <- file.path(root, "data", "derived", "predosing_comparison_mice.csv")
  output <- file.path(root, "data", "derived", "compound_descriptive.csv")

  d <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
  d$Tx <- as.integer(d$Tx); d$dead <- as.integer(d$dead)
  d$matched_quartile <- as.integer(d$matched_quartile)
  raw <- sub("^[^|]+\\|", "", d$experiment)
  d$compound <- ifelse(grepl("^17aE2", raw), "17alpha-estradiol",
    ifelse(grepl("^ACA(_|$)", raw), "Acarbose",
    ifelse(grepl("^Cana", raw), "Canagliflozin",
    ifelse(grepl("^[Rr]apa(_|$)", raw), "Rapamycin",
    ifelse(raw == "CC", "Candesartan",
    ifelse(raw == "DMF_9", "Dimethyl fumarate",
    ifelse(raw == "Fis_Cyc", "Fisetin",
    ifelse(raw == "GGA", "Geranylgeranylacetone",
    ifelse(raw == "Gly", "Glycine",
    ifelse(raw == "MIF098", "MIF098",
    ifelse(raw == "Mec", "Meclizine",
    ifelse(raw == "NDGA", "NDGA",
    ifelse(raw == "OH-EST", "16-hydroxyestradiol", NA)))))))))))))
  # Exclude multicomponent formulations from the single-compound subset.
  d$compound[raw == "HBX"] <- "2-(2-Hydroxyphenyl)benzoxazole"
  d$compound[grepl("^Res(_|$)", raw)] <- "Resveratrol"
  stopifnot(all(!is.na(d$compound) | raw %in% c("Prot", "MetRapa", "RaAc_9", "RaAc_16", "PB125")))
  d <- d[!is.na(d$compound), ]

  weighted_median <- function(z) {
    fit <- survfit(Surv(entry, lifespan_days, dead) ~ 1, data = z,
                   weights = standardization_weight, conf.type = "none")
    as.numeric(summary(fit)$table["median"])
  }
  estimate <- function(z) {
    cmed <- weighted_median(z[z$Tx == 0, ])
    tmed <- weighted_median(z[z$Tx == 1, ])
    c(control = cmed, treated = tmed, gain_pct = 100 * (tmed - cmed) / cmed)
  }

  rows <- lapply(sort(unique(d$compound)), function(nm) {
    z <- d[d$compound == nm, ]; lo <- estimate(z[z$matched_quartile == 1, ])
    ov <- estimate(z); hi <- estimate(z[z$matched_quartile == 4, ])
    data.frame(compound = nm, experiments = length(unique(z$experiment)), unique_mice = length(unique(z$mouse)),
      lowest_control_median = lo["control"], lowest_treated_median = lo["treated"], lowest_gain_pct = lo["gain_pct"],
      overall_control_median = ov["control"], overall_treated_median = ov["treated"], overall_gain_pct = ov["gain_pct"],
      highest_control_median = hi["control"], highest_treated_median = hi["treated"], highest_gain_pct = hi["gain_pct"])
  })
  out <- do.call(rbind, rows)
  out <- out[order(-out$overall_gain_pct), ]
  lo <- estimate(d[d$matched_quartile == 1, ]); ov <- estimate(d); hi <- estimate(d[d$matched_quartile == 4, ])
  out <- rbind(out, data.frame(compound = "Pooled single-agent set", experiments = length(unique(d$experiment)),
    unique_mice = length(unique(d$mouse)), lowest_control_median = lo["control"], lowest_treated_median = lo["treated"],
    lowest_gain_pct = lo["gain_pct"], overall_control_median = ov["control"], overall_treated_median = ov["treated"],
    overall_gain_pct = ov["gain_pct"], highest_control_median = hi["control"], highest_treated_median = hi["treated"],
    highest_gain_pct = hi["gain_pct"]))
  stopifnot(nrow(out) == 15, tail(out$experiments, 1) == 26,
            tail(out$unique_mice, 1) == 7238,
            !"Protandim" %in% out$compound,
            sum(out$highest_gain_pct[-nrow(out)] > out$lowest_gain_pct[-nrow(out)]) == 13)
  write.csv(out, output, row.names = FALSE)
  print(out, row.names = FALSE, digits = 4)

}


fit_compound_moderation <- function() {
  # Common treatment-by-weight interaction in the pre-treatment single-compound set.
  suppressPackageStartupMessages(library(survival))

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")
  root <- normalizePath(file.path(dirname(script), "..", ".."), winslash = "/")
  input <- file.path(root, "data", "derived", "predosing_comparison_mice.csv")
  output <- file.path(root, "data", "derived", "compound_common_interaction.csv")

  d <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
  d$Tx <- as.integer(d$Tx); d$dead <- as.integer(d$dead)
  raw <- sub("^[^|]+\\|", "", d$experiment)
  d$compound <- ifelse(grepl("^17aE2", raw), "17alpha-estradiol",
    ifelse(grepl("^ACA(_|$)", raw), "Acarbose",
    ifelse(grepl("^Cana", raw), "Canagliflozin",
    ifelse(grepl("^[Rr]apa(_|$)", raw), "Rapamycin",
    ifelse(raw == "CC", "Candesartan",
    ifelse(raw == "DMF_9", "Dimethyl fumarate",
    ifelse(raw == "Fis_Cyc", "Fisetin",
    ifelse(raw == "GGA", "Geranylgeranylacetone",
    ifelse(raw == "Gly", "Glycine",
    ifelse(raw == "MIF098", "MIF098",
    ifelse(raw == "Mec", "Meclizine",
    ifelse(raw == "NDGA", "NDGA",
    ifelse(raw == "OH-EST", "16-hydroxyestradiol", NA)))))))))))))
  # Exclude multicomponent formulations from the single-compound subset.
  d$compound[raw == "HBX"] <- "2-(2-Hydroxyphenyl)benzoxazole"
  d$compound[grepl("^Res(_|$)", raw)] <- "Resveratrol"
  stopifnot(all(!is.na(d$compound) | raw %in% c("Prot", "MetRapa", "RaAc_9", "RaAc_16", "PB125")))
  d <- d[!is.na(d$compound), ]
  d$experiment <- factor(d$experiment)
  d$Txz <- d$Tx * d$weight_z_matched

  base_formula <- Surv(entry, lifespan_days, dead) ~ 0 +
    Tx:experiment + weight_z_matched:experiment + Txz + strata(stratum) + cluster(mouse)
  fit <- coxph(base_formula, data = d, ties = "efron", singular.ok = TRUE)
  s <- summary(fit)$coefficients["Txz", ]
  ci <- summary(fit)$conf.int["Txz", ]

  out <- data.frame(
    analysis = "pooled_common_interaction",
    estimate = unname(ci["exp(coef)"]),
    lower_95 = unname(ci["lower .95"]),
    upper_95 = unname(ci["upper .95"]),
    p = unname(s["Pr(>|z|)"]),
    compounds = length(unique(d$compound)), experiments = length(unique(d$experiment)),
    unique_mice = length(unique(d$mouse))
  )
  print(out, row.names = FALSE, digits = 7)
  stopifnot(nrow(out) == 1, out$compounds[1] == 14, out$experiments[1] == 26,
            out$unique_mice[1] == 7238, !"C2011|Prot" %in% as.character(d$experiment))
  write.csv(out, output, row.names = FALSE)

}


fit_matched_quartiles()


fit_quartile_robustness()


summarize_quartile_sensitivities()


summarize_compound_quartiles()


fit_compound_moderation()
