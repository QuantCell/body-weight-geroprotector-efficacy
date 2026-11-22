suppressPackageStartupMessages(library(survival))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
package <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

results <- list()
keep <- function(name, table) results[[name]] <<- as.data.frame(table, stringsAsFactors = FALSE)
write_results <- function(tables, folder, combined = list()) {
  write_table <- function(table, name) {
    connection <- file(file.path(folder, paste0(name, ".csv")), "wb")
    write.csv(table, connection, row.names = FALSE, na = "")
    close(connection)
  }
  for (name in setdiff(names(tables), unlist(combined))) write_table(tables[[name]], name)
  for (shared in names(combined)) {
    parts <- tables[combined[[shared]]]
    columns <- unique(unlist(lapply(parts, names)))
    stacked <- do.call(rbind, lapply(names(parts), function(name) {
      table <- parts[[name]]
      for (column in setdiff(columns, names(table))) table[[column]] <- NA
      cbind(model_set = name, table[columns], stringsAsFactors = FALSE)
    }))
    write_table(stacked, shared)
  }
}

local({
  input <- file.path(package, "data", "derived", "quartile_comparisons.csv")

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

  keep("matched_quartile_survival", quartiles)
  keep("quartile_trend", trend)
  print(quartiles, row.names = FALSE, digits = 5)
  print(trend, row.names = FALSE, digits = 5)

  # standardized control mortality by day 500
  dead_by_500 <- function(z) {
    fit <- survfit(Surv(entry, lifespan_days, dead) ~ 1, data = z,
                   weights = standardization_weight, conf.type = "none")
    100 * (1 - summary(fit, times = 500, extend = TRUE)$surv)
  }
  day500 <- do.call(rbind, lapply(1:4, function(q) {
    controls <- d[d$matched_quartile == q & d$Tx == 0, ]
    data.frame(quartile = q, control_dead_by_day_500_pct = dead_by_500(controls))
  }))
  keep("matched_quartile_early_mortality", day500)
  print(day500, row.names = FALSE, digits = 5)

  # the same trend with treatment and quartile effects allowed to vary across age bands
  sb <- survSplit(Surv(entry, lifespan_days, dead) ~ ., data = d, cut = c(450, 600, 750, 900, 1050),
                  start = "t0", end = "t1", event = "dead", episode = "band")
  extra <- c()
  for (j in setdiff(sort(unique(sb$band)), 1)) {
    sb[[paste0("tb", j)]] <- sb$Tx * (sb$band == j); extra <- c(extra, paste0("tb", j))
    for (z in 2:4) { nm <- paste0("q", z, "b", j); sb[[nm]] <- (sb$matched_quartile == z) * (sb$band == j); extra <- c(extra, nm) }
  }
  trend_age_s <- summary(coxph(as.formula(paste(
    "Surv(t0, t1, dead) ~ 0 + experiment:Tx + q_factor + tx_q +", paste(extra, collapse = " + "),
    "+ strata(stratum) + cluster(mouse)")), data = sb, ties = "efron", singular.ok = TRUE))$coefficients["tx_q", ]
  trend_age <- data.frame(
    analysis = "matched quartile trend with treatment and quartile effects varying by age",
    experiments = length(unique(d$experiment)),
    comparison_rows = nrow(d),
    unique_mice = length(unique(d$mouse)),
    HR_per_quartile = exp(trend_age_s["coef"]),
    lo = exp(trend_age_s["coef"] - 1.96 * trend_age_s["robust se"]),
    hi = exp(trend_age_s["coef"] + 1.96 * trend_age_s["robust se"]),
    p = trend_age_s["Pr(>|z|)"]
  )
  keep("quartile_trend_age_varying", trend_age)
  print(trend_age, row.names = FALSE, digits = 5)

  controls <- d[d$Tx == 0, ]
  cuts <- lapply(split(controls$body_weight, as.character(controls$stratum)),
                 quantile, probs = c(0.25, 0.5, 0.75), type = 7)
  common <- d
  common$matched_quartile <- vapply(seq_len(nrow(common)), function(i)
    1L + sum(common$body_weight[i] > cuts[[as.character(common$stratum[i])]]), integer(1))
  filled <- tapply(paste(common$Tx, common$matched_quartile), as.character(common$stratum),
                   function(x) length(unique(x)) == 8)
  common <- droplevels(common[common$stratum %in% names(filled)[filled], ])
  common$q_factor <- factor(common$matched_quartile, levels = 1:4, labels = paste0("Q", 1:4))
  common$tx_q <- common$Tx * (common$matched_quartile - 1)
  common_fit <- coxph(
    Surv(entry, lifespan_days, dead) ~ 0 + experiment:Tx + q_factor + tx_q +
      strata(stratum) + cluster(mouse),
    data = common,
    ties = "efron",
    singular.ok = TRUE
  )
  common_s <- summary(common_fit)$coefficients["tx_q", ]
  common_trend <- data.frame(
    analysis = "control-derived cut-offs in both groups",
    experiments = length(unique(common$experiment)),
    experiment_site_cells = length(unique(common$stratum)),
    comparison_rows = nrow(common),
    unique_mice = length(unique(common$mouse)),
    HR_per_quartile = exp(common_s["coef"]),
    lo = exp(common_s["coef"] - 1.96 * common_s["robust se"]),
    hi = exp(common_s["coef"] + 1.96 * common_s["robust se"]),
    p = common_s["Pr(>|z|)"]
  )
  keep("quartile_trend_control_cutoffs", common_trend)
  print(common_trend, row.names = FALSE, digits = 5)

  horizon <- 1000
  replicates <- 2000L
  entry <- unique(d$entry)
  stopifnot(length(entry) == 1L, all(d$lifespan_days > entry))
  mice <- unique(as.character(d$mouse))
  prepare <- function(z) {
    z$cell <- as.integer(factor(paste(z$stratum, z$Tx, z$matched_quartile)))
    z$mouse_index <- match(as.character(z$mouse), mice)
    z
  }
  survival_summary <- function(time, dead, w, include_median) {
    times <- sort(unique(time))
    at <- match(time, times)
    removed <- as.numeric(rowsum(w, at, reorder = TRUE))
    deaths <- as.numeric(rowsum(w * dead, at, reorder = TRUE))
    survival <- cumprod(1 - deaths / rev(cumsum(rev(removed))))
    inside <- times < horizon
    mean_days <- sum(diff(c(entry, times[inside], horizon)) * c(1, survival[inside]))
    if (!include_median) return(c(rmst = mean_days))

    k <- which(survival <= 0.5 + 1e-12)[1]
    median_days <- times[k]
    if (!is.na(k) && abs(survival[k] - 0.5) < 1e-12) {
      next_death <- which(deaths > 0 & seq_along(times) > k)[1]
      if (!is.na(next_death)) median_days <- mean(times[c(k, next_death)])
    }
    c(rmst = mean_days, median = median_days)
  }
  # Control followed by treated within each quartile.
  estimate <- function(z, multiplier, include_median) {
    raw <- multiplier[z$mouse_index]
    w <- raw / ave(raw, z$cell, FUN = sum)
    arm <- 2L * (z$matched_quartile - 1L) + z$Tx + 1L
    arms <- do.call(cbind, lapply(1:8, function(a) {
      i <- arm == a
      survival_summary(z$lifespan_days[i], z$dead[i], w[i], include_median)
    }))
    as.numeric(t(arms))
  }
  datasets <- list("separate ranks" = prepare(d), "control-derived cut-offs" = prepare(common))
  names_to_fit <- names(datasets)
  point <- setNames(lapply(names_to_fit, function(name)
    estimate(datasets[[name]], rep(1, length(mice)), name == "separate ranks")), names_to_fit)
  draws <- lapply(point, function(p) matrix(NA_real_, replicates, length(p)))
  set.seed(20261005)
  for (b in seq_len(replicates)) {
    multiplier <- rexp(length(mice))
    for (name in names_to_fit) {
      draws[[name]][b, ] <- estimate(datasets[[name]], multiplier, name == "separate ranks")
    }
  }
  interval <- function(x) as.numeric(quantile(x, c(0.025, 0.975)))
  rmst_rows <- list()
  contrast_rows <- list()
  for (name in names(datasets)) {
    z <- datasets[[name]]
    p <- point[[name]]
    b <- draws[[name]]
    for (q in 1:4) {
      control <- z[z$matched_quartile == q & z$Tx == 0, ]
      gain <- interval(b[, 2 * q] - b[, 2 * q - 1])
      rmst_rows[[length(rmst_rows) + 1]] <- data.frame(
        grouping = name, quartile = q, experiments = length(unique(z$experiment)),
        experiment_site_cells = length(unique(z$stratum)), unique_mice = length(unique(z$mouse)),
        mean_control_weight_g = mean(tapply(control$body_weight, as.character(control$stratum), mean)),
        control_rmst_days = p[2 * q - 1], treated_rmst_days = p[2 * q],
        gain_days = p[2 * q] - p[2 * q - 1], gain_lo = gain[1], gain_hi = gain[2]
      )
    }
    contrast <- interval((b[, 8] - b[, 7]) - (b[, 2] - b[, 1]))
    contrast_rows[[length(contrast_rows) + 1]] <- data.frame(
      grouping = name, horizon_age_days = horizon, replicates = replicates,
      heaviest_minus_lightest_days = (p[8] - p[7]) - (p[2] - p[1]), lo = contrast[1], hi = contrast[2]
    )
  }
  rmst_rows <- do.call(rbind, rmst_rows)
  contrast_rows <- do.call(rbind, contrast_rows)
  keep("restricted_mean_survival", rmst_rows)
  keep("restricted_mean_contrasts", contrast_rows)

  medians <- point[["separate ranks"]][9:16]
  median_draws <- draws[["separate ranks"]][, 9:16]
  median_rows <- do.call(rbind, lapply(1:4, function(q) {
    control <- 2 * q - 1
    treated <- 2 * q
    days <- interval(median_draws[, treated] - median_draws[, control])
    pct <- interval(100 * (median_draws[, treated] / median_draws[, control] - 1))
    data.frame(
      quartile = q, replicates = replicates,
      control_median = medians[control], treated_median = medians[treated],
      gain_days = medians[treated] - medians[control], gain_days_lo = days[1], gain_days_hi = days[2],
      gain_pct = 100 * (medians[treated] / medians[control] - 1),
      gain_pct_lo = pct[1], gain_pct_hi = pct[2]
    )
  }))
  keep("median_gain_intervals", median_rows)
  print(rmst_rows, row.names = FALSE, digits = 5)
  print(contrast_rows, row.names = FALSE, digits = 5)
})

# single-compound quartile gains
local({

  input <- file.path(package, "data", "derived", "quartile_comparisons.csv")

  d <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
  d$Tx <- as.integer(d$Tx); d$dead <- as.integer(d$dead)
  d$matched_quartile <- as.integer(d$matched_quartile)
  raw <- sub("^[^|]+\\|", "", d$experiment)
  compound_rules <- c("^17aE2" = "17alpha-estradiol", "^ACA(_|$)" = "Acarbose",
    "^Cana" = "Canagliflozin", "^[Rr]apa(_|$)" = "Rapamycin", "^CC$" = "Candesartan",
    "^DMF_9$" = "Dimethyl fumarate", "^Fis_Cyc$" = "Fisetin", "^GGA$" = "Geranylgeranylacetone",
    "^Gly$" = "Glycine", "^HBX$" = "2-(2-hydroxyphenyl)benzoxazole", "^MIF098$" = "MIF098",
    "^Mec$" = "Meclizine", "^NDGA$" = "NDGA", "^OH-EST$" = "16-hydroxyestradiol",
    "^Res(_|$)" = "Resveratrol")
  d$compound <- NA_character_
  for (pattern in names(compound_rules)) {
    d$compound[is.na(d$compound) & grepl(pattern, raw)] <- compound_rules[[pattern]]
  }
  stopifnot(all(!is.na(d$compound) | raw %in% c("MetRapa", "RaAc_9", "RaAc_16", "Prot")))
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
  stopifnot(nrow(out) == 15, tail(out$unique_mice, 1) == 7238)
  keep("compound_survival", out)
  print(out, row.names = FALSE, digits = 4)
})

# common treatment-by-weight interaction across compounds
local({

  input <- file.path(package, "data", "derived", "quartile_comparisons.csv")

  d <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
  d$Tx <- as.integer(d$Tx); d$dead <- as.integer(d$dead)
  raw <- sub("^[^|]+\\|", "", d$experiment)
  compound_rules <- c("^17aE2" = "17alpha-estradiol", "^ACA(_|$)" = "Acarbose",
    "^Cana" = "Canagliflozin", "^[Rr]apa(_|$)" = "Rapamycin", "^CC$" = "Candesartan",
    "^DMF_9$" = "Dimethyl fumarate", "^Fis_Cyc$" = "Fisetin", "^GGA$" = "Geranylgeranylacetone",
    "^Gly$" = "Glycine", "^HBX$" = "2-(2-hydroxyphenyl)benzoxazole", "^MIF098$" = "MIF098",
    "^Mec$" = "Meclizine", "^NDGA$" = "NDGA", "^OH-EST$" = "16-hydroxyestradiol",
    "^Res(_|$)" = "Resveratrol")
  d$compound <- NA_character_
  for (pattern in names(compound_rules)) {
    d$compound[is.na(d$compound) & grepl(pattern, raw)] <- compound_rules[[pattern]]
  }
  stopifnot(all(!is.na(d$compound) | raw %in% c("MetRapa", "RaAc_9", "RaAc_16", "Prot")))
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
            out$unique_mice[1] == 7238)
  keep("compound_common_interaction", out)
})

# descriptive sensitivity analyses
local({
  km_median <- function(z) {
    fit <- survfit(Surv(entry, lifespan_days, dead) ~ 1, data = z, conf.type = "none")
    as.numeric(summary(fit)$table["median"])
  }

  empirical_bins <- function(v, k) {
    cut(v, breaks = quantile(v, probs = 0:k / k, type = 7),
        include.lowest = TRUE, labels = FALSE)
  }

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

  weighted_median <- function(z) {
    fit <- survfit(Surv(entry, lifespan_days, dead) ~ 1, data = z, weights = w, conf.type = "none")
    as.numeric(summary(fit)$table["median"])
  }

  standardized_bins <- function(d, analysis, k) {
    cell <- interaction(d$stratum, d$Tx, drop = TRUE)
    d$bin <- ave(d$body_weight, cell, FUN = function(v)
      pmin(pmax(ceiling(rank(v, ties.method = "first") / length(v) * k), 1), k))
    d$w <- 1 / ave(d$body_weight, cell, d$bin, FUN = length)
    rows <- lapply(seq_len(k), function(q) {
      z <- d[d$bin == q, ]; c <- z[z$Tx == 0, ]; t <- z[z$Tx == 1, ]
      mc <- weighted_median(c); mt <- weighted_median(t)
      data.frame(analysis = analysis, bins = k, bin = q,
                 mean_control_weight = weighted.mean(c$body_weight, c$w),
                 control_median = mc, treated_median = mt,
                 gain_pct = 100 * (mt - mc) / mc)
    })
    structure(do.call(rbind, rows), bins = d$bin)
  }

  comparisons <- read.csv(file.path(package, "data", "derived", "quartile_comparisons.csv"),
                          stringsAsFactors = FALSE)
  lower <- comparisons[!duplicated(comparisons$mouse), ]
  nonlower <- read.csv(file.path(package, "data", "derived", "quartile_nonlowering.csv"),
                       stringsAsFactors = FALSE)

  quarters <- standardized_bins(comparisons, "weight-lowering standardized quarters", 4)
  stopifnot(identical(as.integer(attr(quarters, "bins")), as.integer(comparisons$matched_quartile)))
  rows <- list(
    summarize_bins(lower, "weight-lowering conventional quartiles", 4),
    summarize_bins(nonlower, "non-lowering conventional quartiles", 4),
    standardized_bins(comparisons, "weight-lowering standardized thirds", 3),
    quarters,
    standardized_bins(comparisons, "weight-lowering standardized fifths", 5)
  )
  rows <- lapply(rows, function(x) { attr(x, "bins") <- NULL; x })
  out <- do.call(rbind, rows)
  keep("quartile_sensitivities", out)
  print(out, row.names = FALSE, digits = 6)
})

write_results(results, file.path(package, "data", "outputs"), combined = list(quartile_tests = c(
  "quartile_trend", "quartile_trend_control_cutoffs", "compound_common_interaction",
  "quartile_trend_age_varying")))
