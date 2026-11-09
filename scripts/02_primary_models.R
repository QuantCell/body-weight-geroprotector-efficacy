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

# all-experiment treatment-by-weight interaction
local({
  input <- file.path(package, "data", "derived", "within_experiment_long.csv")

  d <- read.csv(input, stringsAsFactors = FALSE)
  d$dead <- as.integer(d$dead); d$Tx <- as.numeric(d$Tx); d$wz <- as.numeric(d$wz)
  d$interaction <- d$Tx * d$wz
  rows <- list()
  for (sx in c("m", "f")) {
    x <- d[d$sex == sx, ]
    x$experimentsex <- factor(x$experimentsex); x$stratum <- factor(x$stratum); x$mouse <- factor(x$mouse)
    fit <- coxph(
      Surv(entry, lifespan_days, dead) ~ 0 + experimentsex:Tx + experimentsex:wz + interaction +
        strata(stratum) + cluster(mouse),
      data = x, ties = "efron", singular.ok = TRUE
    )
    s <- summary(fit)$coefficients["interaction", ]
    rows[[length(rows) + 1]] <- data.frame(
      sex = sx, experiments = length(unique(x$experimentsex)), rows = nrow(x),
      unique_mice = length(unique(x$mouse)), deaths = sum(x$dead),
      HR = exp(s["coef"]), lo = exp(s["coef"] - 1.96 * s["robust se"]),
      hi = exp(s["coef"] + 1.96 * s["robust se"]), p = s["Pr(>|z|)"]
    )
  }
  out <- do.call(rbind, rows)
  keep("primary_interaction", out)
  print(out, row.names = FALSE, digits = 6)
})

# site consistency of the primary interaction
local({

  input <- file.path(package, "data", "derived", "within_experiment_long.csv")

  d <- read.csv(input, stringsAsFactors = FALSE)
  d$dead <- as.integer(d$dead); d$Tx <- as.numeric(d$Tx); d$wz <- as.numeric(d$wz)
  d$interaction <- d$Tx * d$wz
  d$site <- vapply(strsplit(d$stratum, "|", fixed = TRUE), function(z) z[1], "")
  sites <- sort(unique(d$site))
  rows <- list()
  for (sx in c("m", "f")) {
    for (st in sites) {
      x <- d[d$sex == sx & d$site == st, ]
      x$experimentsex <- factor(x$experimentsex); x$stratum <- factor(x$stratum); x$mouse <- factor(x$mouse)
      fit <- coxph(
        Surv(entry, lifespan_days, dead) ~ 0 + experimentsex:Tx + experimentsex:wz + interaction +
          strata(stratum) + cluster(mouse),
        data = x, ties = "efron", singular.ok = TRUE
      )
      s <- summary(fit)$coefficients["interaction", ]
      rows[[length(rows) + 1]] <- data.frame(
        sex = sx, site = st, experiments = length(unique(x$experimentsex)), unique_mice = length(unique(x$mouse)),
        HR = exp(s["coef"]), lo = exp(s["coef"] - 1.96 * s["robust se"]),
        hi = exp(s["coef"] + 1.96 * s["robust se"]), p = s["Pr(>|z|)"]
      )
    }
    x <- d[d$sex == sx, ]
    x$experimentsite <- factor(paste(x$experimentsex, x$site)); x$stratum <- factor(x$stratum); x$mouse <- factor(x$mouse)
    terms <- paste0("interaction_", sites)
    for (i in seq_along(sites)) x[[terms[i]]] <- x$interaction * (x$site == sites[i])
    fit <- coxph(
      as.formula(paste("Surv(entry, lifespan_days, dead) ~ 0 + experimentsite:Tx + experimentsite:wz +",
                       paste(terms, collapse = " + "), "+ strata(stratum) + cluster(mouse)")),
      data = x, ties = "efron", singular.ok = TRUE
    )
    b <- coef(fit)[terms]; V <- vcov(fit)[terms, terms]
    contrast <- cbind(1, -diag(length(terms) - 1))
    wald <- as.numeric(t(contrast %*% b) %*% solve(contrast %*% V %*% t(contrast)) %*% (contrast %*% b))
    rows[[length(rows) + 1]] <- data.frame(
      sex = sx, site = "heterogeneity", experiments = NA, unique_mice = NA,
      HR = NA, lo = NA, hi = NA, p = pchisq(wald, nrow(contrast), lower.tail = FALSE)
    )
  }
  out <- do.call(rbind, rows)
  keep("site_consistency", out)
  print(out, row.names = FALSE, digits = 4)
})

# weight-lowering and non-lowering experiment classes
local({
  d <- read.csv(file.path(package, "data", "derived", "within_experiment_long.csv"), stringsAsFactors = FALSE)
  w <- read.csv(file.path(package, "data", "derived", "experiment_weight_effects.csv"), stringsAsFactors = FALSE)
  meta <- read.csv(file.path(package, "data", "derived", "experiment_interactions.csv"), stringsAsFactors = FALSE)

  w$post_dosing_flag <- tolower(as.character(w$post_dosing)) == "true"
  w <- w[w$post_dosing_flag, ]
  w$weighted_diff <- w$diff_g * w$n_tx
  agg <- aggregate(cbind(weighted_diff, n_tx) ~ cohort + group, data = w, sum)
  agg$mean_diff_g <- agg$weighted_diff / agg$n_tx
  agg$lowering <- agg$mean_diff_g < -1

  meta$experimentsex <- paste(meta$cohort, meta$group, meta$sex, sep = "|")
  meta <- merge(meta, unique(w[c("cohort", "group", "init")]), by = c("cohort", "group"), all.x = TRUE)
  meta$predosing <- is.finite(meta$init) & meta$weight_age < meta$init

  d$cohort <- sub("\\|.*$", "", d$experimentsex)
  d$remainder <- sub("^[^|]+\\|", "", d$experimentsex)
  d$group <- sub("\\|[^|]+$", "", d$remainder)
  d <- merge(d, agg[c("cohort", "group", "lowering")], by = c("cohort", "group"), all.x = TRUE)
  d$lowering[is.na(d$lowering)] <- FALSE
  d <- merge(d, meta[c("experimentsex", "predosing")], by = "experimentsex", all.x = TRUE)
  d$dead <- as.integer(d$dead)
  d$Tx <- as.numeric(d$Tx)
  d$wz <- as.numeric(d$wz)
  d$interaction <- d$Tx * d$wz

  fit_subset <- function(x, label) {
    x$experimentsex <- factor(x$experimentsex)
    x$stratum <- factor(x$stratum)
    x$mouse <- factor(x$mouse)
    f <- coxph(
      Surv(entry, lifespan_days, dead) ~ 0 + experimentsex:Tx + experimentsex:wz +
        interaction + strata(stratum) + cluster(mouse),
      data = x, ties = "efron", singular.ok = TRUE
    )
    s <- summary(f)$coefficients["interaction", ]
    data.frame(
      analysis = label,
      experiments = length(unique(x$experimentsex)),
      rows = nrow(x),
      unique_mice = length(unique(x$mouse)),
      HR = exp(s["coef"]),
      lo = exp(s["coef"] - 1.96 * s["robust se"]),
      hi = exp(s["coef"] + 1.96 * s["robust se"]),
      p = s["Pr(>|z|)"]
    )
  }

  fit_contrast <- function(x, label) {
    x$experimentsex <- factor(x$experimentsex)
    x$stratum <- factor(x$stratum)
    x$mouse <- factor(x$mouse)
    x$lowering_num <- as.numeric(x$lowering)
    f <- coxph(
      Surv(entry, lifespan_days, dead) ~ 0 + experimentsex:Tx + experimentsex:wz +
        interaction + interaction:lowering_num + strata(stratum) + cluster(mouse),
      data = x, ties = "efron", singular.ok = TRUE
    )
    s <- summary(f)$coefficients["interaction:lowering_num", ]
    data.frame(
      analysis = label,
      experiments = length(unique(x$experimentsex)),
      rows = nrow(x),
      unique_mice = length(unique(x$mouse)),
      HR = exp(s["coef"]),
      lo = exp(s["coef"] - 1.96 * s["robust se"]),
      hi = exp(s["coef"] + 1.96 * s["robust se"]),
      p = s["Pr(>|z|)"]
    )
  }

  male <- d[d$sex == "m", ]
  female <- d[d$sex == "f", ]
  male_pre <- male[male$predosing %in% TRUE, ]

  out <- do.call(rbind, list(
    fit_subset(male[male$lowering, ], "male_weight_lowering"),
    fit_subset(male[!male$lowering, ], "male_nonlowering"),
    fit_contrast(male, "male_class_contrast"),
    fit_subset(female[female$lowering, ], "female_weight_lowering"),
    fit_subset(female[!female$lowering, ], "female_nonlowering"),
    fit_subset(male_pre, "male_predosing_all"),
    fit_subset(male_pre[male_pre$lowering, ], "male_predosing_weight_lowering"),
    fit_subset(male_pre[!male_pre$lowering, ], "male_predosing_nonlowering"),
    fit_contrast(male_pre, "male_predosing_class_contrast"),
    fit_subset(female[female$predosing %in% TRUE, ], "female_predosing_all")
  ))

  keep("weight_class", out)
  print(out, row.names = FALSE, digits = 6)
})

# successful-experiment landmark models
local({
  input <- file.path(package, "data", "derived", "landmark_long.csv")
  d <- read.csv(input, stringsAsFactors = FALSE)
  d$dead <- as.integer(d$dead); d$Tx <- as.numeric(d$Tx)
  d$positive_m <- as.logical(d$positive_m); d$positive_f <- as.logical(d$positive_f)
  d$interaction <- d$Tx * d$weight_z

  fit_pool <- function(x, label) {
    x$experimentstage <- factor(x$experiment_stage)
    x$stratum <- factor(paste(x$experiment_stage, x$site, x$sex, sep = "|"))
    x$mouse <- factor(x$mouse)
    fit <- coxph(
      Surv(entry, lifespan_days, dead) ~ 0 + experimentstage:Tx + experimentstage:weight_z + interaction +
        strata(stratum) + cluster(mouse),
      data = x, ties = "efron", singular.ok = TRUE
    )
    s <- summary(fit)$coefficients["interaction", ]
    data.frame(
      analysis = label, landmark = unique(x$landmark), sex = unique(x$sex),
      experiments = length(unique(x$experiment)), rows = nrow(x), unique_mice = length(unique(x$mouse)),
      HR = exp(s["coef"]), lo = exp(s["coef"] - 1.96 * s["robust se"]),
      hi = exp(s["coef"] + 1.96 * s["robust se"]), p = s["Pr(>|z|)"]
    )
  }

  fit_grid <- function(data, lag_days = 0) {
    if (lag_days > 0) {
      data <- data[data$lifespan_days > data$entry + lag_days, ]
      data$entry <- data$entry + lag_days
    }
    rows <- list()
    for (sx in c("m", "f")) for (lm in sort(unique(data$landmark))) {
      positive <- if (sx == "m") data$positive_m else data$positive_f
      x <- data[data$sex == sx & data$landmark == lm & positive, ]
      good <- names(which(tapply(x$Tx, x$experiment_stage, function(z) length(unique(z)) == 2 && min(table(z)) >= 15)))
      x <- x[x$experiment_stage %in% good, ]
      if (length(unique(x$experiment)) >= 1) rows[[length(rows) + 1]] <- fit_pool(x, ifelse(lag_days, "90-day lag", "published extenders"))
    }
    do.call(rbind, rows)
  }

  main <- fit_grid(d, 0)
  lag90 <- fit_grid(d, 90)
  keep("landmark_models", main)
  keep("landmark_lag90", lag90)

  x0 <- d[d$sex == "m" & d$landmark == 12 & d$positive_m, ]
  loo <- list()
  for (drop in c("none", sort(unique(x0$cohort)))) {
    x <- if (drop == "none") x0 else x0[x0$cohort != drop, ]
    good <- names(which(tapply(x$Tx, x$experiment_stage, function(z) length(unique(z)) == 2 && min(table(z)) >= 15)))
    x <- x[x$experiment_stage %in% good, ]
    z <- fit_pool(x, paste("drop", drop)); z$drop_cohort <- drop
    loo[[length(loo) + 1]] <- z
  }
  loo <- do.call(rbind, loo)
  keep("landmark_leave_one_cohort_out", loo)
  print(main, row.names = FALSE, digits = 5)
})

write_results(results, file.path(package, "data", "outputs"), combined = list(treatment_moderation = c(
  "primary_interaction", "weight_class", "landmark_models", "landmark_lag90",
  "landmark_leave_one_cohort_out")))
