# Treatment-by-weight moderation and sensitivity analyses.


fit_site_consistency <- function() {
  suppressPackageStartupMessages(library(survival))
  args <- grep("^--file=", commandArgs(FALSE), value=TRUE)
  root <- normalizePath(file.path(dirname(sub("^--file=", "", args[1])), "..", ".."), winslash="/")
  data <- read.csv(file.path(root, "data/derived/primary_model_mice.csv"), stringsAsFactors=FALSE)
  estimates <- list(); tests <- list()
  for (sx in c("f", "m")) {
    d <- data[data$sex == sx, ]
    d$site <- sub("\\|.*$", "", d$stratum)
    results <- lapply(sort(unique(d$site)), function(site) {
      x <- d[d$site == site, ]
      x$experiment_sex <- factor(x$experiment_sex); x$stratum <- factor(x$stratum); x$mouse <- factor(x$mouse)
      fit <- coxph(Surv(entry, lifespan_days, dead) ~ 0 + experiment_sex:Tx + experiment_sex:wz +
                   interaction + strata(stratum) + cluster(mouse), data=x, ties="efron")
      s <- summary(fit)$coefficients["interaction", ]
      data.frame(sex=sx, site=site, experiments=length(unique(x$experiment_sex)), mice=length(unique(x$mouse)),
                 coefficient=s["coef"], se=s["robust se"], HR=exp(s["coef"]),
                 lo=exp(s["coef"]-1.96*s["robust se"]), hi=exp(s["coef"]+1.96*s["robust se"]),
                 p=s["Pr(>|z|)"])
    })
    results <- do.call(rbind, results)
    # Sites contain disjoint mice; site-specific nuisance terms make covariance block diagonal.
    L <- rbind(c(-1, 1, 0), c(-1, 0, 1))
    difference <- L %*% results$coefficient
    V <- L %*% diag(results$se^2) %*% t(L)
    statistic <- drop(t(difference) %*% solve(V) %*% difference)
    estimates[[sx]] <- results
    tests[[sx]] <- data.frame(sex=sx, test="between_site_heterogeneity", chi_square=statistic, df=2,
                             p=pchisq(statistic, 2, lower.tail=FALSE))
  }
  estimates <- do.call(rbind, estimates); tests <- do.call(rbind, tests)
  write.csv(estimates, file.path(root,"data/derived/site_moderation.csv"), row.names=FALSE)
  write.csv(tests, file.path(root,"data/derived/site_heterogeneity.csv"), row.names=FALSE)
  print(estimates, row.names=FALSE); print(tests, row.names=FALSE)
}


fit_primary_moderation <- function() {
  suppressPackageStartupMessages(library(survival))

  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
  package <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
  input <- file.path(package, "data", "derived", "primary_model_mice.csv")
  output <- file.path(package, "data", "derived", "primary_moderation.csv")

  d <- read.csv(input, stringsAsFactors = FALSE)
  d$dead <- as.integer(d$dead); d$Tx <- as.numeric(d$Tx); d$wz <- as.numeric(d$wz)
  d$interaction <- d$Tx * d$wz
  rows <- list()
  for (sx in c("m", "f")) {
    x <- d[d$sex == sx, ]
    x$experiment_sex <- factor(x$experiment_sex); x$stratum <- factor(x$stratum); x$mouse <- factor(x$mouse)
    fit <- coxph(
      Surv(entry, lifespan_days, dead) ~ 0 + experiment_sex:Tx + experiment_sex:wz + interaction +
        strata(stratum) + cluster(mouse),
      data = x, ties = "efron", singular.ok = TRUE
    )
    s <- summary(fit)$coefficients["interaction", ]
    rows[[length(rows) + 1]] <- data.frame(
      sex = sx, experiments = length(unique(x$experiment_sex)), rows = nrow(x),
      unique_mice = length(unique(x$mouse)), deaths = sum(x$dead),
      HR = exp(s["coef"]), lo = exp(s["coef"] - 1.96 * s["robust se"]),
      hi = exp(s["coef"] + 1.96 * s["robust se"]), p = s["Pr(>|z|)"]
    )
  }
  out <- do.call(rbind, rows)
  write.csv(out, output, row.names = FALSE)
  print(out, row.names = FALSE, digits = 6)

}


fit_intervention_classes <- function() {
  suppressPackageStartupMessages(library(survival))

  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
  package <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
  d <- read.csv(file.path(package, "data", "derived", "primary_model_mice.csv"), stringsAsFactors = FALSE)
  w <- read.csv(file.path(package, "data", "derived", "experiment_weight_effects.csv"), stringsAsFactors = FALSE)
  meta <- read.csv(file.path(package, "data", "derived", "experiment_interactions.csv"), stringsAsFactors = FALSE)

  w$post_dosing_flag <- tolower(as.character(w$post_dosing)) == "true"
  w <- w[w$post_dosing_flag, ]
  w$weighted_diff <- w$diff_g * w$n_tx
  agg <- aggregate(cbind(weighted_diff, n_tx) ~ cohort + group, data = w, sum)
  agg$mean_diff_g <- agg$weighted_diff / agg$n_tx
  agg$lowering <- agg$mean_diff_g < -1

  meta$experiment_sex <- paste(meta$cohort, meta$group, meta$sex, sep = "|")
  meta <- merge(meta, unique(w[c("cohort", "group", "init")]), by = c("cohort", "group"), all.x = TRUE)
  meta$predosing <- is.finite(meta$init) & meta$weight_age < meta$init

  d$cohort <- sub("\\|.*$", "", d$experiment_sex)
  d$remainder <- sub("^[^|]+\\|", "", d$experiment_sex)
  d$group <- sub("\\|[^|]+$", "", d$remainder)
  d <- merge(d, agg[c("cohort", "group", "lowering")], by = c("cohort", "group"), all.x = TRUE)
  d$lowering[is.na(d$lowering)] <- FALSE
  d <- merge(d, meta[c("experiment_sex", "predosing")], by = "experiment_sex", all.x = TRUE)
  d$dead <- as.integer(d$dead)
  d$Tx <- as.numeric(d$Tx)
  d$wz <- as.numeric(d$wz)
  d$interaction <- d$Tx * d$wz

  fit_subset <- function(x, label) {
    x$experiment_sex <- factor(x$experiment_sex)
    x$stratum <- factor(x$stratum)
    x$mouse <- factor(x$mouse)
    f <- coxph(
      Surv(entry, lifespan_days, dead) ~ 0 + experiment_sex:Tx + experiment_sex:wz +
        interaction + strata(stratum) + cluster(mouse),
      data = x, ties = "efron", singular.ok = TRUE
    )
    s <- summary(f)$coefficients["interaction", ]
    data.frame(
      analysis = label,
      experiments = length(unique(x$experiment_sex)),
      rows = nrow(x),
      unique_mice = length(unique(x$mouse)),
      HR = exp(s["coef"]),
      lo = exp(s["coef"] - 1.96 * s["robust se"]),
      hi = exp(s["coef"] + 1.96 * s["robust se"]),
      p = s["Pr(>|z|)"]
    )
  }

  fit_contrast <- function(x, label) {
    x$experiment_sex <- factor(x$experiment_sex)
    x$stratum <- factor(x$stratum)
    x$mouse <- factor(x$mouse)
    x$lowering_num <- as.numeric(x$lowering)
    f <- coxph(
      Surv(entry, lifespan_days, dead) ~ 0 + experiment_sex:Tx + experiment_sex:wz +
        interaction + interaction:lowering_num + strata(stratum) + cluster(mouse),
      data = x, ties = "efron", singular.ok = TRUE
    )
    s <- summary(f)$coefficients["interaction:lowering_num", ]
    data.frame(
      analysis = label,
      experiments = length(unique(x$experiment_sex)),
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
    fit_contrast(male_pre, "male_predosing_class_contrast")
  ))

  write.csv(out, file.path(package, "data", "derived", "weight_lowering_split.csv"), row.names = FALSE)
  print(out, row.names = FALSE, digits = 6)

}


fit_landmark_moderation <- function() {
  suppressPackageStartupMessages(library(survival))

  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
  package <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
  input <- file.path(package, "data", "derived", "landmark_model_mice.csv")
  outdir <- file.path(package, "data", "derived")
  d <- read.csv(input, stringsAsFactors = FALSE)
  d$dead <- as.integer(d$dead); d$Tx <- as.numeric(d$Tx)
  d$positive_m <- as.logical(d$positive_m); d$positive_f <- as.logical(d$positive_f)
  d$interaction <- d$Tx * d$weight_z

  fit_pool <- function(x, label) {
    x$armstage <- factor(x$experiment_stage)
    x$stratum <- factor(paste(x$experiment_stage, x$site, x$sex, sep = "|"))
    x$mouse <- factor(x$mouse)
    fit <- coxph(
      Surv(entry, lifespan_days, dead) ~ 0 + armstage:Tx + armstage:weight_z + interaction +
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
  write.csv(main, file.path(outdir, "successful_landmarks.csv"), row.names = FALSE)
  write.csv(lag90, file.path(outdir, "successful_landmarks_lag90.csv"), row.names = FALSE)

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
  write.csv(loo, file.path(outdir, "successful_landmarks_loo12.csv"), row.names = FALSE)
  print(main, row.names = FALSE, digits = 5)

}


fit_site_consistency()


fit_primary_moderation()


fit_intervention_classes()


fit_landmark_moderation()
