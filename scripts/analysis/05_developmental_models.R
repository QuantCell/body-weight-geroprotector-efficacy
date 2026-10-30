# Developmental weights, adult gain and mortality gradients.


fit_control_weight_gradients <- function() {
  suppressPackageStartupMessages(library(survival))

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")
  root <- normalizePath(file.path(dirname(script), "..", ".."), winslash = "/")
  input <- file.path(root, "data", "derived", "control_weight_gradient_input.csv")
  output <- file.path(root, "data", "derived", "control_weight_gradient_by_age.csv")

  d <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
  d$dead <- as.integer(d$dead)
  d$lifespan_days <- as.numeric(d$lifespan_days)
  d <- d[d$group == "Control" & is.finite(d$dead) & is.finite(d$lifespan_days), ]

  extract <- function(fit, term) {
    estimate <- summary(fit)$coefficients[term, ]
    c(
      HR = unname(exp(estimate["coef"])),
      lo = unname(exp(estimate["coef"] - 1.96 * estimate["se(coef)"])),
      hi = unname(exp(estimate["coef"] + 1.96 * estimate["se(coef)"])),
      p = unname(estimate["Pr(>|z|)"])
    )
  }

  rows <- list()
  for (month in c(6, 12, 18, 24)) {
    column <- paste0("bw_", month)
    entry <- month * 30.4
    d[[column]] <- as.numeric(d[[column]])
    d[[column]][d[[column]] < 10 | d[[column]] > 100] <- NA
    for (sex in c("m", "f")) {
      frame <- d[d$sex == sex & is.finite(d[[column]]) & d$lifespan_days > entry, ]
      frame$weight_z <- ave(frame[[column]], frame$cohort, frame$site, FUN = function(values) {
        (values - mean(values)) / sd(values)
      })
      frame <- frame[is.finite(frame$weight_z), ]
      # Stratify baseline hazards by cohort and site; retain all eligible controls.
      frame$cohort_site <- interaction(frame$cohort, frame$site, drop = TRUE)
      fit <- coxph(
        Surv(rep(entry, nrow(frame)), lifespan_days, dead) ~
          weight_z + strata(cohort_site),
        data = frame, ties = "efron"
      )
      estimate <- extract(fit, "weight_z")
      rows[[length(rows) + 1]] <- data.frame(
        weight_age = month, sex = sex, n = nrow(frame),
        HR_per_SD = estimate["HR"], lo = estimate["lo"],
        hi = estimate["hi"], p = estimate["p"]
      )
    }
  }
  result <- do.call(rbind, rows)
  write.csv(result, output, row.names = FALSE)
  print(result, row.names = FALSE, digits = 5)

}


fit_puberty_moderation <- function() {
  suppressPackageStartupMessages(library(survival))

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")
  root <- normalizePath(file.path(dirname(script), "..", ".."), winslash = "/")
  input <- file.path(root, "data", "inputs", "itp_gn_puberty_weight_data.csv")
  output <- file.path(root, "data", "derived", "puberty_moderation.csv")

  d <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
  d$Tx <- as.numeric(d$Tx)
  d$dead <- as.integer(d$dead)
  d$positive_m <- tolower(as.character(d$positive_m)) %in% c("true", "1")
  d$positive_f <- tolower(as.character(d$positive_f)) %in% c("true", "1")
  reference <- d[d$Tx == 0 & !duplicated(d$mouse), ]
  reference$key <- paste(reference$cohort, reference$site, reference$sex, sep = "|")
  means <- tapply(reference$weight_42d, reference$key, mean)
  sds <- tapply(reference$weight_42d, reference$key, sd)
  key <- paste(d$cohort, d$site, d$sex, sep = "|")
  d$weight_42d_z <- (d$weight_42d - means[key]) / sds[key]
  stopifnot(all(is.finite(d$weight_42d_z)))
  write.csv(data.frame(reference_group = names(means),
                      controls = as.integer(table(reference$key)[names(means)]),
                      mean_weight_g = as.numeric(means), sd_weight_g = as.numeric(sds)),
            file.path(root, "data", "derived", "puberty_weight_reference.csv"),
            row.names = FALSE)
  d$interaction <- d$Tx * d$weight_42d_z

  fit_sex <- function(x, sex) {
    x$experiment_sex <- factor(paste(x$experiment, x$sex, sep = "|"))
    x$stratum <- factor(paste(x$experiment, x$site, x$sex, sep = "|"))
    x$mouse <- factor(x$mouse)
    fit <- coxph(
      Surv(entry, age, dead) ~ 0 + experiment_sex:Tx + experiment_sex:weight_42d_z +
        interaction + strata(stratum) + cluster(mouse),
      data = x, ties = "efron", singular.ok = TRUE
    )
    estimate <- summary(fit)$coefficients["interaction", ]
    data.frame(
      label = paste("pubertal weight", sex), sex = sex,
      experiments = length(unique(x$experiment)), unique_mice = length(unique(x$mouse)),
      HR = exp(estimate["coef"]),
      lo = exp(estimate["coef"] - 1.96 * estimate["robust se"]),
      hi = exp(estimate["coef"] + 1.96 * estimate["robust se"]),
      p = estimate["Pr(>|z|)"]
    )
  }

  rows <- list()
  for (sex in c("m", "f")) {
    eligible <- if (sex == "m") d$positive_m else d$positive_f
    frame <- d[d$sex == sex & eligible, ]
    valid_experiments <- names(which(tapply(frame$Tx, frame$experiment, function(values) {
      length(unique(values)) == 2 && min(table(values)) >= 15
    })))
    rows[[length(rows) + 1]] <- fit_sex(frame[frame$experiment %in% valid_experiments, ], sex)
  }
  result <- do.call(rbind, rows)
  write.csv(result, output, row.names = FALSE)
  print(result, row.names = FALSE, digits = 4)

}


fit_adult_weight_change <- function() {
  suppressPackageStartupMessages(library(survival))

  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
  package <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
  input <- file.path(package, "data", "derived", "adult_weight_change_mice.csv")
  output <- file.path(package, "data", "derived", "early_weight_change.csv")
  d <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
  d$is_control <- tolower(as.character(d$is_control)) == "true"

  pooled_within_sd <- function(x, cell) {
    centered <- x - ave(x, cell, FUN = mean)
    sqrt(sum(centered^2) / (length(x) - length(unique(cell))))
  }

  prepare <- function(sex_value, timing) {
    w <- d[d$sex == sex_value, ]
    if (timing == "primary") {
      w <- w[w$is_control | (!is.na(w$init_months) & w$init_months < 12), ]
    } else if (timing == "treated_primary") {
      w <- w[!w$is_control & !is.na(w$init_months) & w$init_months < 12, ]
    } else if (timing == "after_6") {
      w <- w[w$is_control | (!is.na(w$init_months) & w$init_months > 6), ]
    } else if (timing == "at_or_after_12") {
      w <- w[w$is_control | (!is.na(w$init_months) & w$init_months >= 12), ]
    } else {
      stop(timing)
    }
    event_count <- tapply(w$dead, w$cell, sum)
    valid <- names(event_count)[event_count >= 12]
    w <- w[w$cell %in% valid, ]
    w$change_z <- (w$pct_6_12 - ave(w$pct_6_12, w$cell, FUN = mean)) /
      pooled_within_sd(w$pct_6_12, w$cell)
    w$w6_z <- (w$bw_6 - ave(w$bw_6, w$cell, FUN = mean)) /
      pooled_within_sd(w$bw_6, w$cell)
    w
  }

  fit_one <- function(sex_value, timing, adjusted) {
    w <- prepare(sex_value, timing)
    form <- if (adjusted) {
      Surv(lifespan_days - 365, dead) ~ change_z + w6_z + strata(cell) + cluster(cell)
    } else {
      Surv(lifespan_days - 365, dead) ~ change_z + strata(cell) + cluster(cell)
    }
    fit <- coxph(form, data = w, ties = "efron", model = FALSE, x = FALSE, y = FALSE)
    tab <- summary(fit)$coefficients
    ci <- summary(fit)$conf.int
    row <- data.frame(
      sex = sex_value,
      timing = timing,
      adjusted_for_w6 = adjusted,
      n_mice = nrow(w),
      events = sum(w$dead),
      cells = length(unique(w$cell)),
      cohorts = length(unique(w$cohort)),
      hr = unname(ci["change_z", "exp(coef)"]),
      ci_low = unname(ci["change_z", "lower .95"]),
      ci_high = unname(ci["change_z", "upper .95"]),
      p = unname(tab["change_z", "Pr(>|z|)"]),
      w6_hr = NA_real_,
      w6_ci_low = NA_real_,
      w6_ci_high = NA_real_,
      w6_p = NA_real_
    )
    if (adjusted) {
      row$w6_hr <- unname(ci["w6_z", "exp(coef)"])
      row$w6_ci_low <- unname(ci["w6_z", "lower .95"])
      row$w6_ci_high <- unname(ci["w6_z", "upper .95"])
      row$w6_p <- unname(tab["w6_z", "Pr(>|z|)"])
    }
    print(row)
    row
  }

  specs <- list(
    c("m", "primary", FALSE),
    c("m", "primary", TRUE),
    c("f", "primary", FALSE),
    c("f", "primary", TRUE),
    c("m", "treated_primary", TRUE),
    c("m", "after_6", TRUE),
    c("m", "at_or_after_12", TRUE)
  )

  rows <- lapply(specs, function(s) fit_one(s[[1]], s[[2]], as.logical(s[[3]])))
  result <- do.call(rbind, rows)
  write.csv(result, output, row.names = FALSE)
  cat("Saved", output, "\n")

}


fit_maternal_covariates <- function() {
  # Maternal-age and parity adjustment analysis.
  suppressPackageStartupMessages(library(survival))

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")
  root <- normalizePath(file.path(dirname(script), "..", ".."), winslash = "/")
  input <- file.path(root, "data", "inputs", "itp_gn_earlylife_controls.csv")
  output <- file.path(root, "data", "derived", "maternal_adjustment.csv")

  d <- read.csv(input, check.names = FALSE)
  names(d) <- tolower(names(d))

  required <- c("sex", "age", "dead", "site", "cohort", "weight_6m", "maternal_age_days", "parity")
  missing <- setdiff(required, names(d))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

  extract <- function(fit, term) {
    b <- unname(coef(fit)[term])
    se <- unname(sqrt(vcov(fit)[term, term]))
    c(HR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)))
  }

  rows <- list()
  for (sx in c("m", "f")) {
    z <- d[d$sex == sx & complete.cases(d[, required]), ]
    z <- z[z$age > 183, ]
    z$w6z <- ave(z$weight_6m, z$cohort, z$site, FUN = function(v) as.numeric(scale(v)))
    z$maz <- as.numeric(scale(z$maternal_age_days))
    z$pz <- as.numeric(scale(z$parity))
    f0 <- coxph(Surv(rep(183, nrow(z)), age, dead) ~ w6z + strata(cohort, site), data = z, ties = "efron")
    f1 <- coxph(Surv(rep(183, nrow(z)), age, dead) ~ w6z + maz + pz + strata(cohort, site), data = z, ties = "efron")
    a <- extract(f0, "w6z")
    b <- extract(f1, "w6z")
    rows[[length(rows) + 1]] <- data.frame(
      sex = sx, n = nrow(z),
      HR_weight = a["HR"], lo = a["lo"], hi = a["hi"], p = a["p"],
      HR_adjusted_maternal = b["HR"], adj_lo = b["lo"], adj_hi = b["hi"], adj_p = b["p"]
    )
  }
  out <- do.call(rbind, rows)
  write.csv(out, output, row.names = FALSE)
  print(out)

}


fit_developmental_weights <- function() {
  suppressPackageStartupMessages(library(survival))

  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
  root <- normalizePath(file.path(dirname(script), "../.."), winslash="/")
  input <- file.path(root, "data/inputs/itp_gn_earlylife_controls.csv")
  output <- file.path(root, "data/derived/developmental_figure_estimates.csv")
  x <- read.csv(input, stringsAsFactors=FALSE)

  fit_row <- function(data, formula, term, panel, label, sex, fit = NULL) {
    f <- if (is.null(fit)) coxph(formula, data = data, ties = "efron") else fit
    s <- summary(f)$coefficients[term, ]
    data.frame(
      panel = panel, label = label, sex = sex, n = nrow(data),
      HR = exp(s["coef"]),
      lo = exp(s["coef"] - 1.96 * s["se(coef)"]),
      hi = exp(s["coef"] + 1.96 * s["se(coef)"]),
      p = s["Pr(>|z|)"]
    )
  }

  rows <- list()
  for (sx in c("m", "f")) {
    z <- subset(x, sex == sx & is.finite(weight_42d) & age > 183)
    z$w14z <- ave(z$weight_42d, z$cohort, z$site,
                  FUN = function(v) (v - mean(v)) / sd(v))
    z <- subset(z, is.finite(w14z))
    rows[[length(rows) + 1]] <- fit_row(
      z, Surv(rep(183, nrow(z)), age, dead) ~ w14z + strata(cohort, site),
      "w14z", "age", "1.4", sx
    )
  }

  z <- subset(x, sex == "m" & is.finite(weight_42d) & is.finite(weight_6m))
  z$w14z <- ave(z$weight_42d, z$cohort, z$site,
                FUN = function(v) (v - mean(v)) / sd(v))
  z$w6z <- ave(z$weight_6m, z$cohort, z$site,
               FUN = function(v) (v - mean(v)) / sd(v))
  z <- subset(z, age > 183 & is.finite(w14z) & is.finite(w6z))

  rows[[length(rows) + 1]] <- fit_row(
    z, Surv(rep(183, nrow(z)), age, dead) ~ w14z + strata(cohort, site),
    "w14z", "decomposition", "1.4-month alone", "m"
  )
  rows[[length(rows) + 1]] <- fit_row(
    z, Surv(rep(183, nrow(z)), age, dead) ~ w6z + strata(cohort, site),
    "w6z", "decomposition", "6-month alone", "m"
  )
  joint <- coxph(Surv(rep(183, nrow(z)), age, dead) ~ w14z + w6z + strata(cohort, site), data=z, ties="efron")
  joint_s <- summary(joint)$coefficients
  joint_ci <- exp(confint(joint))
  joint_out <- data.frame(term=c("w42z", "w6z"), n=nrow(z), HR=exp(joint_s[, "coef"]),
                        lo=joint_ci[, 1], hi=joint_ci[, 2], p=joint_s[, "Pr(>|z|)"])
  write.csv(joint_out, file.path(root, "data/derived/puberty_joint.csv"), row.names=FALSE)

  rows[[length(rows) + 1]] <- fit_row(
    z, Surv(rep(183, nrow(z)), age, dead) ~ w14z + w6z + strata(cohort, site),
    "w14z", "decomposition", "1.4-month adjusted", "m", fit=joint
  )
  rows[[length(rows) + 1]] <- fit_row(
    z, Surv(rep(183, nrow(z)), age, dead) ~ w14z + w6z + strata(cohort, site),
    "w6z", "decomposition", "6-month adjusted", "m", fit=joint
  )

  out <- do.call(rbind, rows)
  write.csv(out, output, row.names = FALSE)
  print(out, row.names = FALSE, digits = 5)

}


fit_control_weight_gradients()


fit_puberty_moderation()


fit_adult_weight_change()


fit_maternal_covariates()


fit_developmental_weights()
