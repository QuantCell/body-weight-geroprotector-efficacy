# Individual-level external cohort models.


fit_restriction_models <- function() {
  suppressPackageStartupMessages(library(survival))

  script_path <- sub(
    "^--file=", "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
  )
  root <- normalizePath(file.path(dirname(script_path), "..", ".."))
  d <- read.csv(file.path(root, "data", "derived", "drido_180_landmark_pooled_cr.csv"))

  km_median <- function(frame) {
    unname(summary(survfit(Surv(entry, SurvDays, event) ~ 1, data = frame))$table["median"])
  }

  quartile_rows <- lapply(1:4, function(q) {
    al <- d[d$group == "AL" & d$matched_quarter == q, ]
    cr <- d[d$group == "CR" & d$matched_quarter == q, ]
    al_median <- km_median(al)
    cr_median <- km_median(cr)
    combined <- rbind(al, cr)
    data.frame(
      quarter = q,
      mean_weight_g = mean(combined$weight_g),
      al_n = nrow(al),
      cr_n = nrow(cr),
      al_km_median_days = al_median,
      cr_km_median_days = cr_median,
      gain_days = cr_median - al_median,
      gain_pct = 100 * (cr_median - al_median) / al_median
    )
  })
  quartiles <- do.call(rbind, quartile_rows)

  d$quarter_c <- d$matched_quarter - 1
  d$cr20 <- as.numeric(d$diet == "20")
  d$cr40 <- as.numeric(d$diet == "40")
  d$any_cr_quarter <- as.numeric(d$group == "CR") * d$quarter_c
  fit <- coxph(
    Surv(entry, SurvDays, event) ~ cr20 + cr40 + quarter_c + any_cr_quarter,
    data = d,
    ties = "efron"
  )
  term <- summary(fit)$coefficients["any_cr_quarter", ]
  ci <- summary(fit)$conf.int["any_cr_quarter", ]
  interaction <- data.frame(
    analysis = "Pooled CR treatment-by-weight-quarter trend",
    n = nrow(d),
    events = sum(d$event),
    hr = unname(ci["exp(coef)"]),
    ci_low = unname(ci["lower .95"]),
    ci_high = unname(ci["upper .95"]),
    p = unname(term["Pr(>|z|)"])
  )

  gradient_rows <- lapply(c("AL", "CR"), function(g) {
    frame <- d[d$group == g, ]
    if (g == "CR") {
      model <- coxph(Surv(entry, SurvDays, event) ~ weight_z + factor(diet), data = frame, ties = "efron")
    } else {
      model <- coxph(Surv(entry, SurvDays, event) ~ weight_z, data = frame, ties = "efron")
    }
    coef_row <- summary(model)$coefficients["weight_z", ]
    ci_row <- summary(model)$conf.int["weight_z", ]
    data.frame(
      group = g,
      n = nrow(frame),
      events = sum(frame$event),
      hr = unname(ci_row["exp(coef)"]),
      ci_low = unname(ci_row["lower .95"]),
      ci_high = unname(ci_row["upper .95"]),
      p = unname(coef_row["Pr(>|z|)"])
    )
  })
  gradients <- do.call(rbind, gradient_rows)

  out <- file.path(root, "data", "derived")
  write.csv(quartiles, file.path(out, "drido_pooled_cr_quartiles.csv"), row.names = FALSE)
  write.csv(interaction, file.path(out, "drido_pooled_cr_interaction.csv"), row.names = FALSE)
  write.csv(gradients, file.path(out, "drido_weight_mortality_by_group.csv"), row.names = FALSE)

}


fit_olfactory_models <- function() {
  suppressPackageStartupMessages(library(survival))

  # Olfactory-exposure cohort. Mice were exposed to no odour (water, ZU),
  # adult-female odour (FU) or adult-male odour (MU) from day 3 to day 60. Six-month weight is standardized
  # within sex x exposure x genotype among mice alive beyond the landmark, with follow-up from six months.
  # Control gradients use the no-odour mice, stratified by genotype; the six-week versus six-month model
  # pools all females, stratified by exposure and genotype. Intervals are model-based.

  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
  package <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/",
                           mustWork = TRUE)
  input <- file.path(package, "data", "derived", "olfactory_cohort_mice.csv")
  output <- file.path(package, "data", "derived", "olfactory_cohort_weight_mortality.csv")

  m <- read.csv(input, stringsAsFactors = FALSE)
  entry <- 182.4
  standardize <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v)) / sd(v))
  row <- function(analysis, sex, fit, term, frame) {
    s <- summary(fit)$coefficients
    se <- if ("robust se" %in% colnames(s)) s[term, "robust se"] else s[term, "se(coef)"]
    data.frame(analysis = analysis, sex = sex, n = nrow(frame), deaths = sum(frame$Dead),
               HR_per_SD = exp(s[term, "coef"]), lo = exp(s[term, "coef"] - 1.96 * se),
               hi = exp(s[term, "coef"] + 1.96 * se), p = 2 * pnorm(-abs(s[term, "coef"] / se)),
               dead_by_day_500_pct = NA_real_)
  }

  six <- m[!is.na(m$W6m) & m$Age > entry, ]
  six$weight_z <- standardize(six$W6m, interaction(six$Sex, six$Rx, six$Gene, drop = TRUE))
  rows <- list()
  for (sex in c("F", "M")) {
    controls <- six[six$Sex == sex & six$Rx == "ZU", ]
    fit <- coxph(Surv(rep(entry, nrow(controls)), Age, Dead) ~ weight_z + strata(Gene),
                 data = controls, ties = "efron")
    rows[[length(rows) + 1]] <- row("control_6m", sex, fit, "weight_z", controls)

    all_controls <- m[m$Sex == sex & m$Rx == "ZU", ]
    km <- summary(survfit(Surv(Age, Dead) ~ 1, data = all_controls), times = 500)
    rows[[length(rows) + 1]] <- data.frame(
      analysis = "control_dead_by_day_500", sex = sex, n = nrow(all_controls),
      deaths = sum(all_controls$Dead), HR_per_SD = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_,
      dead_by_day_500_pct = 100 * (1 - km$surv))
  }

  females <- m[m$Sex == "F" & !is.na(m$W6w) & !is.na(m$W6m) & m$Age > entry, ]
  group <- interaction(females$Rx, females$Gene, drop = TRUE)
  females$z6w <- standardize(females$W6w, group)
  females$z6m <- standardize(females$W6m, group)
  fit <- coxph(Surv(rep(entry, nrow(females)), Age, Dead) ~ z6w + z6m + strata(Rx, Gene),
               data = females, ties = "efron")
  rows[[length(rows) + 1]] <- row("females_6w_adjusted_for_6m", "F", fit, "z6w", females)
  rows[[length(rows) + 1]] <- row("females_6m_adjusted_for_6w", "F", fit, "z6m", females)



  result <- do.call(rbind, rows)
  write.csv(result, output, row.names = FALSE)

  print(result, digits = 4)

  # Estimate all twelve early-adult treatment-by-weight interactions.
  interactions <- list()
  for (landmark in c("12weeks", "6months", "12months")) {
    column <- switch(landmark, "12weeks" = "W12w", "6months" = "W6m", "12months" = "W12m")
    start <- switch(landmark, "12weeks" = 84, "6months" = 182.4, "12months" = 364.8)
    for (sex in c("F", "M")) for (exposure in c("FU", "MU")) {
      z <- m[m$Sex == sex & m$Rx %in% c("ZU", exposure) & !is.na(m[[column]]) &
             m[[column]] >= 10 & m[[column]] <= 100 & m$Age > start, ]
      z$Tx <- as.integer(z$Rx == exposure)
      z$weight_z <- standardize(z[[column]], interaction(z$Rx, z$Gene, drop = TRUE))
      stopifnot(!anyDuplicated(z$UMN))
      {
        fit <- coxph(Surv(rep(start, nrow(z)), Age, Dead) ~ Tx * weight_z + strata(Gene),
                     data = z, ties = "efron", robust = FALSE)
        s <- summary(fit)$coefficients["Tx:weight_z", ]
        b <- s["coef"]; se <- s["se(coef)"]
        interactions[[length(interactions) + 1]] <- data.frame(
          weight_age = landmark, sex = sex, exposure = exposure,
          mice = nrow(z), deaths = sum(z$Dead),
          HR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se),
          p = 2 * pnorm(-abs(b / se)))
      }
    }
  }
  write.csv(do.call(rbind, interactions), file.path(package, "data/derived/olfactory_cohort_interactions.csv"), row.names = FALSE)

}


fit_restriction_models()


fit_olfactory_models()
