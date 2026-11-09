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

# olfactory-exposure cohort
local({

  input <- file.path(package, "data", "derived", "olfactory_mice.csv")

  m <- read.csv(input, stringsAsFactors = FALSE)
  entry <- 182
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
  keep("olfactory_weight_models", result)
  print(result, digits = 4)

  ages <- c(W12w = 84, W6m = 183, W12m = 365)
  labels <- c(W12w = "12 weeks", W6m = "6 months", W12m = "12 months")
  inter <- list()
  for (w in names(ages)) for (sex in c("F", "M")) for (exposure in c("FU", "MU")) {
    d <- m[m$Sex == sex & m$Rx %in% c(exposure, "ZU") & !is.na(m[[w]]) & m$Age > ages[[w]], ]
    d$tx <- as.integer(d$Rx == exposure)
    d$wz <- standardize(d[[w]], interaction(d$Rx, d$Gene, drop = TRUE))
    fit <- coxph(Surv(rep(ages[[w]], nrow(d)), Age, Dead) ~ tx * wz + strata(Gene), data = d, ties = "efron")
    s <- summary(fit)$coefficients["tx:wz", ]
    inter[[length(inter) + 1]] <- data.frame(
      weight_age = labels[[w]], sex = c(F = "Female", M = "Male")[[sex]],
      exposure = c(FU = "Adult-female odour", MU = "Adult-male odour")[[exposure]], mice = nrow(d), deaths = sum(d$Dead),
      interaction_hr = exp(s[["coef"]]), interaction_lo = exp(s[["coef"]] - 1.96 * s[["se(coef)"]]),
      interaction_hi = exp(s[["coef"]] + 1.96 * s[["se(coef)"]]), interaction_p = s[["Pr(>|z|)"]])
  }
  interactions <- do.call(rbind, inter)
  keep("olfactory_interactions", interactions)
  print(interactions, digits = 3)
})

write_results(results, file.path(package, "data", "outputs"))
