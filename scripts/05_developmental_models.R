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

# early-life weight moderation and maternal adjustment
local({
  indir <- file.path(package, "data", "inputs")
  d <- read.csv(file.path(indir,"itp_gn_puberty_weight_data.csv"),
                stringsAsFactors=FALSE,check.names=FALSE)
  d$Tx <- as.numeric(d$Tx); d$dead <- as.integer(d$dead)
  d$positive_m <- tolower(as.character(d$positive_m)) %in% c("true","1")
  d$positive_f <- tolower(as.character(d$positive_f)) %in% c("true","1")
  d$interaction <- d$Tx*d$weight_42d_z

  fit_pool <- function(x,label) {
    x$experimentsex <- factor(paste(x$experiment,x$sex,sep="|"))
    x$stratum <- factor(paste(x$experiment,x$site,x$sex,sep="|"))
    x$mouse <- factor(x$mouse)
    f <- coxph(Surv(entry,age,dead) ~ 0 + experimentsex:Tx + experimentsex:weight_42d_z +
                 interaction + strata(stratum) + cluster(mouse),
               data=x,ties="efron",singular.ok=TRUE)
    s <- summary(f)$coefficients["interaction",]
    data.frame(label=label,sex=unique(x$sex),experiments=length(unique(x$experiment)),
      unique_mice=length(unique(x$mouse)),HR=exp(s["coef"]),
      lo=exp(s["coef"]-1.96*s["robust se"]),
      hi=exp(s["coef"]+1.96*s["robust se"]),p=s["Pr(>|z|)"])
  }
  o <- list()
  for(sx in c("m","f")) {
    pos <- if(sx=="m") d$positive_m else d$positive_f
    x <- d[d$sex==sx & pos,]
    good <- names(which(tapply(x$Tx,x$experiment,function(z)
      length(unique(z))==2 && min(table(z))>=15)))
    x <- x[x$experiment %in% good,]
    o[[length(o)+1]] <- fit_pool(x,paste("pubertal weight",sx))
  }
  o <- do.call(rbind,o)
  keep("puberty_moderation", o)

  ctl <- read.csv(file.path(indir,"itp_gn_earlylife_controls.csv"),
                  stringsAsFactors=FALSE,check.names=FALSE)
  ctl$dead <- as.integer(ctl$dead)

  mout <- list()
  for(sx in c("m","f")) {
    x <- ctl[ctl$sex==sx & is.finite(ctl$weight_6m) &
               is.finite(ctl$maternal_age_days) & is.finite(ctl$parity),]
    x <- x[x$age>183,]
    x$w6z <- ave(x$weight_6m, x$cohort, x$site, FUN=function(v) as.numeric(scale(v)))
    x$maz <- as.numeric(scale(x$maternal_age_days))
    x$pz <- as.numeric(scale(x$parity))
    f0 <- coxph(Surv(rep(183,nrow(x)),age,dead) ~ w6z + strata(cohort,site),
                data=x,ties="efron")
    f1 <- coxph(Surv(rep(183,nrow(x)),age,dead) ~ w6z + maz + pz +
                  strata(cohort,site),
                data=x,ties="efron")
    get <- function(f) {
      s <- summary(f)$coefficients["w6z",]
      c(HR=unname(exp(s["coef"])),lo=unname(exp(s["coef"]-1.96*s["se(coef)"])),
        hi=unname(exp(s["coef"]+1.96*s["se(coef)"])),p=unname(s["Pr(>|z|)"]))
    }
    a <- get(f0); b <- get(f1)
    mout[[length(mout)+1]] <- data.frame(sex=sx,n=nrow(x),
      HR_weight=a["HR"],lo=a["lo"],hi=a["hi"],p=a["p"],
      HR_adjusted_maternal=b["HR"],adj_lo=b["lo"],adj_hi=b["hi"],adj_p=b["p"])
  }
  mout <- do.call(rbind,mout)
  keep("maternal_adjustment", mout)

  cat("Pubertal-weight treatment moderation:\n"); print(o,row.names=FALSE,digits=4)
  cat("\nMaternal-age/parity sensitivity:\n"); print(mout,row.names=FALSE,digits=4)
})

# early-life and six-month weight, separately and mutually adjusted
local({
  x <- read.csv(file.path(package, "data", "inputs", "itp_gn_earlylife_controls.csv"), stringsAsFactors = FALSE)

  fit_row <- function(data, formula, term, panel, label, sex) {
    f <- coxph(formula, data = data, ties = "efron")
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
  rows[[length(rows) + 1]] <- fit_row(
    z, Surv(rep(183, nrow(z)), age, dead) ~ w14z + w6z + strata(cohort, site),
    "w14z", "decomposition", "1.4-month adjusted", "m"
  )
  rows[[length(rows) + 1]] <- fit_row(
    z, Surv(rep(183, nrow(z)), age, dead) ~ w14z + w6z + strata(cohort, site),
    "w6z", "decomposition", "6-month adjusted", "m"
  )

  out <- do.call(rbind, rows)
  keep("developmental_estimates", out)
  print(out, row.names = FALSE, digits = 5)
})

# weight change between six and twelve months
local({
  input <- file.path(package, "data", "derived", "early_weight_change.csv")
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
  keep("adult_weight_change", result)
})

# control weight-mortality gradient by age at weighing
local({

  input <- file.path(package, "data", "inputs", "control_weight_gradient_input.csv")

  d <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
  d$dead <- as.integer(d$dead)
  d$lifespan_days <- as.numeric(d$lifespan_days)
  d <- d[d$group == "Control" & is.finite(d$dead) & is.finite(d$lifespan_days), ]

  rows <- list()
  for (month in c(6, 12, 18, 24)) {
    column <- paste0("bw_", month)
    entry <- month * 30.4
    d[[column]] <- as.numeric(d[[column]])
    d[[column]][d[[column]] < 10 | d[[column]] > 100] <- NA
    for (sex in c("m", "f")) {
      frame <- d[d$sex == sex & is.finite(d[[column]]) & d$lifespan_days > entry, ]
      frame$weight_z <- ave(frame[[column]], frame$cohort, frame$site,
                            FUN = function(v) (v - mean(v)) / sd(v))
      frame <- frame[is.finite(frame$weight_z), ]
      frame$cohort_site <- interaction(frame$cohort, frame$site, drop = TRUE)
      fit <- coxph(Surv(rep(entry, nrow(frame)), lifespan_days, dead) ~
                     weight_z + strata(cohort_site),
                   data = frame, ties = "efron")
      estimate <- summary(fit)$coefficients["weight_z", ]
      rows[[length(rows) + 1]] <- data.frame(
        weight_age = month, sex = sex, n = nrow(frame),
        HR_per_SD = unname(exp(estimate["coef"])),
        lo = unname(exp(estimate["coef"] - 1.96 * estimate["se(coef)"])),
        hi = unname(exp(estimate["coef"] + 1.96 * estimate["se(coef)"])),
        p = unname(estimate["Pr(>|z|)"])
      )
    }
  }
  result <- do.call(rbind, rows)
  keep("control_weight_gradients", result)
  print(result, row.names = FALSE, digits = 5)
})

write_results(results, file.path(package, "data", "outputs"), combined = list(developmental_weights = c(
  "developmental_estimates", "puberty_moderation")))
