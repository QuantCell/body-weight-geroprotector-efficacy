# Comparator models use cohort × site × intervention cells, analysed separately by sex.


summarize_comparator_splits <- function() {
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
  package <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
  input <- file.path(package, "data", "derived", "comparator_effect_cells.csv")
  output <- file.path(package, "data", "derived", "comparator_split.csv")
  d <- read.csv(input, stringsAsFactors = FALSE)

  # Form empirical quantile bins of control-group mean weight separately by sex.
  rows <- list()
  for (sx in c("m", "f")) {
    x <- d[d$sex == sx, ]
    for (k in c(3, 4, 5, 10)) {
      breaks <- unique(quantile(x$ctl_bw6, probs = 0:k / k, type = 7, na.rm = TRUE))
      if (length(breaks) != k + 1) stop("Non-unique empirical quantile cut points")
      x$bin <- cut(x$ctl_bw6, breaks = breaks, include.lowest = TRUE,
                   labels = paste0("Q", seq_len(k)))
      for (b in levels(x$bin)) {
        z <- x[x$bin == b, ]; w <- 1 / z$se^2
        beta <- sum(z$log_hr * w) / sum(w); se <- sqrt(1 / sum(w))
        rows[[length(rows) + 1]] <- data.frame(
          sex = sx, split = k, bin = b, cells = nrow(z),
          mean_control_weight = mean(z$ctl_bw6), HR = exp(beta),
          reduction_pct = 100 * (1 - exp(beta)),
          lo = exp(beta - 1.96 * se), hi = exp(beta + 1.96 * se)
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  write.csv(out, output, row.names = FALSE)
  print(out[out$split == 3 & out$bin %in% c("Q1", "Q3"), ], row.names = FALSE, digits = 6)

}


fit_comparator_gradients <- function() {
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
  package <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
  input <- file.path(package, "data", "derived", "comparator_effect_cells.csv")
  output <- file.path(package, "data", "derived", "comparator_gradients.csv")

  d <- read.csv(input, stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(d[c("cohort", "site", "sex", "group")]))
  d$cluster <- paste(d$cohort, d$site, sep = "|")

  fit_gradient <- function(z, cohort_fe = FALSE, clustered = FALSE) {
    x <- z$ctl_bw6
    y <- z$log_hr
    w <- 1 / z$se^2
    X <- cbind(intercept = 1, weight = x)
    if (cohort_fe) {
      X <- cbind(X, model.matrix(~ factor(cohort), data = z)[, -1, drop = FALSE])
    }
    bread <- solve(crossprod(X, w * X))
    beta <- bread %*% crossprod(X, w * y)
    residual <- as.vector(y - X %*% beta)
    if (clustered) {
      groups <- unique(z$cluster)
      meat <- matrix(0, ncol(X), ncol(X))
      for (group in groups) {
        selected <- z$cluster == group
        score <- crossprod(X[selected, , drop = FALSE], w[selected] * residual[selected])
        meat <- meat + score %*% t(score)
      }
      covariance <- bread %*% meat %*% bread * length(groups) / (length(groups) - 1)
      degrees <- length(groups) - 1
      critical <- qt(0.975, degrees)
      p <- 2 * pt(-abs(beta[2] / sqrt(covariance[2, 2])), degrees)
    } else {
      covariance <- bread
      critical <- qnorm(0.975)
      p <- 2 * pnorm(-abs(beta[2] / sqrt(covariance[2, 2])))
    }
    slope_se <- sqrt(covariance[2, 2])
    weight_sd <- sd(x)
    data.frame(
      cells = nrow(z), control_clusters = length(unique(z$cluster)),
      hr_per_sd = exp(beta[2] * weight_sd),
      ci_low = exp((beta[2] - critical * slope_se) * weight_sd),
      ci_high = exp((beta[2] + critical * slope_se) * weight_sd),
      p = as.numeric(p)
    )
  }

  rows <- list()
  for (sx in c("m", "f")) {
    z <- d[d$sex == sx, ]
    rows[[length(rows) + 1]] <- cbind(sex = sx, analysis = "within_cohort_shared_control_clustered",
                                      fit_gradient(z, cohort_fe = TRUE, clustered = TRUE))
  }
  result <- do.call(rbind, rows)
  write.csv(result, output, row.names = FALSE)
  print(result, row.names = FALSE, digits = 6)

}


summarize_comparator_splits()


fit_comparator_gradients()
