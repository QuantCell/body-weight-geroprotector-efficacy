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

# continuous comparator gradient
local({
  input <- file.path(package, "data", "derived", "comparator_cells.csv")

  d <- read.csv(input, stringsAsFactors = FALSE)
  stopifnot(nrow(d) == 534, sum(d$sex == "m") == 273, sum(d$sex == "f") == 261,
            sum(d$cohort == "C2021" & d$site == "TJL") == 12)
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
    rows[[length(rows) + 1]] <- cbind(sex = sx, analysis = "within_cohort_clustered",
                                      fit_gradient(z, cohort_fe = TRUE, clustered = TRUE))
  }
  result <- do.call(rbind, rows)
  keep("comparator_gradients", result)

  print(result, row.names = FALSE, digits = 6)
})

write_results(results, file.path(package, "data", "outputs"))
