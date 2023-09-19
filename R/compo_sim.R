#' @title composite event generator for simulation
#' @description composite event generator for simulation by non-homogenuous Poisson process
#' @return a list of simulated data and parameters
#' @param n_patients  number of patients
#' @param n_vars  number of covariates
#' @param vars_cate  vector of "continuous", "binary"
#' @param true_beta  True parameters
#' @param lambda  baseline death hazard for entire population
#' @param constant_baseline_hazard  whether to use constant baseline hazard; default is TRUE;
#'        otherwise, use Weibull baseline hazard with shape parameter a_shape_weibull
#'        and scale parameter sigma_scale_weibull
#' @param a_shape_weibull  shape parameter of Weibull distribution
#' @param sigma_scale_weibull  scale parameter of Weibull distribution
#' @param sigma_scale_gamma  scale parameter of gamma distribution
#' @param non_linear_hazard  whether to use non-linear hazard function
#' @param non_linear_function  the non-linear function to generate hazard
#' @param seed  seed for random number generation
#' @param verbose  whether to print out summary of the number of recurrent events per patients
compo_sim <- function(
    n_patients = 1000,
    n_vars = 10,
    vars_cate = c(rep("binary", 6), rep("continuous", 4)),
    true_beta = c(0, 0, 0, 0.6, 0, 0, 0.8, 0, 0.7, 0),
    lambda = 1e-2,
    constant_baseline_hazard = TRUE,
    a_shape_weibull = 4,
    sigma_scale_weibull = 4,
    sigma_scale_gamma = 0.05,
    non_linear_hazard = FALSE,
    non_linear_function = NULL,
    seed = 926,
    verbose = FALSE) {
  assertthat::assert_that(assertthat::are_equal(n_vars, length(vars_cate)))
  assertthat::assert_that(assertthat::are_equal(length(true_beta), length(vars_cate)))
  assertthat::assert_that(assertthat::see_if(n_patients > 0, msg = "n must be greater than 0"))
  assertthat::assert_that(
    assertthat::see_if(
      sum(vars_cate %in% c("continuous", "binary")) == n_vars,
      msg = "vars_cate must be a character vector of \"continuous\" and \"binary\" "
    )
  )

  set.seed(seed)

  z <- sapply(1:n_vars, function(t) {
    if (vars_cate[t] == "binary") {
      sample(
        x = c(0, 1),
        size = n_patients,
        replace = TRUE,
        prob = c(0.5, 0.5)
      )
    } else {
      runif(n_patients, min = 0, max = 1)
    }
  })

  ## frailty term
  kesi <- rgamma(n_patients, shape = 1 / sigma_scale_gamma, scale = sigma_scale_gamma)

  ## hazard generation
  if (non_linear_hazard) {
    if (is.null(non_linear_function)) {
      stop("non_linear_function must be specified when non_linear_hazard is TRUE")
    }
    lambdaZ <- apply(z, 1, non_linear_function)
  } else {
    if (n_vars == 1) {
      lambdaZ <- exp(true_beta * z)
    } else {
      lambdaZ <- c(exp(z %*% true_beta))
    }
  }

  lambdaA <- lambda * kesi
  lambdaZZ <- lambdaZ * kesi

  ## generate stopping time for each individual
  t_stop <- sapply(lambdaA, function(x) {
    rexp(1, rate = x)
  })

  ## helper function
  if (constant_baseline_hazard) {
    # In order to have the realastic number of recurrent event, we need to find
    # a suitable constant baseline hazard to make that happen
    baseline_hazard <- lambda
    Lambda0_t <- function(t) {
      t * baseline_hazard
    }

    Finv <- function(u, t_max) {
      u * t_max
    }
  } else {
    baseline_hazard <- 1
    Lambda0_t <- function(t) {
      pweibull(t, shape = a_shape_weibull, scale = sigma_scale_weibull)
    }

    Finv <- function(u, t_max) {
      sigma_scale_weibull * (-log(1 - u * Lambda0_t(t_max)))^(1 / a_shape_weibull)
    }
  }

  ## cumulative hazard at stopping time for each patient
  n_recurrents_per_patients <- sapply(Lambda0_t(t_stop) * lambdaZZ, function(x) rpois(1, lambda = x))

  if (verbose) {
    cat("\nTable of the number of recurrent event per patients:\n")
    print(table(n_recurrents_per_patients))
    cat("\nSummary of the number of recurrent event per patients:\n")
    print(summary(n_recurrents_per_patients))
  }

  ## n_recurrents_per_patients could be 0
  ## then t_per_patients is 0, meaning no composite event, only terminal events
  t_per_patients <- purrr::map2(
    n_recurrents_per_patients,
    t_stop,
    function(x, y) {
      if (x == 0) {
        return(y)
      }
      Finv(sort(runif(x)), y)
    }
  )

  ## generate status for each patient
  status_per_patients <- sapply(n_recurrents_per_patients, function(x) {
    if (x == 0) {
      return(0)
    } else if (x == 1) {
      return(1)
    } else {
      return(c(rep(2, x - 1), 1))
    }
  })

  ## although n_recurrents_per_patients could be 0
  ## it is still in the final return dataset, need to consider this
  n_observed_events_per_patients <- sapply(
    n_recurrents_per_patients, function(x) {
      if (x == 0) {
        return(1)
      }
      return(x)
    }
  )
  colnames(z) <- paste(vars_cate, 1:n_vars, sep = "")
  df_true <- data.frame(
    Id = rep(1:n_patients, times = n_observed_events_per_patients),
    Time = unlist(t_per_patients),
    Status = unlist(status_per_patients),
    z[rep(1:n_patients, times = n_observed_events_per_patients), ]
  )

  return(list(
    dataset = df_true,
    true_beta = true_beta,
    t_stop = t_stop,
    n_patients = n_patients,
    n_recurrents_per_patients = n_recurrents_per_patients,
    constant_baseline_hazard = constant_baseline_hazard,
    baseline_hazard = baseline_hazard,
    lambda = lambda,
    lambdaA = lambdaA,
    lambdaZ = lambdaZ,
    kesi = kesi,
    sigma_scale_gamma = sigma_scale_gamma,
    a_shape_weibull = a_shape_weibull,
    sigma_scale_weibull = sigma_scale_weibull
  ))
}
