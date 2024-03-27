#' @title composite event generator for simulation
#' @export
#' @description composite event generator for simulation by non-homogeneous Poisson process
#' @return a list of simulated data and parameters
#' @param n_patients  number of patients
#' @param n_vars  number of covariates
#' @param vars_cate  vector of "continuous", "binary"
#' @param true_beta  True parameters
#' @param lambda  baseline death hazard for entire population
#' @param a_shape_weibull  shape parameter of Weibull distribution
#' @param sigma_scale_weibull  scale parameter of Weibull distribution
#' @param sigma_scale_gamma  scale parameter of gamma distribution
#' @param non_linear_hazard  whether to use non-linear hazard function
#' @param non_linear_function  the non-linear function to generate hazard
#' @param seed  seed for random number generation
#' @param verbose  whether to print out summary of the number of recurrent events per patients
#' @examples
#' # example
#' library(doParallel)
#' registerDoParallel(cores = 16)
#' rst <- foreach(i = 1:48) %dopar%{
#'   data_list <- compo_sim(
#'     n_patients = 200,
#'     seed = i,
#'     verbose = FALSE
#'   )
#'   dim(data_list[[1]])
#'   library(dplyr)
#'   df_train <- data_list[[1]] %>%
#'     dplyr::mutate(X = Time) %>%
#'     dplyr::select(-c("Time"))
#'   estimate_list <- wcompo_est(
#'     data = df_train,
#'     weight = c(1, 1)
#'   )
#'   estimate_list$beta
#' }
#' df_rst <- t(do.call("cbind", rst))
#' colMeans(df_rst)
#' boxplot(df_rst)
compo_sim <- function(
    n_patients = 1000,
    n_vars = 10,
    vars_cate = c(rep("binary", 6), rep("continuous", 4)),
    true_beta = c(0, 0, 0, 0.6, 0, 0, 0.8, 0, 0.7, 0),
    lambda = 1e-2,
    a_shape_weibull = 4,
    sigma_scale_weibull = 4,
    sigma_scale_gamma = 0.05,
    non_linear_hazard = FALSE,
    non_linear_function = NULL,
    seed = 926,
    verbose = FALSE) {
  assertthat::assert_that(assertthat::are_equal(n_vars, length(vars_cate)))
  assertthat::assert_that(assertthat::are_equal(length(true_beta), length(vars_cate)))
  assertthat::assert_that(assertthat::see_if(n_patients > 0, msg = "n_patients must be greater than 0"))
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

    Lambda0_t <- function(t) {
      pweibull(t, shape = a_shape_weibull, scale = sigma_scale_weibull)
    }

    Finv <- function(u, t_max) {
      sigma_scale_weibull * (-log(1 - u * Lambda0_t(t_max)))^(1 / a_shape_weibull)
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
    lambda = lambda,
    lambdaA = lambdaA,
    lambdaZ = lambdaZ,
    kesi = kesi,
    sigma_scale_gamma = sigma_scale_gamma,
    a_shape_weibull = a_shape_weibull,
    sigma_scale_weibull = sigma_scale_weibull
  ))
}

#' @import dplyr
#' @title composite event generator for simulation
#' @description composite event generator for simulation, this simulator is constructed based on Mao. L(2015), it assumes constant baseline hazard.
#' @return a list of simulated data and parameters
#' @param n_patients  number of patients
#' @param n_vars  number of covariates
#' @param vars_cate  vector of "continuous", "binary"
#' @param true_beta  True parameters
#' @param s  number of events one patient can up to have, default is 2000
#' @param non_linear_hazard  whether to use non-linear hazard function
#' @param non_linear_function  the non-linear function to generate hazard
#' @param sigma_scale_gamma  variance of frailty term
#' @param seed  seed for random number generation
#' @export
#' @examples
#' # example code
#' library(doParallel)
#' registerDoParallel(cores = 16)
#' rst <- foreach(i = 1:48) %dopar%{
#'   data_list <- compo_sim_mao(
#'     n_patients = 200,
#'     seed = i
#'   )
#'   library(dplyr)
#'   df_train <- manual_censoring(data_list[[1]], 0.8)
#'
#'   estimate_list <- wcompo_est(
#'     data = df_train,
#'     weight = c(1, 1)
#'   )
#'
#'   estimate_list$beta
#' }
#' df_rst <- t(do.call("cbind", rst))
#' colMeans(df_rst)
#' boxplot(df_rst)
compo_sim_mao <- function(
    n_patients = 1000,
    n_vars = 10,
    vars_cate = c(rep("binary", 6), rep("continuous", 4)),
    true_beta = c(0, 0, 0, 0.6, 0, 0, 0.8, 0, 0.7, 0),
    non_linear_hazard = FALSE,
    non_linear_function = NULL,
    sigma_scale_gamma = 0.25,
    seed = 926,
    s = 2000) {
  assertthat::assert_that(assertthat::are_equal(n_vars, length(vars_cate)))
  assertthat::assert_that(assertthat::are_equal(length(true_beta), length(vars_cate)))
  assertthat::assert_that(assertthat::see_if(n_patients > 0, msg = "n must be greater than 0"))
  assertthat::assert_that(
    assertthat::see_if(
      sum(
        vars_cate %in% c("continuous", "binary")
      ) == n_vars,
      msg = "vars_cate must be a character vector of \"continuous\" and \"binary\" "
    )
  )

  set.seed(seed)

  # covariate --> n Bernoulli trials/uniform distribution/norm noise variable
  z <- sapply(1:n_vars, function(t) {
    if (vars_cate[t] == "binary") {
      sample(
        x = c(0, 1),
        size = n_patients,
        replace = TRUE,
        prob = c(0.5, 0.5)
      )
    } else {
      runif(n_patients, min = 0.2, max = 1)
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
  # print(lambdaZ)
  ## lambda asterisk, to be consistent with the paper simulation details

  temp <- rep(0, n_vars)
  temp[vars_cate == "continuous"] <- 0.2
  lambda <- c(0.3 * exp(temp %*% true_beta))

  ## lambdaA and lambdaZ and add the frailty term
  lambdaA <- lambda * kesi
  lambdaZZ <- lambdaZ * kesi

  U <- runif(n_patients)

  ## this step is trying to accommodate the minimum of F_{\tau^\tilde|Z}
  ## is lambdaA/lambdaZ by L'Hoptial's rule, solving for \tau^\tilde

  f_t_tilde <- function(x) {
    if (U[x] <= lambdaA[x] / (lambdaZZ[x])) {
      return(0)
    } else {
      return(uniroot(
        function(t) {
          (1 - exp(-lambdaA[x] * t)) / (1 - exp(-lambdaZZ[x] * t)) - U[x]
        },
        c(10^-7, 10^7),
        tol = .Machine$double.eps^10
      )$root)
    }
  }

  ## find \tau^\tilde
  t_tilde <- sapply(1:n_patients, function(x) {
    f_t_tilde(x)
  })

  # generate 100 Composite event times individually
  t_comp_matrix <-
    t(sapply(1:s, function(x) {
      rexp(n_patients, rate = lambdaZZ)
    }))

  # calculate the cumulative time for every individual
  # (each column is one patient)
  t_comp_matrix_cumsum <- apply(t_comp_matrix, 2, cumsum)

  # chose the maxmium of \tau^\tilde and t_comp_matrix[1,]
  # as Stopping time for each individual
  t_ast <- t(as.matrix(pmax(t_tilde, t_comp_matrix[1, ])))

  # check if composite event censor by stopping time
  # if so, discard the composite event afterwards
  t_ast_matrix <- t_ast[rep(1, s), ]
  censor_stopping <- t_comp_matrix_cumsum > t_ast_matrix
  t_comp_matrix_cumsum[censor_stopping] <- NA
  Time <- as.vector(t_comp_matrix_cumsum)

  ## create design matrix(including duplicates)
  if (n_vars == 1) {
    temp <- rep(z, each = s)
  } else {
    temp <- z[rep(1:n_patients, each = s), ]
  }

  ## data frame to hold all info
  Id <- rep(1:n_patients, each = s)
  data <- cbind.data.frame(Id, temp, Time, Status = 2)
  colnames(data)[1:n_vars + 1] <- paste(vars_cate, 1:n_vars, sep = "")

  # remove the event for each subject when it is censored by stopping time
  temp <- data %>%
    dplyr::group_by(Id) %>%
    dplyr::summarise(death_index = sum(!is.na(Time))) %>%
    dplyr::mutate(terminal = ifelse(death_index == s, 2, 1))

  data_na_rm <- data[!is.na(data$Time), ]
  data_na_rm$Status[cumsum(temp$death_index)] <- temp$terminal

  if (mean(temp$terminal == 1) != 1) {
    warning(paste0("The event rate should be 1 instead of ", mean(temp$terminal == 1), " in the simulate dataset! Retry with larger s."))
  }

  return(list(
    dataset = data_na_rm,
    lambda = lambda,
    lambdaA = lambdaA,
    lambdaZ = lambdaZ,
    kesi = kesi,
    sigma_scale_gamma = sigma_scale_gamma
  ))
}

