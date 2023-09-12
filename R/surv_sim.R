#' Survival data generation using either weibull or inverse CDF method
#'
#' @param n_patients number of patients
#' @param n_vars number of covariates
#' @param vars_cate vector of "continuous", "binary"
#' @param scale_weibull scale parameter for baseline hazard of weibull distribution
#' @param shape_weibull shape parameter of weibull, shape_weibull < 1 decreasing hazard, shape_weibull > 1 increasing hazard, shape_weibull=1 constant hazard
#' @param beta Fixed effect parameter, could control noise variable by set beta[i] = 0
#' @param rateC censoring rate,
#' @param seed seed
#' @param method "InverseCDF" or "Weibull"
#'
#' @return data frame with columns: Id, HazardWOBaseline, Time, Censor, X, Status, x1, x2, ...
#' @export surv_sim
#' @examples
#' set.seed(926)
#' nsim <- 100
#' beta_true <- c(0, 0.8, -0.8, -1.5, 0, 3, 0, -1, 0.8, 1.5, 3, 0, -0.8, 0, -3)
#' betaHat <- matrix(nrow = nsim, ncol = 15)
#' for (k in 1:nsim) {
#'   df_surv <- surv_sim(
#'     n_patients = 1000,
#'     scale_weibull = 1,
#'     shape_weibull = 1,
#'     beta = beta_true,
#'     rateC = 0.3,
#'     n_vars = 15,
#'     vars_cate = c(rep("binary", 7), rep("continous", 8)),
#'     method = "InverseCDF",
#'     seed = 926 + k
#'   )
#'   fit <- coxph(Surv(X, Status) ~ ., data = df_surv[, -c(1:4)])
#'   betaHat[k, ] <- fit$coef
#' }
#' colMeans(betaHat) - beta_true
surv_sim <- function(
    n_patients,
    n_vars,
    vars_cate,
    scale_weibull,
    shape_weibull,
    beta,
    rateC,
    seed = 926,
    method = c("InverseCDF", "Weibull")) {
  assertthat::assert_that(assertthat::are_equal(n_vars, length(vars_cate)))
  assertthat::assert_that(assertthat::are_equal(length(beta), length(vars_cate)))
  assertthat::assert_that(assertthat::see_if(n_patients > 0, msg = "`n_patients` must be integer greater than 0"))
  assertthat::assert_that(assertthat::see_if(shape_weibull > 0, msg = "`shape_weibull` must be greater than 0"))
  assertthat::assert_that(assertthat::see_if(scale_weibull > 0, msg = "`scale_weibull` must be greater than 0"))
  assertthat::assert_that(assertthat::see_if(rateC > 0, msg = "`rateC` must between 0 and 1"))
  assertthat::assert_that(assertthat::see_if(sum(vars_cate %in% c("continous", "binary")) == n_vars,
                                             msg = "vars_cate must be a character vector of \"continous\" and \"binary\" "
  ))

  method <- match.arg(method)
  set.seed(seed)

  x <- sapply(1:n_vars, function(t) {
    if (beta[t] == 0) {
      rnorm(n_patients)
    } else if (vars_cate[t] == "binary") {
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

  if (method == "InverseCDF") {
    # Weibull latent event times
    v <- runif(n = n_patients)
    event_time <- (-log(v) / (scale_weibull * exp(c(x %*% beta))))^(1 / shape_weibull)
  } else {
    # An alternative draw for event times
    lambda_wiki <- scale_weibull^(-1 / shape_weibull) # change definition of scale_weibull to Wikipedia's
    lambda_prime <- lambda_wiki / exp(c(x %*% beta) / shape_weibull) # re-scale according to beta
    event_time <- rweibull(n_patients, shape = shape_weibull, scale = lambda_prime)
  }

  # censoring rate control
  rateExp <- uniroot(
    function(t) {
      C <- rexp(n = n_patients, rate = t)
      C_percent <- 1 - mean(event_time <= C)
      return(C_percent - rateC)
    },
    interval = c(1e-10, 100)
  )$root

  # censoring times
  C <- rexp(n = n_patients, rate = rateExp)

  # follow-up times and event indicators
  time <- pmin(event_time, C)
  status <- as.numeric(event_time <= C)
  colnames(x) <- paste(vars_cate, 1:n_vars, sep = "")

  # data set
  rst <- cbind.data.frame(
    Id = 1:n_patients,
    HazardWOBaseline = as.numeric(format(exp(c(x %*% beta)), scientific = FALSE)),
    Time = event_time,
    Censor = C,
    X = time,
    Status = status,
    x
  )

  return(rst)
}
