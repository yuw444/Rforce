#' @title Fast implementation of Lu Mao's Wcompo method by Kim So Young
#' @description This function is a fast implementation of Lu Mao's Wcompo method by Kim So Young
#'
#' @import survival
#' @import dplyr
#' @param data data frame including Id, X, Status, and covariates
#' @param weight weight for death and recurrent events
#' @return list, including beta estimate, y, and t
#' @export
#' @examples
#' data_list <- compo_sim(
#'   n_patients = 500,
#'   constant_baseline_hazard = FALSE,
#'   seed = 926,
#'   verbose = FALSE
#' )
#' library(dplyr)
#' df_train <- data_list[[1]] %>%
#'   dplyr::mutate(X = Time) %>%
#'   dplyr::select(-c("Time"))
#' estimate_list <- wcompo_est(
#'   data = df_train,
#'   weight = c(1, 1)
#' )
#' estimate_list$beta
wcompo_est <- function(
    data,
    weight) {
  assertthat::assert_that(
    assertthat::see_if(
      !is.null(data$Id) & !is.null(data$X) & !is.null(data$Status)
    ),
    msg = "`data` must have `Id`, `X` and `Status` column"
  )
  n_status <- length(unique(data$Status))
  assertthat::assert_that(
    assertthat::see_if(
      n_status - 1 == length(weight)
    ),
    msg = "The weight must have the same length as the number of unique `Status` in `data`"
  )

  id <- data$Id
  time_temp <- data$X
  status_temp <- data$Status
  z_temp <- data %>% dplyr::select(-c("Id", "X", "Status"))

  unique.id <- unique(id)
  n <- length(unique.id)
  time_unique <- unique(time_temp)
  fail <- time_unique[order(time_unique)]

  m <- length(fail)

  death <- ifelse(status_temp == 1, 1, 0)
  recurrent <- ifelse(status_temp == 2, 1, 0)
  event <- death * weight[1] + recurrent * weight[2]
  E <- cbind(id, recurrent, time_temp, death, z_temp, event)
  last.obs <- E[E[, 2] == 0, ]

  p <- dim(z_temp)[2]
  time1 <- last.obs[, 3]
  delta1 <- last.obs[, 4]
  z1 <- last.obs[, 5:(p + 4)]
  dN1 <- matrix(0, n, m)

  for (i in 1:n) {
    time.event.id <- E[E$id == unique.id[i], 3]
    event.w <- E[E$id == unique.id[i], (p + 5)]
    n.event <- length(event.w)

    if (n.event == 1) {
      dN1[i, ] <- (match(fail, time.event.id)) * event.w
      dN1[i, is.na(dN1[i, ])] <- 0
    } else {
      event.w1 <- event.w[-n.event]
      event.w2 <- event.w[n.event]
      dN1[i, ] <- (match(fail, time.event.id) != n.event) * event.w1[1] + (match(fail, time.event.id) == n.event) * event.w2
      dN1[i, is.na(dN1[i, ])] <- 0
    }
  }

  sorting <- order(time1)
  time <- time1[sorting]
  delta <- delta1[sorting]
  z <- as.matrix(z1[sorting, ])
  p <- dim(z)[2]
  dN <- dN1[sorting, ]
  censor <- 1 - delta

  fit.cox <- survival::coxph(Surv(time, censor) ~ z) ###
  gg3 <- (matrix(rep(time, m), n, m) >= matrix(rep(t(fail), each = n), n, m)) ## index matrix Ti >= t

  Gcweight <- matrix(0, n, m)
  gamma1 <- summary(fit.cox)$coef[, 1]
  lambda0 <- survival::basehaz(fit.cox, center = F)

  for (i in 1:n) {
    surv.i <- exp(-lambda0[, 1] %*% t(exp(z[i, ] %*% gamma1)))
    kmest <- stepfun(lambda0$time, c(1, surv.i))
    de <- (matrix(rep(kmest(time), m), n, m))
    de <- de + (de == 0)
    foo.0 <- matrix(rep(t(kmest(fail)), each = n), n, m) / de ## matrix of G(t)/G(Ti) |Z=0
    foo.0[gg3] <- 1 ## keep the last matrix as 1 if Ti >= t
    if (delta[i] == 0) {
      Gcweight[i, ] <- (fail < time[i]) * foo.0[i, ]
    } else {
      Gcweight[i, ] <- foo.0[i, ]
    }
  }

  ###### beta estimate ########

  beta <- rep(0, p)
  delta0 <- 4
  step <- 0
  while (delta0 > 0.000001 & (step <= 20)) {
    step <- step + 1
    beta <- as.vector(beta)
    expz <- exp(z %*% beta) # n*1
    zexpz <- matrix(rep(expz, p), n, p) * z # n*p
    temp0 <- t(expz) %*% (Gcweight) # 1*L1;  ##adding cohort weight to the S0
    S0hat <- temp0 + (temp0 == 0)
    S1hat <- t(zexpz) %*% (Gcweight) ## adding cohort weight to the S1
    S1overS0hat <- S1hat / (matrix(rep(S0hat, each = p), p, m)) # p*m

    ## calculate information next
    Gcweight2 <- matrix(kronecker(rep(1, p), Gcweight), n, m * p) # n x mp
    S2hat <- t(z) %*% (matrix(rep(zexpz, m), n, m * p) * Gcweight2) # (p x n) x (n x mp) = p x mp   t(z2expz)%*%(Gcweight)  #p*n x n*(m*p) = p* (m*p)
    S2overS0hat <- S2hat / matrix(kronecker(S0hat, matrix(rep(1, p * p), p, p)), p, m * p) # p x mp

    # Score function
    tempU <- NULL
    for (pp in 1:p) {
      matrixU <- (matrix(rep(z[, pp], m), n, m) - matrix(rep(S1overS0hat[pp, ], n), n, m, byrow = T)) * dN
      tempU[pp] <- sum(matrixU)
    }

    dN.bar <- colSums(dN) # 1*m
    IpartIhat <- S2overS0hat %*% kronecker(dN.bar, diag(1, p, p)) # p*p
    IpartIIhat <- (t(matrix(rep(dN.bar, p), m, p)) * S1overS0hat) %*% t(S1overS0hat) # p*p


    # Information matrix
    tempI <- IpartIhat - IpartIIhat


    ## newton Raphlson method
    iI <- solve(tempI)
    beta <- beta + iI %*% tempU # calculating beta estimate
    delta0 <- max(abs(iI %*% tempU))
  }
  expz1 <- exp(z %*% beta)
  temp01 <- t(expz1) %*% (Gcweight)
  S0hat1 <- temp01 + (temp01 == 0)
  mu0.t <- cumsum(colSums(dN) / S0hat1)

  return(
    list(
      beta = beta,
      y = mu0.t,
      t = fail
    )
  )
}


