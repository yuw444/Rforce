#' fit kaplan meier survival function
#' @param time: time to event
#' @param status: event status, 0 is censored, 1 is terminal event
#' @return a list of unique time and survival probability

km_fit <- function(time, status) {
  unique_time <- sort(unique(time[status == 1]))
  if (unique_time[1] != 0) {
    unique_time <- c(0, unique_time)
  }
  if (unique_time[length(unique_time)] != max(time)) {
    unique_time <- c(unique_time, max(time))
  }

  n_at_risks <- sapply(unique_time, function(x) {
    sum(x <= time)
  })

  n_deaths <- sapply(unique_time, function(x) {
    sum(x == time)
  })

  surv_prop <- rep(-1, length(unique_time))

  surv_prop[1] <- 1

  for (i in 2:length(unique_time)) {
    surv_prop[i] <- surv_prop[i - 1] * (1 - n_deaths[i] / n_at_risks[i])
  }

  return(list(unique_time = unique_time, surv_prop = surv_prop))
}

#' calculate pseudo risk time at the given interval
#' @import assertthat
#' @param Gt: a list of unique time and survival probability of Censoring
#' @param X: the time point associated with the Status
#' @param Status: the event status at the time point X
#' @param low: the lower bound of the interval
#' @param high: the upper bound of the interval
#' @return: the pseudo risk time at the given interval
#' @examples
#' # example code
#'  Gt <- list(unique_time = c(0, 1, 1.8, 3.5, 4),
#'             surv_prop = c(1, 0.8,0.6, 0.3, 0.1))
#' pseudo_risk_time(Gt, X =1.5 , Status = 1, low =0.5, high=1.9)

pseudo_risk_time <- function(Gt, X, Status = 0, low = 1, high = 3) {
  assertthat::assert_that(
    assertthat::see_if(
      low < high &
        low >= 0 &
        high >= 0 &
        X >= 0,
      msg = "recheck the range of \"low\", \"high\" and \"X\" "
    )
  )
  assertthat::assert_that(assertthat::see_if(
    Status %in% c(0, 1),
    msg = " \"Status\" must be either 0 or 1"
  ))
  assertthat::assert_that(assertthat::see_if(
    sum(
      names(Gt) %in% c("unique_time", "surv_prop")
    ) ==
      2,
    msg = "The format of Gt is incorrect"
  ))

  if (Status == 0 && low >= X) {
    return(0)
  }

  if (high <= X) {
    return(high - low)
  }

  if (Status == 0 && low < X) {
    return(X - low)
  }

  if (Status == 1) {
    wt <-
      data.frame(
        t = c(
          0,
          Gt$unique_time[(sum(Gt$unique_time < X) + 1):length(Gt$unique_time)]
        ),
        w = Gt$surv_prop[(sum(Gt$unique_time < X)):length(Gt$unique_time)] /
          Gt$surv_prop[sum(Gt$unique_time < X)]
      )

    low_order <- sum(low >= wt$t)

    wt_used <- wt[low_order:nrow(wt), ]

    wt_used[1, 1] <- low

    high_order <- sum(high >= wt_used$t)

    wt_updated <- wt_used[1:high_order, ]

    # par(mfrow=c(1,2))
    # plot(wt, type = "s", xlim = c(0,5),ylim=c(0,1))
    # plot(wt_updated, type = "s", xlim = c(0,5),ylim=c(0,1))
    # dev.off()

    return(sum(
      diff(c(
        wt_updated$t,
        high
      )) *
        wt_updated$w
    ))
  }
}

#' calculate observed risk time at the given interval
#' @param Gt: a list of unique time and survival probability of Censoring
#' @param X: the time point associated with the Status
#' @param Status: the event status at the time point X
#' @param low: the lower bound of the interval
#' @param high: the upper bound of the interval
#' @return: the observed risk time at the given interval
#'
observed_risk_time <- function(Gt, X, Status = 0, low = 1, high = 3) {
  assertthat::assert_that(
    assertthat::see_if(
      low < high &
        low >= 0 &
        high >= 0 &
        X >= 0,
      msg = "recheck the range of \"low\", \"high\" and \"X\" "
    )
  )
  assertthat::assert_that(assertthat::see_if(
    Status %in% c(0, 1),
    msg = " \"Status\" must be either 0 or 1"
  ))
  assertthat::assert_that(assertthat::see_if(
    sum(
      names(Gt) %in% c("unique_time", "surv_prop")
    ) ==
      2,
    msg = "The format of Gt is incorrect"
  ))

  if (low >= X) {
    return(0)
  }
  if (high <= X) {
    return(high - low)
  }
  if (X >= low) {
    return(X - low)
  }
}

#' break the given length into several intervals with the given intervals
#' @import assertthat
#' @export
#' @param length: the length to be broken
#' @param interval: the given intervals
#' @return: a vector of length of each interval
#' @examples
#' # example code
#'  break_length_by_interval(100, 20:26)
break_length_by_interval <- function(length, interval) {
  assertthat::assert_that(
    assertthat::see_if(
      sum(interval) >= length,
      msg = "provided break points is not suitable for the interval!"
    )
  )
  out <- interval[cumsum(interval) < length]

  if (sum(out) < length) {
    out <- c(out, length - sum(out))
  }

  return(out)
}

#' calculate the number of events in each interval
#' @import assertthat
#' @param event_times: a vector of event times
#' @param ids: a vector of ids
#' @param status: a vector of status
#' @param weights_by_status: a vector of weights by status
#' @param interval: a vector of intervals
#' @return a matrix of number of events in each interval by id, row is id, column is interval
#'
counts_by_interval_and_id <- function(
  event_times,
  ids,
  status,
  weights_by_status,
  interval
) {
  assertthat::assert_that(
    assertthat::see_if(
      length(unique(status)) == length(weights_by_status),
      msg = "The length of weights_by_status must be same with length of unique status!"
    ),
    assertthat::are_equal(length(event_times), length(ids)),
    assertthat::are_equal(length(event_times), length(status)),
    assertthat::are_equal(length(ids), length(status))
  )

  n_interval <- length(interval)
  break_points <- cumsum(interval)
  n_ids <- length(unique(ids))
  unique_status <- sort(unique(status))

  for (s in 1:length(unique_status)) {
    status[status == unique_status[s]] <- weights_by_status[s]
  }

  out <- matrix(nrow = n_ids, ncol = n_interval)

  for (i in 1:n_ids) {
    event_times_per_id <- event_times[ids == i]
    status_per_id <- status[ids == i]
    for (j in 1:n_interval) {
      out[i, j] <-
        sum(status_per_id[event_times_per_id < break_points[j]])
    }
  }

  return(out)
}

validate <- function(object, ...) {
  UseMethod("validate")
}

#' validate data.frame object
#' @import assertthat
#' @param object: a data.frame object
#' @param required_cols: a vector of required column names
#' @return: TRUE if the object is valid
validate.data.frame <- function(
  object,
  required_cols = NULL
) {
  assert_that(is.data.frame(object), msg = "The object must be a data frame!")
  if (!is.null(required_cols)) {
    for (col in required_cols) {
      assert_that(
        col %in% colnames(object),
        msg = paste0("The data frame must contain the column: ", col)
      )
    }
  }
  return(TRUE)
}

#' convert the recorded event time per patient to the number of events per interval per patient
#' @import dplyr assertthat tidyr
#' @export
#' @param data_to_convert: a data frame with columns of Id, X, Status
#' @param units_of_cpiu: a vector of units of CPIU
#' @param weights_by_status: a vector of weights by status, default is c(0,1,1) for censoring (status: 0), terminal event(status: 1) and recurrent event(status: 1)
#' @param pseudo_risk: a boolean value indicating whether to use pseudo risk time
#' @return: a dataframe with number of events in each interval by id, row is id, column is interval
#'
patients_to_cpius <- function(
  data_to_convert,
  units_of_cpiu,
  weights_by_status = c(0, 1, 1),
  pseudo_risk = TRUE,
  wide_format = TRUE
) {
  validate(data_to_convert, required_cols = c("Id", "X", "Status"))

  df_cov <- data_to_convert %>%
    dplyr::select(-c(`Id`, `X`, `Status`))

  lst_cov_summary <- lapply(data, function(col) {
    if (is.factor(col) || is.character(col)) {
      return(
        list(
          type = "categorical",
          levels = unique(col)
        )
      )
    } else {
      return(
        list(
          type = "numeric",
          mean = mean(col, na.rm = TRUE),
          sd = sd(col, na.rm = TRUE),
          min = min(col, na.rm = TRUE),
          max = max(col, na.rm = TRUE)
        )
      )
    }
  })

  list_status <- sort(unique(data_to_convert$Status))

  assert_that(
    length(setdiff(c("0", "1"), as.character(list_status))) == 0,
    msg = 'In `data_to_convert`, Column `Status` must at least contain both elements "0"(censor) and "1"(terminal)!'
  )

  assertthat::assert_that(
    length(list_status) == length(weights_by_status),
    msg = "The length of weights_by_status must be same with length of unique statuses!"
  )

  assertthat::assert_that(
    sum(units_of_cpiu) >= max(data_to_convert$X),
    msg = "units couldn't cover the entire observed time"
  )

  data_to_convert$Events <- data_to_convert$Status

  for (i in seq_along(list_status)) {
    data_to_convert$Events[data_to_convert$Status == list_status[i]] <-
      weights_by_status[i]
  }

  # 0. size of CPIUs, number of patients
  size_cpius <- length(units_of_cpiu)
  n_patients <- length(unique(data_to_convert$Id))

  # create a template based on given units
  # 1. grab the maximum of observed time(X) for each Id
  df_terminal <- data_to_convert %>%
    dplyr::group_by(Id) %>%
    dplyr::filter(row_number() == n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(X, Status)

  if (any(!(df_terminal$Status %in% c(0, 1)))) {
    warning(
      "One or more patients have the last event not being the censoring or terminal event! Considering the last event as the censoring point for IPCW calculation!"
    )
    df_terminal$Status[!(df_terminal$Status %in% c(0, 1))] <- 0
  }

  Gt <- Rforce:::km_fit(df_terminal$X, 1 - df_terminal$Status)

  interval_breaks <- c(0, cumsum(units_of_cpiu))

  # 2. create template and nthInterval and riskTime
  template_for_convert <- data_to_convert %>%
    dplyr::group_by(`Id`) %>%
    dplyr::slice_tail() %>%
    dplyr::ungroup() %>%
    dplyr::slice(rep(seq_len(n()), each = size_cpius)) %>%
    dplyr::mutate(nthInterval = rep(1:size_cpius, times = n_patients))

  if (pseudo_risk == TRUE) {
    risk_time <- sapply(1:n_patients, function(n) {
      sapply(1:size_cpius, function(b) {
        Rforce:::pseudo_risk_time(
          Gt,
          X = df_terminal$X[n],
          Status = df_terminal$Status[n],
          low = interval_breaks[b],
          high = interval_breaks[b + 1]
        )
      })
    })
  } else {
    risk_time <- sapply(1:n_patients, function(n) {
      sapply(1:size_cpius, function(b) {
        observed_risk_time(
          Gt,
          X = df_terminal$X[n],
          Status = df_terminal$Status[n],
          low = interval_breaks[b],
          high = interval_breaks[b + 1]
        )
      })
    })
  }

  # 3. number of events calculation

  df_event <- data_to_convert %>%
    dplyr::group_by(Id) %>%
    dplyr::filter(X <= interval_breaks[2]) %>%
    dplyr::tally(Events) %>%
    dplyr::ungroup()

  if (size_cpius >= 2) {
    for (i in 2:size_cpius) {
      df_event <- data_to_convert %>%
        dplyr::group_by(Id) %>%
        dplyr::filter(X <= interval_breaks[i + 1]) %>%
        dplyr::tally(Events) %>%
        dplyr::ungroup() %>%
        dplyr::left_join(., df_event, "Id")
    }
  }

  colnames(df_event)[1:size_cpius + 1] <-
    paste0("nEvents", size_cpius:1)

  df_event <- df_event %>%
    replace(is.na(.), 0)

  n_events <-
    apply(df_event[, 1:size_cpius + 1], 1, function(x) {
      diff(c(0, rev(x)))
    })

  # 4. data frame to return
  template_for_convert <- template_for_convert %>%
    dplyr::mutate(pseudo_risk_time = c(risk_time), nEvents = c(n_events)) %>%
    dplyr::filter(pseudo_risk_time > 0)

  if (wide_format == TRUE) {
    template_to_return <- template_for_convert %>%
      tidyr::pivot_wider(
        names_from = nthInterval,
        values_from = c(pseudo_risk_time, nEvents)
      ) %>%
      dplyr::mutate(dplyr::across(
        starts_with("pseudo_risk_time"),
        ~ tidyr::replace_na(.x, 0)
      )) %>%
      dplyr::mutate(dplyr::across(
        starts_with("nEvents"),
        ~ tidyr::replace_na(.x, 0)
      ))

    designMatrix_Y <- template_to_return %>%
      dplyr::select(!c(`Id`, `X`, `Status`, `Events`)) %>%
      dplyr::select(!dplyr::starts_with("pseudo_risk_time"))

    auxiliaryFeatures <- template_to_return %>%
      dplyr::select(c(Id, dplyr::starts_with("pseudo_risk_time"))) %>%
      dplyr::mutate(X = df_terminal$X, Status = df_terminal$Status)

    return(structure(
      list(
        data = data_to_convert,
        unitsOfCPIUs = units_of_cpiu,
        nIntervals = size_cpius,
        eventTypes = list_status,
        eventWeights = weights_by_status,
        pseudoRisk = pseudo_risk,
        wideFormat = wide_format,
        designMatrix_Y = designMatrix_Y,
        auxiliaryFeatures = auxiliaryFeatures,
        nPatients = n_patients,
        variableNameOrignal = colnames(df_cov),
        variableSummaryOrignal = lst_cov_summary,
        variableUsed = NULL,
        variableIds = NULL
      ),
      class = "CPIU"
    ))
  }

  designMatrix_Y <- template_for_convert %>%
    dplyr::select(
      !c(`Id`, `X`, `Status`, `Events`, `pseudo_risk_time`)
    )

  auxiliaryFeatures <- template_for_convert %>%
    dplyr::select(c(`Id`, `pseudo_risk_time`))

  return(
    structure(
      list(
        data = data_to_convert,
        unitsOfCPIUs = units_of_cpiu,
        nIntervals = size_cpius,
        eventTypes = list_status,
        eventWeights = weights_by_status,
        pseudoRisk = pseudo_risk,
        wideFormat = wide_format,
        designMatrix_Y = designMatrix_Y,
        auxiliaryFeatures = auxiliaryFeatures,
        nPatients = n_patients,
        variableNameOrignal = colnames(df_cov),
        variableSummaryOrignal = lst_cov_summary,
        isDummy = FALSE,
        variableUsed = NULL,
        variableIds = NULL
      ),
      class = "CPIU"
    )
  )
}

#' validate CPIU object
#' @import assertthat
#' @param object: a CPIU object
#' @return: TRUE if the object is valid
validate.CPIU <- function(object) {
  assert_that(
    inherits(object, "CPIU"),
    msg = "The object must be of class 'CPIU'!"
  )
  validate(object$data, required_cols = c("Id", "X", "Status"))
  if (object$wideFormat == TRUE) {
    assert_that(
      object$nPatients == nrow(object$auxiliaryFeatures) &&
        object$nPatients == nrow(object$designMatrix_Y) &&
        object$nPatients == length(unique(object$auxiliaryFeatures$Id)),
      msg = "In `CPIU` class, .$nPatients must be equal to the number of rows of .$auxiliaryFeatures and .$designMatrix_Y!"
    )
  } else {
    assert_that(
      object$nPatients == length(unique(object$auxiliaryFeatures$Id))
    )
  }
  assert_that(
    is.numeric(object$unitsOfCPIUs),
    msg = "In `CPIU` class, .$unitsOfCPIUs must be a numeric vector!"
  )
  assert_that(
    is.numeric(object$nIntervals),
    msg = "In `CPIU` class, .$nIntervals must be a numeric scalar!"
  )
  assert_that(
    is.numeric(object$eventWeights),
    msg = "In `CPIU` class, .$eventWeights must be a numeric vector!"
  )
  assert_that(
    is.logical(object$pseudoRisk),
    msg = "In `CPIU` class, .$pseudoRisk must be a logical scalar!"
  )
  assert_that(
    is.logical(object$wideFormat),
    msg = "In `CPIU` class, .$wideFormat must be a logical scalar!"
  )
  assert_that(
    is.data.frame(object$designMatrix_Y),
    msg = "In `CPIU` class, .$designMatrix_Y must be a data frame!"
  )

  if (object$wideFormat == TRUE) {
    validate(
      object$designMatrix_Y,
      required_cols = paste0("nEvents_", 1:object$nIntervals)
    )
    validate(object$auxiliaryFeatures, required_cols = c("Id", "X", "Status"))
    assert_that(
      ncol(object$auxiliaryFeatures) == object$nIntervals + 3,
      msg = "In `CPIU` class, the number of columns of .$auxiliaryFeatures must be equal to .$nIntervals + 3!"
    )
    assert_that(
      length(setdiff(object$auxiliaryFeatures$Status, object$eventTypes)) == 0,
      msg = "In `CPIU` class, .$eventTypes must be the same as the unique values of .$auxiliaryFeatures$Status!"
    )
  } else {
    validate(object$designMatrix_Y, required_cols = c("nEvents"))
    validate(
      object$auxiliaryFeatures,
      required_cols = c("Id", "pseudo_risk_time")
    )
  }

  assert_that(
    length(object$unitsOfCPIUs) == object$nIntervals,
    msg = "In `CPIU` class, the length of .$unitsOfCPIUs must be equal to .$nIntervals!"
  )

  assert_that(
    length(object$eventWeights) == length(object$eventTypes),
    msg = "In `CPIU` class, the length of .$eventWeights must be equal to the length of .$eventTypes!"
  )

  assert_that(
    nrow(object$designMatrix_Y) == nrow(object$auxiliaryFeatures),
    msg = "In `CPIU` class, the number of rows of .$designMatrix_Y must be equal to the number of rows of .$auxiliaryFeatures!"
  )

  assert_that(
    length(setdiff(c("0", "1"), as.character(object$eventTypes))) == 0,
    msg = 'In `CPIU` class, .$eventTypes must at least contain both "0"(censor) and "1"(terminal)!'
  )

  return(TRUE)
}

#' convert the recorded event time per patient to the number of events and wt at the given time point per patient
#' @param data_to_convert: a data frame with columns of Id, X, Status
#' @param weights_by_status: a vector of weights by status, default is \code{c(0,1,1)} for censoring (status: 0), terminal event(status: 1) and recurrent event(status: 1)
#' @param time_to_evaluate a scalar, or vector with length equal to \code{length(unique(data_to_convert$Id))}
#'              the time point to calculate the number of events and wt for each patient;
#'
#' @return a data.frame with number of events and wt at the given time point per patient

add_wt <- function(data_to_convert, weights_by_status, time_to_evaluate) {
  validate(data_to_convert, required_cols = c("Id", "X", "Status"))

  list_status <- sort(unique(data_to_convert$Status))

  assert_that(
    length(setdiff(c("0", "1"), as.character(list_status))) == 0,
    msg = 'In `data_to_convert`, Column `Status` must at least contain both elements "0"(censor) and "1"(terminal)!'
  )

  assertthat::assert_that(
    length(list_status) == length(weights_by_status),
    msg = "The length of weights_by_status must be same with length of unique statuses!"
  )

  patient_table <- c(table(data_to_convert$Id))

  if (length(time_to_evaluate) == 1) {
    data_to_convert$time_of_interest <- time_to_evaluate
  } else {
    if (length(time_to_evaluate) == length(patient_table)) {
      data_to_convert$time_of_interest <-
        rep(time_to_evaluate, patient_table)
    } else {
      stop(
        "The length of time_to_evaluate should be 1 or equal to the number of patients"
      )
    }
  }

  data_to_convert$Events <- data_to_convert$Status
  for (i in seq_along(list_status)) {
    data_to_convert$Events[data_to_convert$Status == list_status[i]] <-
      weights_by_status[i]
  }
  n_patients <- length(unique(data_to_convert$Id))
  df_terminal <- data_to_convert %>%
    dplyr::group_by(Id) %>%
    dplyr::arrange(X) %>%
    dplyr::slice_tail() %>%
    dplyr::ungroup()
  Gt <- Rforce:::km_fit(df_terminal$X, 1 - df_terminal$Status)
  Gt_step <-
    stepfun(Gt$unique_time, c(1, Gt$surv_prop), right = FALSE)
  df_terminal$Y_observe <- data_to_convert %>%
    dplyr::group_by(Id) %>%
    dplyr::summarise(Y_observe = sum(Events[X <= time_of_interest])) %>%
    dplyr::pull(Y_observe)
  df_terminal$wt <- apply(df_terminal, 1, function(x) {
    if (x["X"] <= x["time_of_interest"]) {
      if (x["Status"] == 0) {
        return(0)
      } else {
        return(1 / Gt_step(x["X"]))
      }
    } else {
      return(1 / Gt_step(x["time_of_interest"]))
    }
  })

  return(df_terminal)
}

#' Calculate the predicted number of events at given time points
#' @param lambda_pred a matrix of predicted hazard rates at each interval for multiple subjects
#' @param interval_cpius a vector of lengths for each CPIU, same length with `ncol(lambda_pred)`
#' @param time_to_evaluate a scalar, time point to evaluate \eqn{\hat Y}
#' @return a matrix of predicted number of events at each `time_points` for multiple subjects
#' @export
#' @examples
#' lambdas <- matrix(1:10, nrow = 2, byrow = TRUE)
#' intervals <- rep(3,5)
#' add_Y_hat(lambdas, intervals, c(3, 6, 9))
add_Y_hat <- function(lambda_pred, length_cpius, time_to_evaluate) {
  assertthat::assert_that(ncol(lambda_pred) == length(length_cpius))
  assertthat::assert_that(length(time_to_evaluate) == 1)
  assertthat::assert_that(!any(time_to_evaluate > sum(length_cpius)))

  temp <- break_length_by_interval(time_to_evaluate, length_cpius)
  if (length(temp) < length(length_cpius)) {
    temp <- c(temp, rep(0, length(length_cpius) - length(temp)))
  }

  y_hat <- (as.matrix(lambda_pred) %*% temp)[, 1]

  return(y_hat)
}

#' Numerical form to calculate the true Y given the true hazard and other parameters
#' @import cubature
#' @param t the time point to check
#' @param constant_baseline_hazard logical; whether constant baseline hazard is used
#' @param baseline_hazard a scalar, the constant baseline hazard
#' @param a_shape_weibull the parameter of weibull distribution when simulate the data
#' @param sigma_scale_weibull the parameter of weibull distribution when simulate the data
#' @param sigma_scale_gamma the parameter of gamma distribution when simulate the frality term
#' @param lambdaZ a vector, recurrent event hazard rate
#' @param lambda a scalar, the hazard rate of the stopping time
#' @return a scalar, the mean number of event at time `t`
#'
true_Y_numerical_form <- function(
  t,
  constant_baseline_hazard = FALSE,
  baseline_hazard = 1,
  a_shape_weibull,
  sigma_scale_weibull,
  sigma_scale_gamma,
  lambdaZ,
  lambda
) {
  if (constant_baseline_hazard) {
    Lambda0_t <- function(x) min(t, x)
  } else {
    Lambda0_t <- function(x) {
      pweibull(min(t, x), shape = a_shape_weibull, scale = sigma_scale_weibull)
    }
  }

  term <- function(y) {
    y[2] *
      Lambda0_t(y[1]) *
      dexp(y[1], rate = y[2] * lambda) *
      dgamma(y[2], shape = 1 / sigma_scale_gamma, scale = sigma_scale_gamma) *
      baseline_hazard
  }

  out <-
    cubature::adaptIntegrate(
      term,
      lowerLimit = c(0, 0),
      upperLimit = c(Inf, Inf)
    )$integral *
    lambdaZ

  return(out)
}

#' Calculate the mean number events for each subject at different time points in the simulated dataset
#' @param compo_sim_list the object from the `compo_sim` or `compo_sim_mao`
#' @param time_points the time points that need to be evaluated
#' @return a data frame (n_subject x length(time_point))
#' @export
true_Y <- function(compo_sim_list, time_points) {
  rst <- sapply(time_points, function(x) {
    true_Y_numerical_form(
      t = x,
      constant_baseline_hazard = ifelse(
        is.null(compo_sim_list$constant_baseline_hazard),
        FALSE,
        compo_sim_list$constant_baseline_hazard
      ),
      baseline_hazard = ifelse(
        is.null(compo_sim_list$baseline_hazard),
        1,
        compo_sim_list$baseline_hazard
      ),
      a_shape_weibull = ifelse(
        is.null(compo_sim_list$a_shape_weibull),
        NA,
        compo_sim_list$a_shape_weibull
      ),
      sigma_scale_weibull = ifelse(
        is.null(compo_sim_list$sigma_scale_weibull),
        NA,
        compo_sim_list$sigma_scale_weibull
      ),
      sigma_scale_gamma = compo_sim_list$sigma_scale_gamma,
      lambdaZ = compo_sim_list$lambdaZ,
      lambda = compo_sim_list$lambda
    )
  })
  return(rst)
}

#' Empirical Mean Number of Events at the Different Time Point
#' @details
#' Calculate the empirical mean number of event at the different
#' time points by using data simulation
#' @param compo_sim_list simulation list either from \code{compo_sim} or \code{compo_sim_mao}
#' @param weibull_baseline which simulation scheme is used, default is \code{compo_sim}
#' @param x patient characteristics to exam
#' @param time_points time point to exam
#' @param n_sims number of simulations to run
#' @param n_size number of patients at each simulation
#'
empirical_Y <- function(
  compo_sim_list,
  weibull_baseline = TRUE,
  x,
  time_points,
  n_sims,
  n_size
) {
  assertthat::assert_that(length(compo_sim_list$true_beta) == length(x))

  hazard_value <- compo_sim_list$hazard_function(x)

  if (any(class(hazard_value) %in% "matrix")) {
    hazard_value <- hazard_value[1, 1]
  }

  hazard_true <- function(y) return(hazard_value)

  y_out <- matrix(NA, nrow = n_sims, ncol = length(time_points))

  if (weibull_baseline) {
    for (sim in 1:n_sims) {
      df <- Rforce::compo_sim(
        n_patients = n_size,
        non_linear_hazard = TRUE,
        non_linear_function = hazard_true,
        n_vars = 1,
        vars_cate = "binary",
        true_beta = 0,
        seed = sim
      )$dataset

      n_events <- lapply(time_points, function(y) {
        df %>%
          dplyr::filter(`Time` < y) %>%
          dplyr::filter(`Status` != 0) %>%
          dplyr::group_by(`Id`) %>%
          dplyr::summarise(`n` = n()) %>%
          dplyr::pull(`n`) %>%
          sum() /
          n_size
      })

      y_out[sim, ] <- unlist(n_events)
    }
  } else {
    for (sim in 1:n_sims) {
      df <- Rforce::compo_sim_mao(
        n_patient = n_size,
        non_linear_hazard = TRUE,
        non_linear_function = hazard_true,
        n_vars = 1,
        vars_cate = "binary",
        true_beta = 0,
        seed = sim
      )$dataset

      n_events <- lapply(time_points, function(y) {
        df %>%
          dplyr::filter(`Time` < y) %>%
          dplyr::filter(`Status` != 0) %>%
          dplyr::group_by(`Id`) %>%
          dplyr::summarise(`n` = n()) %>%
          dplyr::pull(`n`) %>%
          sum() /
          n_size
      })
      y_out[sim, ] <- unlist(n_events)
    }
  }

  colnames(y_out) <- paste("t_", time_points, sep = "")
  rownames(y_out) <- paste("nsim_", 1:n_sims, sep = "")
  return(y_out)
}

#' Calculate the WRSS
#' @details
#' Calculate the WRSS(Gerds T. 2006) to evaluate the performance of `lambda_pred` when the true Y is unknown
#' @param df_test a data frame contains composite event observations
#' @param weights_by_status a vector, the weights assign to each event type
#' @param lambda_pred a matrix contains the predicated hazard rate for each subject at each interval
#' @param length_cpius a vector, the length of each CPIU
#' @param time_point a scalar, the time point to evaluate WRSS
#' @return a data frame contains WRSS for each subject
#' @export
add_wrss <- function(
  df_test,
  weights_by_status,
  lambda_pred,
  length_cpius,
  time_point
) {
  t1 <- add_wt(df_test, weights_by_status, time_point)

  t1$Y_hat <- add_Y_hat(lambda_pred, length_cpius, time_point)

  t1 <- t1 %>%
    dplyr::mutate(`WRSS` = (`Y_observe` - `Y_hat`)^2 * `wt`)

  return(t1)
}


#' Numerical form to calculate the true Y given the true hazard and other parameters
#' @import cubature
#' @param t the time point to check
#' @param constant_baseline_hazard logical; whether constant baseline hazard is used
#' @param baseline_hazard a scalar, the constant baseline hazard
#' @param a_shape_weibull the parameter of weibull distribution when simulate the data
#' @param sigma_scale_weibull the parameter of weibull distribution when simulate the data
#' @param sigma_scale_gamma the parameter of gamma distribution when simulate the frality term
#' @param lambdaZ a vector, recurrent event hazard rate
#' @param lambda a scalar, the hazard rate of the stopping time
#' @return a scalar, the mean number of event at time `t`
#'
Y_hat_numerical_form <- function(
  t,
  constant_baseline_hazard,
  baseline_hazard = 1,
  mu_0,
  sigma_scale_gamma,
  lambdaZ,
  lambda
) {
  if (constant_baseline_hazard) {
    term <- function(x) {
      x[2] *
        min(t, x[1]) *
        dexp(x[1], rate = x[2] * lambda) *
        dgamma(x[2], shape = 1 / sigma_scale_gamma, scale = sigma_scale_gamma) *
        baseline_hazard
    }

    out <-
      cubature::adaptIntegrate(
        term,
        lowerLimit = c(0, 0),
        upperLimit = c(Inf, Inf)
      )$integral *
      lambdaZ

    return(out)
  }

  term <- function(x) {
    x[2] *
      pweibull(
        min(t, x[1]),
        shape = a_shape_weibull,
        scale = sigma_scale_weibull
      ) *
      dexp(x[1], rate = x[2] * lambda) *
      dgamma(x[2], shape = 1 / sigma_scale_gamma, scale = sigma_scale_gamma) *
      baseline_hazard
  }

  out <-
    adaptIntegrate(
      term,
      lowerLimit = c(0, 0),
      upperLimit = c(Inf, Inf)
    )$integral *
    lambdaZ

  return(out)
}
