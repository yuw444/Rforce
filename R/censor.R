#' @title apply administrative censoring in the given dataset
#'
#' @import dplyr purrr
#' @description censoring the given dataset based on the given quantile of observed time
#' @param data_to_convert the dataset to be converted; including the columns of `Id`, `Status` and `Time`
#' @param censoring_quantile the quantile of observed time to be censored
#' @return the censored dataset
#' @export
#' @examples
#' list_data_to_convert <- compo_sim()
#' df_converted <- manual_censoring(
#'     list_data_to_convert$dataset,
#'     0.9
#' )
#' str(df_converted)
#'
admin_censoring <- function(data_to_convert,
                            censoring_quantile) {
  assertthat::assert_that("Id" %in% colnames(data_to_convert))
  assertthat::assert_that("Time" %in% colnames(data_to_convert))
  assertthat::assert_that("Status" %in% colnames(data_to_convert))

  # censoring_time <- quantile(data_to_convert$Time, probs = censoring_quantile)
  censoring_time <- quantile(
    data_to_convert %>%
      dplyr::filter(Status != 0) %>%
      dplyr::select(Time) %>%
      unlist(),
    probs = censoring_quantile
  )
  tt <- data_to_convert %>%
    dplyr::mutate(censored = Time > censoring_time,
                  newTime = pmin(Time, censoring_time)) %>%
    dplyr::group_by(Id) %>%
    dplyr::mutate(nthCensor = cumsum(censored)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(nthCensor == 1 | censored == FALSE)

  Status <- purrr::map2(tt$Status,
                        tt$censored,
                        function(x, y) {
                          if (y == TRUE) {
                            return(0)
                          } else {
                            return(x)
                          }
                        })

  df_out <- data.frame(
    Id = tt$Id,
    X = tt$newTime,
    Status = unlist(Status),
    tt %>% dplyr::select(
      !c("Id", "Time", "Status", "censored", "newTime", "nthCensor")
    )
  )

  return(df_out)
}

#' @title apply random censoring in the given dataset
#'
#' @import dplyr purrr
#' @description censoring the given dataset based on the given event rate
#' @param data_to_convert the dataset to be converted; including the columns of `Id`, `Status` and `Time`
#' @param event_rate the desired event rate
#' @param verbose print out details
#' @return the censored dataset
#' @export
#' @examples
#' list_data_to_convert <- compo_sim()
#' df_converted <- random_censoring(
#'     list_data_to_convert$dataset,
#'     0.8
#' )
#' str(df_converted)
#'

random_censoring <- function(data_to_convert,
                             event_rate,
                             verbose = FALSE) {
  assertthat::assert_that("Id" %in% colnames(data_to_convert))
  assertthat::assert_that("Time" %in% colnames(data_to_convert))
  assertthat::assert_that("Status" %in% colnames(data_to_convert))

  df_terminal <- data_to_convert %>%
    dplyr::group_by(`Id`) %>%
    dplyr::arrange(`Time`) %>%
    dplyr::slice_tail() %>%
    dplyr::ungroup()

  n_patients <- length(unique(data_to_convert$Id))
  n_terminals <- sum(data_to_convert$Status == 1)
  event_rate_current <- n_terminals / n_patients
  if (verbose)
    print(paste0("The actual event rate is ", event_rate_current))

  if (event_rate_current <= event_rate) {
    warning(
      paste0(
        "The actual event rate is lower than the desired event rate,",
        " reset event_rate to the actual event rate ",
        event_rate_current
      )
    )
    event_rate <- event_rate_current
  }

  event_time <- unlist(purrr::map2(df_terminal$Status,
                                   df_terminal$Time,
                                   function(x, y) {
                                     if (x == 0) {
                                       return(0)
                                     } else {
                                       return(y)
                                     }
                                   }))
  upper_bound_censoring <- NULL
  n_times <- 0

  set.seed(NULL)
  while (is.null(upper_bound_censoring) && n_times < 100) {
    seed <- sample(1:100000, 1)
    upper_bound_censoring <- tryCatch(
      uniroot(function(x) {
        set.seed(seed)
        sum(runif(n_patients, min = 0, max = x) > event_time &
              df_terminal$Status == 1) / n_patients - event_rate
      }, interval = c(0, max(event_time) * 1000))$root,
      error = function(e) {
        NULL
      }
    )
    n_times <- n_times + 1
  }

  if(is.null(upper_bound_censoring))
    stop("Try to decrease the desired `event_rate` and rerun")

  if (verbose)
    print(paste0("Upper bound of random censoring is ", upper_bound_censoring))

  n_recurrent_per_patient <- as.vector(table(data_to_convert$Id))

  set.seed(seed)
  censor_time_per_patient <-
    runif(n_patients, min = 0, max = upper_bound_censoring)
  # print(sum(censor_time_per_patient > event_time & df_terminal$Status == 1)/n_patients)
  tt <- data_to_convert %>%
    dplyr::mutate(
      `censoring_time` = rep(censor_time_per_patient, times = n_recurrent_per_patient),
      `censored` = `Time` > `censoring_time`,
      `newTime` = pmin(`Time`, `censoring_time`)
    ) %>%
    dplyr::group_by(Id) %>%
    dplyr::mutate(nthCensor = cumsum(censored)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(nthCensor == 1 | censored == FALSE)

  Status <- purrr::map2(tt$Status,
                        tt$censored,
                        function(x, y) {
                          if (y == TRUE) {
                            return(0)
                          } else {
                            return(x)
                          }
                        })

  df_out <- data.frame(
    Id = tt$Id,
    X = tt$newTime,
    Status = unlist(Status),
    tt %>% dplyr::select(
      !c(
        "Id",
        "Time",
        "Status",
        "censored",
        "newTime",
        "nthCensor",
        "censoring_time"
      )
    )
  )

  return(df_out)
}
