#' @title censoring the given dataset manually
#'
#' @import dplyr
#' @description censoring the given dataset based on the given quantile of observed time
#' @param data_to_convert the dataset to be converted
#' @param censoring_quantile the quantile of observed time to be censored
#' @return the censored dataset
#' @export
#' @examples
#' list_data_to_convert <- compo_sim( )
#' df_converted <- manual_censoring(
#'     list_data_to_convert$dataset,
#'     0.9
#' )
#' str(df_converted)
#'
manual_censoring <- function(
    data_to_convert,
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
    dplyr::mutate(
      cenosored = Time > censoring_time,
      newTime = pmin(Time, censoring_time)
    ) %>%
    dplyr::group_by(Id) %>%
    dplyr::mutate(
      nthCensor = cumsum(cenosored)
    ) %>%
    ungroup() %>%
    dplyr::filter(nthCensor == 1 | cenosored == FALSE)

  Status <- purrr::map2(
    tt$Status,
    tt$cenosored,
    function(x, y) {
      if (y == TRUE) {
        return(0)
      } else {
        return(x)
      }
    }
  )

  df_out <- data.frame(
    Id = tt$Id,
    X = tt$newTime,
    Status = unlist(Status),
    tt %>% dplyr::select(starts_with("bin"), starts_with("con"))
  )

  return(df_out)
}

