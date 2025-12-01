#' Fit a random forest for composite endpoints using CPIU-wide data
#'
#' @description
#' `Rforce()` fits a random forest model for composite endpoints—including
#' recurrent non-fatal events and a terminal event—using the CPIU (composite
#' pseudo–at-risk interval)–wide representation as the working data structure.
#'
#' The function can either:
#' \itemize{
#'   \item take raw patient-level long-formatted \code{data} and a
#'         \code{Surv(id, time, status) ~ covariates} \code{formula}, internally
#'         construct the CPIU-wide design matrix and auxiliary features, and then
#'         fit the forest; or
#'   \item directly accept a pre-computed \code{cpius} object (typically created
#'         by \code{\link{patients_to_cpius}} and \code{\link{cpius_to_dummy}})
#'         that already contains the CPIU-wide matrices and meta-data.
#' }
#'
#' @param data Optional \code{data.frame} containing patient-level data.
#'   Must include at least the patient ID, event time, and event status columns
#'   referenced on the left-hand side of \code{formula}. Ignored if \code{cpius}
#'   is supplied.
#' @param formula A \code{\link[survival]{Surv}} formula of the form
#'   \code{Surv(id, time, status) ~ covariate_1 + covariate_2 + \dots}.
#'   The left-hand side must contain three variables in order: patient ID,
#'   follow-up time, and event status. If a character string is provided, it is
#'   internally converted via \code{as.formula}. Used only when \code{cpius} is
#'   \code{NULL}.
#' @param n_intervals Integer; number of follow-up intervals used to build the
#'   CPIU-wide representation when \code{cpius} is not supplied. Defaults to
#'   \code{10}. The default break points are based on empirical quantiles of the
#'   follow-up time.
#' @param weights_by_status Optional numeric vector of non-negative weights, one
#'   per unique event status in the original patient-level data (including
#'   censoring). If supplied, its length must match the number of unique values
#'   in the status variable. If \code{NULL}, the default is weight 0 for
#'   censoring (status = 0) and weight 1 for all event types.
#' @param cpius Optional list containing pre-computed CPIU-wide data, typically
#'   the output of \code{\link{patients_to_cpius}} (optionally processed by
#'   \code{\link{cpius_to_dummy}}). It is expected to contain at least:
#'   \itemize{
#'     \item \code{design_matrix_Y}: numeric matrix of CPIU-wide outcomes and covariates;
#'     \item \code{auxiliary_features}: numeric matrix of auxiliary CPIU-level
#'           features (e.g., pseudo–at-risk durations, event indicators, weights)
#'           aligned with \code{design_matrix_Y};
#'     \item \code{variableIds}: integer vector of group IDs for dummy-encoded
#'           covariates;
#'     \item \code{unitsOfCPIUs}: numeric vector giving interval widths;
#'     \item \code{isDummy}: logical flag indicating whether covariates are
#'           already dummy-encoded.
#'   }
#'   If \code{cpius} is supplied, the internal data-to-CPIU conversion step is
#'   skipped.
#' @param split_rule Character string specifying the splitting rule:
#'   \itemize{
#'     \item \code{"Rforce-QLR"} — quasi-likelihood ratio–based splitting;
#'     \item \code{"Rforce-GEE"} — GEE-based splitting without interaction;
#'     \item \code{"Rforce-GEEIntr"} — GEE-based splitting with interaction;
#'     \item \code{"RF-SLAM"} — RF-SLAM–style splitting without pseudo–at-risk
#'           durations.
#'   }
#'   Partial matching is not supported; the value must match one of these
#'   strings exactly.
#' @param n_trees Integer; number of trees in the forest. Default is \code{300}.
#' @param mtry Integer; number of candidate variables drawn at random at each
#'   split. If \code{NULL}, a default is chosen internally based on the number
#'   of covariates.
#' @param n_splits Integer; number of candidate split points considered per
#'   variable at each node. If \code{NULL}, a default is chosen internally.
#' @param min_gain Numeric; minimum improvement in the split criterion required
#'   to accept a split. For \code{"Rforce-QLR"} and \code{"RF-SLAM"} this is
#'   interpreted as the percentage of quasi-likelihood increase (e.g. default
#'   \code{0.01} corresponds to 1\%). For \code{"Rforce-GEE"} and
#'   \code{"Rforce-GEEIntr"}, it is treated as a threshold on an adjusted
#'   \eqn{p}-value (e.g. default 0.2).
#' @param min_node_size Integer; minimum number of patients required in a
#'   terminal node (or to attempt further splitting). If \code{NULL}, a default
#'   is chosen internally based on the sample size.
#' @param max_depth Integer; maximum depth of each tree. Use \code{NULL} or a
#'   large value to allow essentially full growth, subject to \code{min_node_size}
#'   and \code{min_gain}. The default is \code{8}.
#' @param seed Integer random seed used to initialize the C back end for
#'   reproducible forests. Defaults to \code{926}.
#'
#' @details
#' When \code{cpius} is \code{NULL}, \code{Rforce()}:
#' \enumerate{
#'   \item checks that \code{data} and \code{formula} are supplied and that the
#'         left-hand side of \code{formula} is of the form
#'         \code{Surv(id, time, status)};
#'   \item validates that all variables appearing on the left- and right-hand
#'         sides of \code{formula} are present in \code{data};
#'   \item computes CPIU interval widths from the empirical distribution of the
#'         follow-up time using \code{n_intervals} empirical quantiles;
#'   \item constructs a long-format working dataset combining covariates and the
#'         \code{Id}, \code{X}, and \code{Status} response variables;
#'   \item calls \code{\link{patients_to_cpius}} to build a CPIU-wide design
#'         matrix and auxiliary features, with pseudo–at-risk durations used for
#'         \code{"Rforce-QLR"}, \code{"Rforce-GEE"}, and
#'         \code{"Rforce-GEEIntr"}, and disabled for \code{"RF-SLAM"};
#'   \item applies \code{\link{cpius_to_dummy}} to dummy-encode factor and
#'         character covariates if they are not already dummy-encoded;
#'   \item passes the resulting CPIU-wide matrices, variable grouping, and
#'         interval widths to the C routine \code{R_Rforce}, together with the
#'         specified splitting rule and hyperparameters.
#' }
#'
#' The fitted forest is returned as an S3 object of class \code{"Rforce"} with
#' components including the CPIU-wide data structure, forest structure, and
#' out-of-bag predictions, along with various meta-data to facilitate downstream
#' methods such as \code{print.Rforce()}, \code{summary.Rforce()}, and
#' \code{predict.Rforce()}.
#'
#' @return
#' An object of class \code{"Rforce"}, a list with elements including (but not
#' limited to):
#' \itemize{
#'   \item \code{cpius} — the CPIU-wide list containing \code{design_matrix_Y},
#'         \code{auxiliary_features}, \code{variableIds}, \code{unitsOfCPIUs},
#'         and related meta-data;
#'   \item \code{bagMatrix} — bootstrap sample indicators for each tree;
#'   \item \code{treePhi} — estimates of the mean structure for each interval in
#'         each tree;
#'   \item \code{forestMatrix} — a compact summary of the forest structure
#'         in matrix form (not intended to be interpreted directly by users);
#'   \item \code{vimpStat} — forest structure and variable importance summaries
#'         returned from the C back end;
#'   \item \code{predicted}, \code{oobPredicted} — fitted and out-of-bag
#'         predictions on the CPIU scale;
#'   \item hyperparameters such as \code{nTrees}, \code{maxDepth},
#'         \code{minNodeSize}, \code{minGain}, \code{mtry}, \code{nsplits},
#'         and \code{seed};
#'   \item \code{_external_forest_C_Ptr} — an external pointer storing the
#'         forest from the Rforce C API. It is valid only for the current R
#'         session unless serialized and restored via \code{\link{saveRforce}}
#'         and \code{\link{loadRforce}}.
#' }
#'
#' @seealso
#' \code{\link{compo_sim}},
#' \code{\link{random_censoring}},
#' \code{\link{patients_to_cpius}},
#' \code{\link{cpius_to_dummy}},
#' \code{\link[survival]{Surv}},
#' \code{\link{predict.Rforce}},
#' \code{\link{summary.Rforce}},
#' \code{\link{vimp.Rforce}},
#' \code{\link{saveRforce}},
#' \code{\link{loadRforce}}.
#'
#' @examples
#' ## Example using simulated data
#' library(Rforce)
#'
#' set.seed(926)
#'
#' data_list <- compo_sim()
#'
#' df_censored <- random_censoring(
#'   data       = data_list$dataset,
#'   event_rate = 0.5
#' )
#'
#' fit <- Rforce(
#'   data         = df_censored,
#'   formula      = Surv(Id, X, Status) ~ .,
#'   n_intervals  = 5,
#'   n_trees      = 50,
#'   max_depth    = 10,
#'   min_node_size = 20,
#'   min_gain     = 0,
#'   mtry         = 2,
#'   n_splits     = 5
#' )
#'
#' fit
#'
#' @export

Rforce <- function(
  data = NULL,
  formula = NULL, # e.g. Surv(Id, X, Status) ~ covariate1 + covariate2
  n_intervals = 10,
  weights_by_status = NULL,
  cpius = NULL,
  split_rule = c("Rforce-QLR", "Rforce-GEE", "Rforce-GEEIntr", "RF-SLAM"),
  n_trees = 300,
  mtry = NULL,
  n_splits = NULL,
  min_gain = NULL,
  min_node_size = NULL,
  max_depth = 8,
  seed = 926
) {
  split_rule <- match.arg(split_rule)

  ## choose splitting rule
  split_rule_index <- switch(
    split_rule,
    "Rforce-QLR" = 0,
    "Rforce-GEE" = 1,
    "Rforce-GEEIntr" = 1,
    "RF-SLAM" = 2,
    -1
  )

  if (split_rule_index == -1) {
    stop(
      'Invalid split_rule, must choose from c("Rforce-QLR", "Rforce-GEE", "Rforce-GEEIntr", "RF-SLAM")'
    )
  }

  gee_interaction <- switch(split_rule, "Rforce-GEEIntr" = 1, 0)

  # validate/normalize 'formula' if provided
  if (!is.null(formula)) {
    if (inherits(formula, "formula")) {
      # already a formula; OK
    } else if (is.character(formula) && length(formula) == 1) {
      tmp_formula <- try(as.formula(formula), silent = TRUE)
      if (inherits(tmp_formula, "try-error")) {
        stop(
          "'formula' is a character string but could not be parsed as a formula."
        )
      }
      formula <- tmp_formula
    } else {
      stop(
        "'formula' must be a formula object or a single string convertible to a formula."
      )
    }
  }

  ## If user provided raw data + formula, build design_matrix_Y and auxiliary_features here
  if (is.null(cpius)) {
    # require data, formula and n_intervals to be provided when not passing precomputed matrices
    if (is.null(data) || is.null(formula) || is.null(n_intervals)) {
      stop(
        "If design_matrix_Y or auxiliary_features are not provided you must supply 'data', 'formula' and 'n_intervals' (e.g. Surv(id, time, status) ~ covariates)."
      )
    }

    formula <- as.formula(formula)

    # Extract time and status variable names from left hand side of the formula
    lhs_vars <- all.vars(formula[[2]])
    rhs_vars <- all.vars(formula[[3]])
    if (length(lhs_vars) != 3) {
      stop(
        "Formula left hand side must be a Surv(id, x, status) with patient id, event time, and status variables, e.g. Surv(id, time, status) ~ covariate_1 + covariate_1."
      )
    }

    validate(data, required_cols = lhs_vars)
    validate(data, required_cols = rhs_vars)

    data_response <- data %>%
      dplyr::select(dplyr::all_of(lhs_vars))

    colnames(data_response) <- c("Id", "X", "Status")

    data_to_convert <- cbind.data.frame(
      data %>% dplyr::select(dplyr::all_of(rhs_vars)),
      data_response
    )

    # compute units_of_cpius if not supplied (default: deciles of follow-up time)
    probs <- seq(0, 1, length.out = n_intervals + 1)
    observed_times <- data_response %>% dplyr::pull(X)
    break_points <- stats::quantile(
      observed_times,
      probs = probs,
      na.rm = TRUE
    )
    break_points[1] <- 0.0 # ensure starting at 0
    break_points[n_intervals + 1] <- max(observed_times, na.rm = TRUE) * 1.001 # ensure ending at max time
    units_of_cpius <- diff(break_points) # compute widths of intervals

    list_status <- sort(unique(data_to_convert %>% dplyr::pull(Status)))
    if (!is.null(weights_by_status)) {
      assertthat::assert_that(
        length(list_status) == length(weights_by_status),
        msg = "The length of weights_by_status must be same with length of unique statuses!"
      )
    } else {
      # default weights: 0 for censoring, 1 for all event types
      weights_by_status <- rep(1, length(list_status))
      weights_by_status[which(list_status == 0)] <- 0
    }

    message("Converting data to design matrix and auxiliary features...")

    # Call patients_to_cpius to produce wide-format design matrix and auxiliary features
    cpius <- patients_to_cpius(
      data_to_convert = data_to_convert,
      units_of_cpiu = units_of_cpius,
      weights_by_status = weights_by_status,
      pseudo_risk = ifelse(split_rule_index == 2, FALSE, TRUE),
      wide_format = TRUE
    )

    cpius <- cpius_to_dummy(cpius)
  } # end internal conversion

  validate(cpius)

  if (!cpius$isDummy) {
    message("Dummy-encoding factor/character covariates...")
    cpius <- cpius_to_dummy(cpius)
  }

  if (is.null(mtry)) {
    mtry <- floor(sqrt(length(unique(cpius$variableIds))))
  }

  if (is.null(n_splits)) {
    # 5 is okay for dummy-encoded data, internally only 1 split can be applied per dummy variable, handled in C code
    n_splits <- 5
  }

  if (is.null(min_node_size)) {
    min_node_size <- max(10, floor(0.01 * length(unique(cpius$design_matrix_Y[, 1]))))
  }

  if (is.null(min_gain)) {
    if (split_rule_index == 0 || split_rule_index == 2) {
      min_gain <- 0.01 # 1% for QLR and RF-SLAM
    } else {
      min_gain <- 0.2 # p-value threshold for GEE-based
    }
  }

  if (is.null(max_depth)) {
    max_depth <- 8 
  }

  if (split_rule_index == 2) {
    message("Fitting RF-SLAM forest...")
  } else if (split_rule_index == 1 && gee_interaction == 1) {
    message("Fitting Rforce-GEEIntr forest...")
  } else if (split_rule_index == 1 && gee_interaction == 0) {
    message("Fitting Rforce-GEE forest...")
  } else {
    message("Fitting Rforce-QLR forest...")
  }

  temp_C_return <- .Call(
    "R_Rforce",
    as.integer(split_rule_index),
    as.integer(gee_interaction),
    matrix(
      as.numeric(unlist(cpius$designMatrix_Y)),
      nrow = dim(cpius$designMatrix_Y)[1],
      ncol = dim(cpius$designMatrix_Y)[2]
    ),
    matrix(
      as.numeric(unlist(cpius$auxiliaryFeatures)),
      nrow = dim(cpius$auxiliaryFeatures)[1],
      ncol = dim(cpius$auxiliaryFeatures)[2]
    ),
    as.integer(cpius$variableIds),
    as.integer(length(unique(cpius$variableIds))),
    as.numeric(cpius$unitsOfCPIUs),
    as.integer(n_trees),
    as.integer(max_depth),
    as.integer(min_node_size),
    as.numeric(min_gain),
    as.integer(mtry),
    as.integer(n_splits),
    as.integer(seed)
  )

  # Extract struct values from the returned list
  rforce_obj <- list(
    cpius = cpius,
    forestMatrix = temp_C_return[["forestMatrix"]],
    predicted = temp_C_return[["predicted"]],
    oobPredicted = temp_C_return[["oobPredicted"]],
    vimpStat = temp_C_return[["vimpStat"]],
    bagMatrix = temp_C_return[["bagMatrix"]],
    treePhi = temp_C_return[["treePhi"]],
    nTrees = n_trees,
    maxDepth = max_depth,
    minNodeSize = min_node_size,
    minGain = min_gain,
    mtry = mtry,
    nsplits = n_splits,
    seed = seed,
    `_external_forest_C_Ptr` = temp_C_return[["_external_forest_C_Ptr"]]
  )

  class(rforce_obj) <- "Rforce"
  return(rforce_obj)
}
