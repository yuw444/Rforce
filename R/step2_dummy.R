#' Dummy-encode factor/character covariates in CPIU object
#'
#' This function takes a CPIU object and dummy-encodes any factor or character covariates
#' in the design matrix for the outcome variable Y. It removes the reference level
#' for each factor to avoid multicollinearity.
#' @importFrom magrittr %>%
#' @importFrom tidyr starts_with
#' @importFrom dplyr select all_of
#' @importFrom fastDummies dummy_cols
#' @param object A CPIU object.
#' @param cols_to_dummy Optional vector of column names to dummy-encode. If NULL, all factor/character columns will be dummy-encoded.
#' @return A CPIU object with dummy-encoded covariates in the `designMatrixY`.

cpius_to_dummy <- function(object, cols_to_dummy = NULL) {
    validate(object)

    df_designMatrix <- object$designMatrix_Y %>%
        dplyr::select(-tidyr::starts_with("nEvents"))

    df_Y <- object$designMatrix_Y %>%
        dplyr::select(tidyr::starts_with("nEvents"))

    if (!is.null(cols_to_dummy)) {
        validate(
            df_designMatrix,
            required_cols = cols_to_dummy
        )

        df_designMatrix <- df_designMatrix %>%
            dplyr::select(all_of(cols_to_dummy))
    }

    if (
        any(unlist(lapply(df_designMatrix, function(col) {
            return(is.character(col) || is.factor(col))
        })))
    ) {
        message("Dummy-encoding factor/character covariates...")
        # dummy-encode factors (keep same approach as tests: remove reference column)
        covariate_df <- fastDummies::dummy_cols(
            df_designMatrix,
            remove_first_dummy = TRUE, # drop reference level per factor
            ignore_na = TRUE,
            remove_selected_columns = TRUE # drop original factor/char columns
        )
        # ensure no rownames and keep as data.frame
        covariate_df <- as.data.frame(covariate_df, stringsAsFactors = FALSE)
        variable_Ids_names <- colnames(covariate_df)
        # remove internal suffixes that to_dummy may add like ".x_1"
        variable_Ids_names <- gsub("\\_\\d+", "", variable_Ids_names)
        unique_vars <- unique(variable_Ids_names)
        variable_Ids_int <- match(variable_Ids_names, unique_vars) - 1
    } else {
        covariate_df <- df_designMatrix
        variable_Ids_int <- as.integer(variable_Ids)
    }

    object$designMatrix_Y <- cbind.data.frame(
        covariate_df,
        df_Y
    )

    object$variableUsed <- colnames(df_designMatrix)
    object$variableIds <- variable_Ids_int
    object$isDummy <- TRUE

    return(object)
}
