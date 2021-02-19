#' calculate_bayes_factors_2dkcde
#'
#' Calculates log Bayes Factors using a 2d kernel cumulative distribution estimate for True and False Positives along the dimensions of log fold change (LFC) and expression level (EL): Pr(LFC, EL | TP)/Pr(LFC, EL | FP)).
#'
#' @param data A data frame containing data for which Bayes Factors will be calculated.
#' @param data_bs A data frame containing bootstrap-sampled data for which a 2D kernel cumulative distribution estimate will be calculated for the data points in `data`.
#' @param fc_column A character string naming the log fold change column.
#' @param expr_column A character string naming the expression level column.
#' @param group_column A character string naming the grouping column containing only `True_Positive` and `False_Positive` entries.
#'
#' @return A list containing the following elements:
#' * `kcde_TP` - an object of class `kcde` for `True_Positives`.
#' * `kcde_FP` - an object of class `kcde` for `False_Positives`.
#' * `cum_prob_BF` - a data frame containing cumulative probability estimates and Bayes Factors for the data in `data`.
#'
#' @md
#' @export
#' @importFrom rlang !! sym
#' @importFrom dplyr select filter mutate %>% rowwise ungroup
#' @importFrom ks kcde
calculate_bayes_factors_2dkcde <- function(data, data_bs, fc_column, expr_column, group_column){
  # Sanity checks.
  if(!group_column %in% colnames(data_bs)) stop(paste("unable to find",group_column,"in 'data'"))
  if(any(!data[,group_column] %in% c("True_Positive","False_Positive")))
    stop(paste("column",group_column,"should only contain 'True_Positive' or 'False_Positive' entries"))

  tryCatch({
    # 1. Prep data.
    data_query <- data %>%
      dplyr::select(!!rlang::sym(fc_column), !!rlang::sym(expr_column))

    data_bs.tp <- data_bs %>%
      dplyr::filter(!!rlang::sym(group_column) == "True_Positive") %>%
      dplyr::select(!!rlang::sym(fc_column), !!rlang::sym(expr_column))

    data_bs.fp <- data_bs %>%
      dplyr::filter(!!rlang::sym(group_column) == "False_Positive") %>%
      dplyr::select(!!rlang::sym(fc_column), !!rlang::sym(expr_column))

    # 2. Estimate 2d cumulative distribution functions for TPs and FPs.
    Fhat_tp <- ks::kcde(as.matrix(data_bs.tp))
    Fhat_fp <- ks::kcde(as.matrix(data_bs.fp))

    # 3. Predict 2d cumulative probabilities for query data.
    cp_tp <- predict(Fhat_tp, x = as.matrix(data))
    cp_fp <- predict(Fhat_fp, x = as.matrix(data))

    # 4. Calculate Bayes Factors.
    fc_median.tp <- median(data_bs.tp[,fc_column], na.rm = T)
    fc_median.fp <- median(data_bs.fp[,fc_column], na.rm = T)
    expr_median.tp <- median(data_bs.tp[,expr_column], na.rm = T)
    expr_median.fp <- median(data_bs.fp[,expr_column], na.rm = T)
    data_bf <- data %>%
      dplyr::mutate(kcum_prob_TP = cp_tp, kcum_prob_FP = cp_fp) %>%
      # Calculate 'lower.tail' cum_prob for data to the right of both distributions.
      dplyr::rowwise() %>%
      dplyr::mutate(kcum_prob_TP = ifelse((!!rlang::sym(fc_column) > fc_median.tp & !!rlang::sym(expr_column) > expr_median.tp),
                                   1 - kcum_prob_TP, kcum_prob_TP),
             kcum_prob_FP = ifelse((!!rlang::sym(fc_column) > fc_median.fp & !!rlang::sym(expr_column) > expr_median.fp),
                                   1 - kcum_prob_FP, kcum_prob_FP)) %>%
      dplyr::ungroup() %>%
      # Calculate Bayes factors.
      dplyr::mutate(log_bayes_factor = log(kcum_prob_TP/kcum_prob_FP))
  },
  error = function(e) stop(paste("unable to calculate bayes factors:",e))
  )

}
