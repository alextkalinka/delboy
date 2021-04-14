#' extract_grid_decision_boundary
#'
#' Reduces a grid of class predictions to a decision boundary.
#'
#' @param grid A data frame with the following columns: `abs_log2FoldChange`, `log10_baseExpr`, and `Predicted_FP`.
#'
#' @return A data frame.
#' @export
#' @importFrom dplyr %>% group_by summarise ungroup filter select
extract_grid_decision_boundary <- function(grid){
  tryCatch({
    db <- grid %>%
      dplyr::group_by(log10_baseExpr) %>%
      dplyr::summarise(max_fp_fc = ifelse(any(Predicted_FP == 1),
                                                      max(abs_log2FoldChange[Predicted_FP == 1]),
                                                      0),
                       min_tp_fc = ifelse(any(Predicted_FP == 0),
                                          min(abs_log2FoldChange[Predicted_FP == 0]),
                                          0),
                       abs_log2FoldChange = ifelse(max_fp_fc > min_tp_fc, min_tp_fc, max_fp_fc)) %>%
      dplyr::filter(abs_log2FoldChange != 0) %>%
      dplyr::ungroup() %>%
      dplyr::select(-max_fp_fc, -min_tp_fc)
    return(db)
  },
  error = function(e) stop(paste("unable to extract grid decision boundary:",e))
  )
}
