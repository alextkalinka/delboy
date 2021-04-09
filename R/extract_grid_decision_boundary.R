#' extract_grid_decision_boundary
#'
#' Reduces a grid of class predictions to a decision boundary.
#'
#' @param grid A data frame with the following columns: `abs_log2FoldChange`, `log10_baseExpr`, and `Predicted_FP`.
#'
#' @return A data frame.
#' @export
#' @importFrom dplyr %>% group_by summarise ungroup
extract_grid_decision_boundary <- function(grid){
  tryCatch({
    db <- grid %>%
      dplyr::group_by(log10_baseExpr) %>%
      dplyr::summarise(abs_log2FoldChange.DB = ifelse(any(Predicted_FP == 1),
                                                      max(abs_log2FoldChange[Predicted_FP == 1]),
                                                      0)) %>%
      dplyr::ungroup()
  },
  error = function(e) stop(paste("unable to extract grid decision boundary:",e))
  )
}
