#' smooth_decision_boundary
#' 
#' Smoothes the fold-change expression decision boundary such that it is a convex, monotonically decreasing function of expression level.
#' 
#' @param db A data frame providing the TP-FP decision boundary for `log10_baseExpr` and `abs_log2FoldChange` columns.
#' @param entry_point Internal point to start smoothing in terms of `log10_baseExpr`. Defaults to 1.25.
#' @return A data frame of smoothed values
#' @export
#' @importFrom dplyr filter arrange desc
smooth_decision_boundary <- function(db, entry_point = 1.25){
  tryCatch({
    db_low <- db %>%
      dplyr::filter(log10_baseExpr <= entry_point) %>%
      dplyr::arrange(dplyr::desc(log10_baseExpr))
    db_high <- db %>%
      dplyr::filter(log10_baseExpr > entry_point) %>%
      dplyr::arrange(log10_baseExpr)
    
    # 1. High end values.
    runn_delta <- 0
    for(i in 1:nrow(db_high)){
      if(i == 1) next
      if(db_high$abs_log2FoldChange[i] > db_high$abs_log2FoldChange[i-1]){
        db_high$abs_log2FoldChange[i] <- db_high$abs_log2FoldChange[i-1]
      }else{
        if(i == 2) next
        if(db_high$abs_log2FoldChange[i] - db_high$abs_log2FoldChange[i-1] < runn_delta){
          
        }
      }
    } 
  },
  error = function(e) stop(paste("unable to smooth decision boundary:",e))
  )
}
