#' smooth_decision_boundary
#' 
#' Smoothes the fold-change expression decision boundary such that it is a convex, monotonically decreasing function of expression level.
#' 
#' @param db A data frame providing the TP-FP decision boundary for `log10_baseExpr` and `abs_log2FoldChange` columns.
#' @param entry_point Internal point to start smoothing in terms of `log10_baseExpr`. Defaults to 1.25.
#' @return A data frame of smoothed values
#' @export
#' @importFrom dplyr filter arrange desc
#' @importFrom stats smooth.spline predict
smooth_decision_boundary <- function(db, entry_point = 1.25){
  tryCatch({
    db_low <- db %>%
      dplyr::filter(log10_baseExpr <= entry_point) %>%
      dplyr::arrange(dplyr::desc(log10_baseExpr))
    db_high <- db %>%
      dplyr::filter(log10_baseExpr > entry_point) %>%
      dplyr::arrange(log10_baseExpr)
    
    # 1. High-end values.
    runn_delta <- 0
    for(i in 1:nrow(db_high)){
      if(i == 1) next
      runn_delta <- db_high$abs_log2FoldChange[i] - db_high$abs_log2FoldChange[i-1]
      # Anti-tonic.
      if(db_high$abs_log2FoldChange[i] > db_high$abs_log2FoldChange[i-1]){
        db_high$abs_log2FoldChange[i] <- db_high$abs_log2FoldChange[i-1]
      }else{
        if(i == 2) next
        # Convex.
        if(db_high$abs_log2FoldChange[i] - db_high$abs_log2FoldChange[i-1] < runn_delta){
          db_high$abs_log2FoldChange[i] <- db_high$abs_log2FoldChange[i-1] - runn_delta
        }
      }
    }
    db$abs_log2FoldChange[match(db_high$log10_baseExpr, db$log10_baseExpr)] <- db_high$abs_log2FoldChange
    
    # 2. Low-end values.
    if(max(db_low$abs_log2FoldChange) < 4){
      spl_fit <- stats::smooth.spline(db_low$abs_log2FoldChange, db_low$log10_baseExpr)
      db_low_expr <- stats::predict(spl_fit, seq(max(db_low$abs_log2FoldChange),4,0.02))
      db %<>%
        rbind(data.frame(log10_baseExpr = db_low_expr$y, 
                         abs_log2FoldChange = db_low_expr$x))
    }
    
    return(db)
  },
  error = function(e) stop(paste("unable to smooth decision boundary:",e))
  )
}
