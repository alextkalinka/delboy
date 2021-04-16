#' smooth_decision_boundary
#' 
#' Smoothes the fold-change expression decision boundary such that it is a convex, monotonically decreasing function of expression level.
#' 
#' @param db A data frame providing the TP-FP decision boundary for `log10_baseExpr` and `abs_log2FoldChange` columns.
#' @param min_log10_exp The minimum log10 expression level for the data.
#' @param entry_point Internal point to start smoothing in terms of `log10_baseExpr`. Defaults to 1.25.
#' @return A data frame of smoothed values
#' @export
#' @importFrom dplyr filter arrange desc
smooth_decision_boundary <- function(db, min_log10_exp, entry_point = 1.25){
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
    diffs <- NULL
    for(i in 1:nrow(db_low)){
      if(i == 1) next
      diffs <- append(diffs, (db_low$abs_log2FoldChange[i] - db_low$abs_log2FoldChange[i-1])/(db_low$log10_baseExpr[i] - db_low$log10_baseExpr[i-1]))
    }
    y_last <- -1*tail(diffs,1)*(tail(db_low$log10_baseExpr,1) + tail(db_low$abs_log2FoldChange,1))
    db <- rbind(db, data.frame(log10_baseExpr = 0, abs_log2FoldChange = y_last))
    
    return(db)
  },
  error = function(e) stop(paste("unable to smooth decision boundary:",e))
  )
}
