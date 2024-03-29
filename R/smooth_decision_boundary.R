# Helper functions.

# Anti-tonic, convex enforcing.
.anti_tonic_convex <- function(data, y){
  # data arranged by x (increasing).
  runn_delta <- 0
  for(i in 1:nrow(data)){
    if(i == 1) next
    runn_delta <- data[i,y] - data[i-1,y]
    # Anti-tonic.
    if(data[i,y] > data[i-1,y]){
      data[i,y] <- data[i-1,y]
    }else{
      if(i == 2) next
      # Convex.
      if(data[i,y] - data[i-1,y] < runn_delta){
        data[i,y] <- data[i-1,y] - runn_delta
      }
    }
  }
  return(data)
}


#' smooth_decision_boundary
#' 
#' Smoothes the fold-change expression decision boundary such that it is a convex, monotonically decreasing function of expression level.
#' 
#' @param db A data frame providing the TP-FP decision boundary for `log10_baseExpr` and `abs_log2FoldChange` columns.
#' @param entry_point Internal point to start smoothing in terms of `log10_baseExpr`. Defaults to 1.25.
#' @return A data frame of smoothed values
#' @export
#' @importFrom dplyr filter arrange desc
#' @importFrom stats smooth.spline predict IQR
#' @importFrom magrittr %<>%
smooth_decision_boundary <- function(db, entry_point = 1.25){
  tryCatch({
    # Check if entry point makes sense.
    if(sum(db$log10_baseExpr <= entry_point) <= 10){
      entry_point <- median(db$log10_baseExpr, na.rm = T)
      .db_message(paste("using DB smoothing entry point of",entry_point),"blue")
    }
    
    db_low <- db %>%
      dplyr::filter(log10_baseExpr <= entry_point) %>%
      dplyr::arrange(dplyr::desc(log10_baseExpr))
    db_high <- db %>%
      dplyr::filter(log10_baseExpr > entry_point) %>%
      dplyr::arrange(log10_baseExpr)
    
    # 1. High-end values.
    db_high <- .anti_tonic_convex(db_high, "abs_log2FoldChange")
    
    db$abs_log2FoldChange[match(db_high$log10_baseExpr, db$log10_baseExpr)] <- db_high$abs_log2FoldChange
    
    # Make sure DB extends across range of interest.
    if(max(db$log10_baseExpr, na.rm=T) < 6){
      db<- rbind(db, data.frame(log10_baseExpr = 6,
                                abs_log2FoldChange = min(db$abs_log2FoldChange, na.rm=T)))
    }
    
    # 2. Low-end values.
    if(max(db_low$abs_log2FoldChange) < 4){
      iqr_fc <- stats::IQR(db_low$abs_log2FoldChange)
      if(iqr_fc <= 0 || is.infinite(iqr_fc) || is.na(iqr_fc) || is.null(iqr_fc))
        stop(paste("IQR of abs_log2FoldChange:",iqr_fc))
      spl_fit <- stats::smooth.spline(db_low$abs_log2FoldChange, db_low$log10_baseExpr)
      db_low_expr <- stats::predict(spl_fit, seq(max(db_low$abs_log2FoldChange),4,0.02))
      db_low <- rbind(db_low,
                      data.frame(log10_baseExpr = db_low_expr$y, 
                           abs_log2FoldChange = db_low_expr$x)) %>%
        dplyr::arrange(abs_log2FoldChange)

      db_low <- .anti_tonic_convex(db_low, "log10_baseExpr")
      
      db %<>%
        rbind(db_low)
    }
    
    return(db)
  },
  error = function(e) stop(paste("unable to smooth decision boundary:",e))
  )
}
