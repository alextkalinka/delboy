# Helper functions.
.predict_FP <- function(expr, log_fc, db){
  db_fc <- (db %>%
    dplyr::filter((abs(log10_baseExpr - expr)) == min(abs(log10_baseExpr - expr))[1]))$abs_log2FoldChange[1]
  return(as.integer(ifelse(db_fc > log_fc, 1, 0)))
}


#' predict_FP_delboy
#' 
#' Predict false positives in expression data.
#' 
#' @param data A data frame with the following columns: `log10_baseExpr` and `abs_log2FoldChange`.
#' @param db A data frame providing a decision boundary with the following columns: `log10_baseExpr` and `abs_log2FoldChange`.
#' 
#' @return A data frame.
#' @export
#' @importFrom dplyr %>% rowwise mutate ungroup filter
predict_FP_delboy <- function(data, db){
  tryCatch({
    fp <- data %>%
      dplyr::rowwise() %>%
      dplyr::mutate(pred_FP = .predict_FP(log10_baseExpr, abs_log2FoldChange, db)) %>%
      dplyr::ungroup()
    return(fp)
  },
  error = function(e) stop(paste("unable to predict FPs in expression data:",e))
  )
}
