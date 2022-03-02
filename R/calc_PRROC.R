#' calc_PRROC
#'
#' Returns AUPrRc, precision, and recall values for a data frame containing TP and a score metric.
#'
#' @param data A data frame containing a 'TP' logical column indicating which genes are true positives.
#' @param score_col A character string naming a column containing a score metric.
#' @param group_col A character string naming the grouping column of 'data'. If `NULL` then no group, defaults to `NULL`.
#'
#' @return A data frame containing the following columns: `Precision`, `Recall`, `Sensitivity_FDR_10pct`, `Sensitivity_FDR_10pct`, and `AUPrRc`.
#' @importFrom dplyr mutate filter select rename
#' @importFrom PRROC pr.curve
#' @importFrom rlang sym
#' @export
calc_PRROC <- function(data, score_col, group_col = NULL){
  if(!"TP" %in% colnames(data))
    stop("get_PRROC: expecting to find a 'TP' column in 'data'")
  if(!score_col %in% colnames(data))
    stop(paste("get_PRROC: expecting to find a",score_col,"column in 'data'"))
  
  tryCatch({
    sc <- rlang::sym(score_col)
    pos_scores <- c((data %>%
                       as.data.frame() %>%
                       dplyr::filter(TP) %>%
                       dplyr::select(!!sc))[,1])
    neg_scores <- c((data %>%
                       as.data.frame() %>%
                       dplyr::filter(!TP) %>%
                       dplyr::select(!!sc))[,1])
    pr <- PRROC::pr.curve(pos_scores, neg_scores, curve = T)
    prc <- as.data.frame(pr$curve) %>%
      dplyr::rename(Recall = V1, Precision = V2, !!sc := V3) %>%
      dplyr::mutate(AUPrRc = pr$auc.davis.goadrich,
                    Sensitivity_FDR_10pct = 100*Recall[which(abs(Precision-0.9)==min(abs(Precision-0.9)))[1]],
                    Sensitivity_FDR_5pct = 100*Recall[which(abs(Precision-0.95)==min(abs(Precision-0.95)))[1]])
    if(!is.null(group_col))
      prc %<>% dplyr::mutate(!!rlang::sym(group_col) := unlist((data %>% dplyr::select(!!rlang::sym(group_col)))[,1])[1])
  },
  error = function(e) stop(paste("get_PRROC: unable to calculate precision-recall curve:",e))
  )
  return(prc)
}
