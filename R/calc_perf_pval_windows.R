#' calc_perf_pval_windows
#' 
#' Calculates sensitivity and FDR across a progressively larger window of DEG p-values.
#' 
#' @param data A data frame of DEG p-values.
#' @param pval_column A character string naming the p-value column.
#' @param gene_column A character string naming the gene column.
#' @param deg_genes A character vector naming DEGs.
#' @return A data frame of performance estimates across the p-val thresholds.
#' @export
#' @importFrom dplyr %>% arrange mutate
#' @importFrom magrittr %<>%
#' @importFrom rlang !! sym
calc_perf_pval_windows <- function(data, pval_column, gene_column, deg_genes){
  tryCatch({
    data %<>%
      dplyr::arrange(!! rlang::sym(pval_column))
    num_tp <- length(deg_genes)
    sens <- fdr <- NULL
    for(i in 1:nrow(data)){
      tp <- as.vector(data[1:i,gene_column]) %in% deg_genes
      sens <- append(sens, 100*sum(tp)/num_tp)
      fdr <- append(fdr, 100*sum(!tp)/i)
    }
    data %<>%
      dplyr::mutate(Sensitivity.percent = sens, FDR.percent = fdr)
    return(data)
  },
  error = function(e) stop(paste("unable to calculate performance across p-value windows:",e))
  )
}
