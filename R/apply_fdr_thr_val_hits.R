#' apply_fdr_thr_val_hits
#'
#' Applies an FDR threshold (for a given metric) to control the FDR for validation hits.
#'
#' @param data A data frame as produced by `delboy::make_delboy_crispr_hit_comparison_table`.
#' @param filter_column A character string naming a column to be filtered.
#' @param filter_thr A numeric threshold for `filter_column`.
#'
#' @return A data frame.
#' @importFrom dplyr filter
#' @importFrom magrittr %<>%
#' @importFrom rlang sym !!
#' @export
apply_fdr_thr_val_hits <- function(data, filter_column, filter_thr){
  tryCatch({
    if(filter_column == "mean_log2FoldChange"){
      data %<>% dplyr::filter((!! rlang::sym(filter_column) > filter_thr | hit_type == "False_Negative"))
    }else if(filter_column == "sd_log2FoldChange"){
      data %<>% dplyr::filter((!! rlang::sym(filter_column) < filter_thr | hit_type == "False_Negative"))
    }else{
      stop(paste("expecting 'filter_column' to be one of: 'mean_log2FoldChange' or 'sd_log2FoldChange', got:",filter_column))
    }
  },
  error = function(e) stop(paste("unable to apply FDR threshold to validation hits:",e))
  )
  return(data)
}
