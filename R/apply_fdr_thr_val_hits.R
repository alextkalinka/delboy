#' apply_fdr_thr_val_hits
#'
#' Applies an FDR threshold (for a given metric) to control the FDR for validation hits.
#'
#' @param data A data frame as produced by `delboy::make_delboy_crispr_hit_comparison_table`.
#' @param filter_column A character string naming a column to be filtered.
#' @param filter_thr A numeric threshold for `filter_column`.
#' @param dir A character string identifying which end of logFC values we are focusing on: `neg` or `pos`.
#' @param input_type A character string naming the input type: `hits` (default), or `all`.
#'
#' @return A data frame.
#' @importFrom dplyr filter
#' @importFrom magrittr %<>%
#' @importFrom rlang sym !!
#' @export
apply_fdr_thr_val_hits <- function(data, filter_column, filter_thr, dir, input_type = "hits"){
  tryCatch({
    if(input_type == "hits"){
      if((filter_column == "mean_log2FoldChange" && dir == "greater")){
        data %<>% dplyr::filter((!! rlang::sym(filter_column) > filter_thr | hit_type == "False_Negative"))
      }else{
        data %<>% dplyr::filter((!! rlang::sym(filter_column) < filter_thr | hit_type == "False_Negative"))
      }
    }else if(input_type == "all"){
      if((filter_column == "mean_log2FoldChange" && dir == "greater")){
        data %<>% dplyr::filter(!(!! rlang::sym(filter_column) < filter_thr & significant_hit))
      }else{
        data %<>% dplyr::filter(!(!! rlang::sym(filter_column) > filter_thr & significant_hit))
      }
    }else{
      stop(paste("'input_type' must be one of 'hits' or 'all', got:",input_type))
    }
  },
  error = function(e) stop(paste("unable to apply FDR threshold to validation hits:",e))
  )
  return(data)
}
