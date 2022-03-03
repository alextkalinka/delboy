#' mark_up_FPs
#'
#' Combines hits from original and validation data for comparison purposes.
#'
#' @param res Output from `delboy::get_crispr_gene_level_hits`.
#' @param filter_column A character string naming the metric against which FPs will be defined.
#' @param filter_thr A numeric threshold value to be applied to `filter_column`.
#'
#' @return A data frame.
#' @export
#' @importFrom dplyr %>% select filter mutate
mark_up_FPs <- function(res, filter_column, filter_thr){
  tryCatch({
    
    if(filter_column == "mean_log2FoldChange"){
      res$hmp_gene_pos %<>% 
        dplyr::filter((!! rlang::sym(filter_column) > filter_thr | hit_type == "False_Negative"))
    }else if(filter_column == "sd_log2FoldChange"){
      res$hmp_gene_pos %<>% 
        dplyr::filter((!! rlang::sym(filter_column) < filter_thr | hit_type == "False_Negative"))
    }else{
      stop(paste("expecting 'filter_column' to be one of: 'mean_log2FoldChange' or 'sd_log2FoldChange', got:",filter_column))
    }
  },
  error = function(e) stop(paste("unable to mark up FPs:",e))
  )
}
