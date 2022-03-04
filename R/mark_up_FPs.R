#' mark_up_FPs
#'
#' Combines hits from original and validation data for comparison purposes.
#'
#' @param res Output from `delboy::combine_harmonic_mean_pvals`.
#' @param filter_column A character string naming the metric against which FPs will be defined.
#' @param filter_thr A numeric threshold value to be applied to `filter_column`.
#' @param dir A character string identifying which end of logFC values we are focusing on: `less` or `greater`.
#'
#' @return A data frame.
#' @export
#' @importFrom dplyr mutate
#' @importFrom magrittr %<>%
#' @importFrom rlang sym !!
mark_up_FPs <- function(res, filter_column, filter_thr, dir){
  tryCatch({
    if((filter_column == "mean_log2FoldChange" && dir == "greater")){
      res %<>% 
        dplyr::mutate(Predicted_False_Positive = (!! rlang::sym(filter_column) < filter_thr & significant_hit))
    }else{
      res %<>% 
        dplyr::mutate(Predicted_False_Positive = (!! rlang::sym(filter_column) > filter_thr & significant_hit))
    }
  },
  error = function(e) stop(paste("unable to mark up FPs:",e))
  )
  return(res)
}
