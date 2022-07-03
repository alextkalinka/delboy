#' mark_up_FPs
#'
#' Combines hits from original and validation data for comparison purposes.
#'
#' @param res Output from `delboy::combine_harmonic_mean_pvals`.
#' @param dir A character string identifying which end of logFC values we are focusing on: `less` or `greater`.
#' @param mlfc_fp_thresh Numeric indicating the scaled log-fold change threshold below which a hit should be considered a False Positive (default = 0.3).
#'
#' @return A data frame.
#' @export
#' @importFrom dplyr mutate group_by ungroup
#' @importFrom magrittr %<>%
#' @importFrom stats median
mark_up_FPs <- function(res, dir, mlfc_fp_thresh = 0.3){
  tryCatch({
    if(dir == "pos"){
      res %<>% 
        dplyr::group_by(significant_hit) %>%
        dplyr::mutate(mean_log2FoldChange_hits.scaled = ifelse(significant_hit, mean_log2FoldChange/median(mean_log2FoldChange,na.rm=T),
                                                               NA),
                      Predicted_False_Positive = (mean_log2FoldChange_hits.scaled < mlfc_fp_thresh & significant_hit)) %>%
        dplyr::ungroup()
    }else{
      res %<>% 
        dplyr::group_by(significant_hit) %>%
        dplyr::mutate(mean_log2FoldChange_hits.scaled = ifelse(significant_hit, -mean_log2FoldChange/median(-mean_log2FoldChange,na.rm=T),
                                                               NA),
                      Predicted_False_Positive = (mean_log2FoldChange_hits.scaled < mlfc_fp_thresh & significant_hit)) %>%
        dplyr::ungroup()
    }
  },
  error = function(e) stop(paste("unable to mark up FPs:",e))
  )
  return(res)
}
