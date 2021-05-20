#' apply_thresholds_val
#' 
#' Applies p-value and abs(LFC) thresholds to validation data.
#' 
#' @param data A data frame containing DEG results for multiple validation samples.
#' @param thresh An object produced by `delboy::find_pval_target_fdr`.
#' @return A modified data frame.
#' @export
#' @importFrom dplyr %>% mutate rowwise ungroup case_when
#' @importFrom magrittr %<>%
apply_thresholds_val <- function(data, thresh){
  tryCatch({
    pval_thresh <- thresh$pvalue_target_FDR
    lfc_thresh <- thresh$abs_Log2FoldChange.argmax_KS_dist
    data %<>%
      dplyr::rowwise() %>%
      dplyr::mutate(hit_type = dplyr::case_when()) %>%
      dplyr::ungroup()
  },
  error = function(e) stop(paste("unable to apply p-value and logFC thresholds to validation data:",e))
  )
}
