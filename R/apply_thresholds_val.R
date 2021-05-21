#' apply_thresholds_val
#' 
#' Applies p-value and abs(LFC) thresholds to validation data.
#' 
#' @param data A data frame containing DEG results for multiple validation samples.
#' @param thresh An object produced by `delboy::find_pval_target_fdr`.
#' @return A modified data frame.
#' @export
#' @importFrom dplyr %>% mutate rowwise ungroup case_when filter
#' @importFrom magrittr %<>%
apply_thresholds_val <- function(data, thresh){
  tryCatch({
    pval_thresh <- thresh$pvalue_target_FDR
    data %<>%
      # Those with signal plus those defined as hits by the pvalue threshold.
      dplyr::filter(signal | pvalue < pval_thresh) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(hit_type = dplyr::case_when((!signal & pvalue < pval_thresh) ~ "False_Positive",
                                                (signal & pvalue >= pval_thresh) ~ "False_Negative",
                                                (signal & pvalue < pval_thresh) ~ "True_Positive",
                                                TRUE ~ "Other")) %>%
      dplyr::filter(hit_type != "Other") %>%
      dplyr::ungroup()
    return(data)
  },
  error = function(e) stop(paste("unable to apply p-value and logFC thresholds to validation data:",e))
  )
}
