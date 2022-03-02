#' calc_perf_stats
#'
#' Returns sensitivity and FDR estimates.
#'
#' @param data A data frame as produced by `delboy::make_delboy_crispr_hit_comparison_table`.
#'
#' @return A data frame.
#' @importFrom dplyr %>% filter summarise n
#' @export
calc_perf_stats <- function(data){
  tryCatch({
    ps <- data %>%
      dplyr::summarise(total = dplyr::n(),
                       tot_tp = sum(all_val_hits$hit_type == "True_Positive"),
                       tot_fp = sum(all_val_hits$hit_type == "False_Positive"),
                       tot_fn = sum(all_val_hits$hit_type == "False_Negative"),
                       Sensitivity = 100*tot_tp/(tot_tp + tot_fn),
                       FDR = 100*tot_fp/(tot_tp + tot_fp),
                       FDR = ifelse(is.nan(FDR),0,FDR),
                       Precision = 1-FDR)
  },
  error = function(e) stop(paste("unable to calculate performance statistics:",e))
  )
  return(ps)
}
