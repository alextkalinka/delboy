#' calculate_perf_stats
#' 
#' Return sensitivity and precision estimates for delboy validation data as produced by `delboy::perf_stats_rnaseq`.
#' 
#' @param perf_stats A data frame of performance stats for individual validation samples.
#' 
#' @return A data frame.
#' @export
#' @importFrom dplyr %>% group_by summarise ungroup
#' @importFrom magrittr %<>%
calculate_perf_stats <- function(perf_stats){
  tryCatch({
    perf_stats %<>%
      dplyr::group_by(Algorithm) %>%
      dplyr::summarise(Total_diff_exp = sum(Total_diff_exp),
                       Num_true_calls = sum(Num_true_calls), Num_false_calls = sum(Num_false_calls),
                       Sensitivity.percent = 100*Num_true_calls/Total_diff_exp,
                       Precision.percent = 100*Num_true_calls/(Num_true_calls + Num_false_calls),
                       FDR.percent = 100 - Precision.percent) %>%
      dplyr::ungroup()
  },
  error = function(e) stop(paste("unable to calculate performance stats for validation data:",e))
  )
}
