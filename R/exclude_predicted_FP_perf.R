#' exclude_predicted_FP_perf
#'
#' Performance stats with predicted false positives removed (can impact both sensitivity and FDR/precision).
#'
#' @param data_svm The `data_svm` element from `delboy::svm_false_positive_classification`.
#' @param perf_stats The `performance_stats` from `delboy::evaluate_performance_rnaseq_calls`.
#' @param num_pos An integer giving the number of true positives in total.
#'
#' @return A data frame.
#' @export
#' @importFrom dplyr %>% mutate
#' @importFrom magrittr %<>%
exclude_predicted_FP_perf <- function(data_svm, perf_stats, num_pos){
  tryCatch({
    updates <- table(data_svm$Predicted_FP, data_svm$False_Positive)[2,]
    perf_stats %<>%
      dplyr::mutate(Num_true_calls = ifelse(Algorithm=="delboy",
                                            Num_true_calls - updates[1],
                                            Num_true_calls),
                    Sensitivity.percent = 100*Num_true_calls/num_pos,
                    Num_false_calls = ifelse(Algorithm=="delboy",
                                             Num_false_calls - updates[2],
                                             Num_false_calls),
                    Precision.percent = 100*Num_true_calls/(Num_true_calls+Num_false_calls),
                    FDR.percent = 100-Precision.percent)
  },
  error = function(e) stop(paste("unable to exclude predicted False Positives:",e))
  )
  return(perf_stats)
}
