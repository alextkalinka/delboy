#' svm_false_positive_classification
#'
#' Performs an SVM to learn the decision boundary between false and true positives for validation data.
#'
#' @param data A data frame as produced by `delboy::make_delboy_hit_comparison_table`.
#' @param kernel A character string naming the kernel type: `linear`, `polynomial`, `radial`.
#'
#' @return A list containing the SVM fit, a data frame containing predicted false positive status, and a data frame produced by `expand.grid` for use in `ggplot2::geom_contour`.
#' @export
#' @importFrom e1071 svm
#' @importFrom dplyr %>% filter select mutate
#' @importFrom magrittr %<>%
#' @importFrom stats predict
svm_false_positive_classification <- function(data, kernel){
  tryCatch({
    data_svm <- data %>%
      dplyr::filter(hit_type != "False_Negative") %>%
      dplyr::mutate(False_Positive = as.factor(ifelse(hit_type=="True_Positive",0,1))) %>%
      dplyr::select(False_Positive, abs_log2FoldChange, log10_baseExpr)
    svm_val <- e1071::svm(False_Positive ~ ., data_svm, kernel = kernel)
    pred_fp <- stats::predict(svm_val, data_svm)
    data_svm %<>%
      dplyr::mutate(Predicted_FP = pred_fp)
    data.grid <- expand.grid(log10_baseExpr = seq(from=0,to=5,length=400),
                             abs_log2FoldChange = seq(from=0,to=4,length=400)) %>%
      dplyr::mutate(Predicted_FP = as.numeric(stats::predict(svm_val, .)))
  },
  error = function(e) stop(paste("unable to perform SVM on validation data:",e))
  )
  return(list(svm_validation_fit = svm_val, data_svm = data_svm, grid = data.grid))
}
