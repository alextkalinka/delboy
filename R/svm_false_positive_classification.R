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
#' @importFrom svmpath svmpath
#' @importFrom dplyr %>% filter select mutate
#' @importFrom magrittr %<>%
#' @importFrom stats predict
svm_false_positive_classification <- function(data, kernel){
  tryCatch({
    data_svm <- data %>%
      dplyr::filter(hit_type != "False_Negative" & !is.na(hit_type) &
                      !is.na(abs_log2FoldChange) &
                      !is.na(log10_baseExpr) & !is.infinite(log10_baseExpr) &
                      !is.infinite(abs_log2FoldChange)) %>%
      dplyr::mutate(False_Positive = as.factor(ifelse(hit_type=="True_Positive",0,1))) %>%
      dplyr::select(False_Positive, abs_log2FoldChange, log10_baseExpr)

    # Try 'svmpath::svmpath' first - due to fully optimized regularization, there is minimal over-fitting.
    # If there is a singularity or other error, then try 'e1071::svm'.
    e1071 <- FALSE
    tryCatch({
      svm_val <- svmpath::svmpath(as.matrix(data_svm %>% select(-False_Positive)),
                                  ifelse(data_svm$False_Positive == 0,-1,1),
                                  param.kernel = 2,
                                  trace = F,
                                  ridge = 1e-12, lambda.min = 1e-3, eps = 1e-9)
      pred_fp <- as.numeric(stats::predict(svm_val, as.matrix(data_svm %>% dplyr::select(-False_Positive)),
                                           lambda = tail(svm_val$lambda,1),
                                           type = "class"))
      pred_fp <- ifelse(pred_fp == -1,0,1)
      data.grid <- expand.grid(abs_log2FoldChange = seq(from=0,to=4,length=400),
                               log10_baseExpr = seq(from=0,to=5,length=400)) %>%
        dplyr::mutate(Predicted_FP = as.numeric(stats::predict(svm_val, as.matrix(.),
                                                               lambda = tail(svm_val$lambda,1),
                                                               type = "class")),
                      Predicted_FP = ifelse(Predicted_FP == -1,0,1))
    },
    error = function(e) e1071 <<- TRUE
    )

    if(e1071){
      svm_val <- e1071::svm(False_Positive ~ ., data_svm, kernel = kernel)
      pred_fp <- stats::predict(svm_val, data_svm)
      data.grid <- expand.grid(abs_log2FoldChange = seq(from=0,to=4,length=400),
                               log10_baseExpr = seq(from=0,to=5,length=400)) %>%
        dplyr::mutate(Predicted_FP = as.numeric(stats::predict(svm_val, as.matrix(.))))
    }

    data_svm %<>%
      dplyr::mutate(Predicted_FP = pred_fp)

    svm_method <- ifelse(e1071,"e1071","svmpath")
  },
  error = function(e) stop(paste("unable to perform SVM on validation data:",e))
  )
  return(list(svm_validation_fit = svm_val, data_svm = data_svm, grid = data.grid, svm.method = svm_method))
}
