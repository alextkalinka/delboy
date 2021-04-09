#' assemble_elnet_hits
#'
#' Assembles Elastic-net hits together with log2 Fold Changes (from `DESeq2`) and an indication of whether a hit is predicted to be a false positive using an SVM fit.
#'
#' @param data Output from `delboy::combine_validation_original_hits`.
#' @param deseq2_res Output from `delboy::run_deseq2`.
#' @param svm_fit Output from `delboy::svm_false_positive_classification`.
#' @param svm_method A character string naming the SVM method. Can be one of `e1071` or `svmpath`.
#'
#' @return A data frame.
#' @export
#' @importFrom dplyr %>% mutate filter select arrange desc
#' @importFrom magrittr %<>%
#' @importFrom stats predict
assemble_elnet_hits <- function(data, deseq2_res, svm_fit, svm_method){
  tryCatch({
    data_svm <- data %>%
      dplyr::filter(hit_type == "Positive") %>%
      dplyr::select(id, log10_baseExpr, abs_log2FoldChange)
    if(svm_method == "e1071"){
      pred.fp <- stats::predict(svm_fit, as.matrix(data_svm %>% select(-id)))
    }else if(svm_method == "svmpath"){
      pred.fp <- stats::predict(svm_fit, as.matrix(data_svm %>% select(-id)),
                                lambda = tail(svm_fit$lambda,1),
                                type = "class")
      pred.fp <- ifelse(pred.fp == -1,0,1)
    }
    data_svm %<>%
      dplyr::mutate(fp = pred.fp)

    elnet_hits <- deseq2_res %>%
      dplyr::filter(id %in% data_svm$id) %>%
      dplyr::select(id, baseMean, log2FoldChange) %>%
      dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
      dplyr::mutate(Predicted_False_Positive = data_svm$fp[match(id, data_svm$id)])
  },
  error = function(e) stop(paste("unable to assemble elnet hits:",e))
  )
  return(elnet_hits)
}
