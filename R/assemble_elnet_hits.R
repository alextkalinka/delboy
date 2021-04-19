#' assemble_elnet_hits
#'
#' Assembles Elastic-net hits together with log2 Fold Changes (from `DESeq2`) and an indication of whether a hit is predicted to be a false positive using an SVM fit.
#'
#' @param data Output from `delboy::combine_validation_original_hits`.
#' @param deseq2_res Output from `delboy::run_deseq2`.
#' @param decision_boundary A data frame providing co-ordinates for the SVM-learned decision boundary for FPs.
#'
#' @return A data frame.
#' @export
#' @importFrom dplyr %>% mutate filter select arrange desc
#' @importFrom magrittr %<>%
#' @importFrom stats predict
assemble_elnet_hits <- function(data, deseq2_res, decision_boundary){
  tryCatch({
    data_svm <- data %>%
      dplyr::filter(hit_type == "Positive") %>%
      dplyr::select(id, log10_baseExpr, abs_log2FoldChange)
    
    pred.fp <- delboy::predict_FP_delboy(data_svm, decision_boundary)$pred_FP
    
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
