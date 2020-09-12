#' combine_validation_original_hits
#'
#' Combine validation hits with hits from unmodified input data to aid analysis of false positives.
#'
#' @param elnet_lr_res An object of class `delboy_elnet`, the output from running `delboy::run_elnet_logistic_reg`.
#' @param deseq2_res The output from `delboy::run_deseq2`.
#' @param delboy_validation_hits The output from running `delboy::make_delboy_hit_comparison_table`.
#'
#' @return A data frame.
#' @export
#' @importFrom dplyr filter select mutate %>%
combine_validation_original_hits <- function(elnet_lr_res, deseq2_res, delboy_validation_hits){
  tryCatch({
    elnet_hits <- c(elnet_lr_res$genes.up, elnet_lr_res$genes.down)
    ret <- deseq2_res %>%
      dplyr::filter(id %in% elnet_hits) %>%
      dplyr::mutate(log10_baseExpr = log10(baseMean),
                    abs_log2FoldChange = abs(log2FoldChange),
                    hit_type = "Positive",
                    data_type = "Original") %>%
      dplyr::select(id,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,
                    log10_baseExpr,abs_log2FoldChange,hit_type,data_type) %>%
      rbind(delboy_validation_hits %>%
              dplyr::mutate(data_type = "Validation")) %>%
      dplyr::mutate(data_type = factor(data_type, levels=c("Validation","Original")))
  },
  error = function(e) stop(paste("unable to combine validation and original hits:",e))
  )
  return(ret)
}
