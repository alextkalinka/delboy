# Helper function.
.collate_stats <- function(lfc, up, down, algorithm){
  true_up <- lfc[names(lfc) %in% up]
  true_up <- true_up[true_up > 0]
  true_down <- lfc[names(lfc) %in% down]
  true_down <- true_down[true_down < 0]
  sens <- 100*(length(true_up) + length(true_down))/length(lfc)
  prec <- 100*(length(true_up) + length(true_down))/(length(up) + length(down))
  fdr <- 100-prec
  return(data.frame(Algorithm = algorithm, Sensitivity.percent = sens, Precision.percent = prec,
                    FDR.percent = fdr))
}


#' perf_stats_rnaseq
#'
#' Extracts performance stats for differential expression calls made by `delboy` and `DESeq2`.
#'
#' @param delboy_res An object of class `delboy_elnet` as produced by `delboy::run_elnet_logistic_reg`.
#' @param deseq2_res The output from `delboy::run_deseq2`.
#' @param lfc_samp A named vector of logFC values where the names are gene names. These are the true positives with their associated logFC values. All other genes are true negatives.
#' @param padj_cutoff A numeric value between 0 and 1 indicating the adjusted p-value cutoff to use for `DESeq2` hits. Defaults to 0.1.
#'
#' @return A data frame of sensitivity and precision estimates.
#' @export
#' @importFrom magrittr %<>%
#' @importFrom dplyr filter %>%
perf_stats_rnaseq <- function(delboy_res, deseq2_res, lfc_samp, padj_cutoff = 0.1){
  tryCatch({
    # delboy stats.
    delboy_stats <- .collate_stats(lfc_samp, delboy_res$genes.up, delboy_res$genes.down, "delboy")

    # deseq2 stats.
    deseq2_res %<>%
      dplyr::filter(padj < padj_cutoff)
    deseq2_up <- deseq2_res %>%
      filter(log2FoldChange > 0)
    deseq2_down <- deseq2_res %>%
      filter(log2FoldChange < 0)

    deseq2_stats <- .collate_stats(lfc_samp, deseq2_up$id, deseq2_down$id, "DESeq2")

    pstats <- rbind(delboy_stats, deseq2_stats)
  },
  error = function(e) stop(paste("unable to extract performance stats for RNAseq data:",e))
  )
  return(pstats)
}
