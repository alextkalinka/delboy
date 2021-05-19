# Helper functions.
# Calculates the argmax KS distance for two distributions.
.argmax_ks_dist <- function(d1, d2){
  # Want to find argmax for d1 - d2 (sign matters) - i.e. d1 is distribution for FPs.
  e_d1 <- stats::ecdf(d1)
  e_d2 <- stats::ecdf(d2)
  range <- c(min(c(d1,d2), na.rm = T), max(c(d1,d2), na.rm = T))
  est <- NULL
  for(i in seq(range[1], range[2], by = 0.05)){
    est <- rbind(est, data.frame(val = i,
                                 # Direction of difference matters.
                                 ecdf_dist = e_d1(i) - e_d2(i)))
  }
  return (est %>%
            dplyr::filter(ecdf_dist == max(ecdf_dist, na.rm = T)) %>%
            # Latest point of maximal difference to remove as many FPs as possible.
            dplyr::filter(val == max(val, na.rm = T)))
}


#' calc_perf_pval_windows
#' 
#' Calculates sensitivity and FDR (after excluding predicted false positives) across a progressively larger window of DEG p-values.
#' 
#' @param data A data frame of DEG p-values.
#' @param pval_column A character string naming the p-value column.
#' @param gene_column A character string naming the gene column.
#' @param abs_lfc_column A character string naming the abs log fold change column.
#' @param deg_genes A character vector naming DEGs.
#' @return A data frame of performance estimates across the p-val thresholds with log fold change values maximising the KS distance at each p-value (`abs_Log2FoldChange.argmax_KS_dist`).
#' @export
#' @importFrom dplyr %>% arrange mutate filter select
#' @importFrom magrittr %<>%
#' @importFrom rlang !! sym
#' @importFrom stats ecdf
calc_perf_pval_windows <- function(data, pval_column, gene_column, abs_lfc_column, deg_genes){
  tryCatch({
    data %<>%
      dplyr::filter(!is.na(!! rlang::sym(pval_column))) %>%
      dplyr::arrange(!! rlang::sym(pval_column))
    num_tp <- length(deg_genes)
    sens <- fdr <- sens.excl <- fdr.excl <- pv <- lfc_mx <- NULL
    for(i in 1:nrow(data)){
      # TPs and FPs.
      tp <- unlist(data[1:i,gene_column]) %in% deg_genes
      tp.genes <- unlist(data[1:i,gene_column])[tp]
      fp.genes <- unlist(data[1:i,gene_column])[!tp]
      if(sum(tp) == 0 || sum(!tp) == 0) next
      sens <- append(sens, 100*sum(tp)/num_tp)
      fdr <- append(fdr, 100*sum(!tp)/i)
      # abs(LFC) that maximises separation between TPs and FPs.
      d1 <- unlist(data[1:i,] %>%
                     dplyr::filter(!! rlang::sym(gene_column) %in% fp.genes) %>%
                     dplyr::select(!! rlang::sym(abs_lfc_column)))
      d2 <- unlist(data[1:i,] %>%
               dplyr::filter(!! rlang::sym(gene_column) %in% tp.genes) %>%
               dplyr::select(!! rlang::sym(abs_lfc_column)))
      lfc_mxdiff <- .argmax_ks_dist(d1, d2)$val
      # Exclude predicted FPs below abs(LFC) threshold.
      td <- data[1:i,] %>%
        dplyr::filter(!! rlang::sym(abs_lfc_column) > lfc_mxdiff)
      tp <- unlist(td[,gene_column]) %in% deg_genes
      sens.excl <- append(sens.excl, 100*sum(tp)/num_tp)
      fdr.excl <- append(fdr.excl, 100*sum(!tp)/nrow(td))
      pv <- append(pv,unlist(data[i,pval_column]))
      lfc_mx <- append(lfc_mx, lfc_mxdiff)
    }
    ret <- data.frame(pvalue = rep(pv,2), Sensitivity.percent = c(sens,sens.excl),
                      FDR.percent = c(fdr,fdr.excl), 
                      abs_Log2FoldChange.argmax_KS_dist = rep(lfc_mx,2),
                      type = c(rep("All",length(pv)),rep("Excl_Pred_FP",length(pv))))
    return(ret)
  },
  error = function(e) stop(paste("unable to calculate performance across p-value windows:",e))
  )
}
