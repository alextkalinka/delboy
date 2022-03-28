#' combine_harmonic_mean_pvals
#'
#' Combines p-values for groups of dependent tests using a harmonic mean approach.
#'
#' @param data A data frame containing p-values and a grouping variable.
#' @param pval_col A character string naming the p-value column.
#' @param group_col A character string naming the grouping column.
#' @param target_fdr A real number (0-1) specifying the FDR.
#'
#' @return An object of class delboy_logfc
#' @export
#' @importFrom harmonicmeanp qharmonicmeanp
#' @importFrom dplyr %>% group_by summarise ungroup mutate
#' @importFrom rlang sym !!
#' @references Wilson, D. J. 2019.The harmonic mean p-value for combining dependent tests. PNAS 116 (4): 1195-1200.
combine_harmonic_mean_pvals <- function(data, pval_col, group_col, target_fdr){
  tryCatch({
    L <- nrow(data)
    alpha.L <- harmonicmeanp::qharmonicmeanp(target_fdr, L)
    pval_sym <- rlang::sym(pval_col)
    group_sym <- rlang::sym(group_col)
    # Combine p-vals using HMP.
    ret <- data %>%
      dplyr::mutate(weight = 1/L) %>%
      dplyr::group_by(!! group_sym) %>%
      dplyr::summarise(# Summary stats for logFC and abundance.
                       mean_log2FoldChange = mean(log2FoldChange, na.rm=T),
                       median_log2FoldChange = median(log2FoldChange, na.rm=T),
                       mean_baseMean = mean(baseMean, na.rm=T),
                       sd_log2FoldChange = sd(log2FoldChange, na.rm=T),
                       sd_baseMean = sd(baseMean, na.rm=T),
                       # HMP calculation.
                       sum_w = sum(weight),
                       pval_mult_test_thresh = sum_w * alpha.L,
                       pvalue.harmonic_mean = sum_w/sum(weight/!!pval_sym)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(significant_hit = pvalue.harmonic_mean < pval_mult_test_thresh)
  },
  error = function(e) stop(paste("unable to combine harmonic mean p-values:",e))
  )
  return(ret)
}
