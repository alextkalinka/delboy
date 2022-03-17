#' adjust_nonnull_lfc_estimates
#'
#' Applies adjustments to estimated non-null logFC estimates to ensure alignment with empirical observations.
#'
#' @param res An object produced by `delboy::get_crispr_gene_level_hits`.
#' @param filter_prop A real number (default = 0.05) specifying what proportion of the lowest-abundance genes to filter out prior to making non-null estimation.
#'
#' @return An object of class delboy_logfc
#' @export
#' @importFrom dplyr %>% between
#' @importFrom magrittr %<>%
#' @importFrom stats quantile median
adjust_nonnull_lfc_estimates <- function(res, filter_prop = 0.05){
  tryCatch({
    if(!dplyr::between(filter_prop,0,1) && filter_prop < 1) 
      stop(paste("''filter_prop' must be in the range [0-1), got:",filter_prop))
    # 1. Filter out low-abundance genes to stabilise non-null estimate.
    lfc <- res$hmp_gene_pos
    ab_thr <- stats::quantile(lfc$mean_baseMean, probs = filter_prop)
    lfc %<>%
      dplyr::filter(mean_baseMean > ab_thr)
    
    # 2. Estimate non-null logFC distribution.
    el <- delboy::estimate_nonnull_logfc_distr(lfc$mean_log2FoldChange, signed = TRUE)
    if(!el$fit_ok) return(el)
    
    # 3. Adjust downwards if median is greater than empirical observation for significant hits.
    mdn_orig.pos <- stats::median((res$hmp_gene_pos %>%
                              dplyr::filter(significant_hit) %>%
                              dplyr::summarise(ml = median(mean_log2FoldChange)))$ml, na.rm=T)
    mdn_orig.neg <- stats::median((res$hmp_gene_neg %>%
                              dplyr::filter(significant_hit) %>%
                              dplyr::summarise(ml = median(mean_log2FoldChange)))$ml, na.rm=T)
    # Sample logFC values from estimated lfc distribution.
    lfc_pos <- sample(el$non_null.pos.lfc, 1e3, 
                      prob = el$non_null.pos.dens/sum(el$non_null.pos.dens), replace = T)
    lfc_neg <- sample(el$non_null.neg.lfc, 1e3, 
                      prob = el$non_null.neg.dens/sum(el$non_null.neg.dens), replace = T)
    # Do we need to adjust logFC to align with empirical obs?
    if(!is.na(mdn_orig.pos)){
      lpd <- stats::median(lfc_pos, na.rm=T) - mdn_orig.pos
      if(lpd > 0) el$non_null.pos.lfc <- el$non_null.pos.lfc - lpd
    }
    if(!is.na(mdn_orig.neg)){
      lnd <- stats::median(lfc_neg, na.rm=T) - mdn_orig.neg
      if(lnd < 0) el$non_null.neg.lfc <- el$non_null.neg.lfc - lnd
    }
  },
  error = function(e) stop(paste("unable to adjust non-null logFC estimates:",e))
  )
  return(el)
}
