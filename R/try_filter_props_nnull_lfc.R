#' try_filter_props_nnull_lfc
#'
#' Applies adjustments to estimated non-null logFC estimates to ensure alignment with empirical observations.
#'
#' @param res An object produced by `delboy::get_crispr_gene_level_hits`.
#' @param filter_prop A real number (default = 0.05) specifying what proportion of the lowest-abundance genes to filter out prior to making non-null estimation.
#'
#' @return An object of class delboy_logfc
#' @export
try_filter_props_nnull_lfc <- function(res, filter_prop){
  tryCatch({
    el <- delboy::adjust_nonnull_lfc_estimates(res, filter_prop = filter_prop)
    if(!el$fit_ok){
      el <- delboy::adjust_nonnull_lfc_estimates(res, filter_prop = 0.025)
      if(!el$fit_ok){
        el <- delboy::adjust_nonnull_lfc_estimates(res, filter_prop = 0.01)
        if(!el$fit_ok){
          el <- delboy::adjust_nonnull_lfc_estimates(res, filter_prop = 0)
        }
      }
    }
  },
  error = function(e) stop(paste("unable to try different filter props for non-null logFC distr estimate:",e))
  )
  return(el)
}
