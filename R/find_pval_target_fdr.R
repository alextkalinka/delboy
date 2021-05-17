#' find_pval_target_fdr
#' 
#' Returns the unadjusted p-value that corresponds to a target FDR using a LOESS fit to the FDR estimates.
#' 
#' @param data A data frame containing FDR estimates and the corresponding unadjusted p-values.
#' @param target_FDR A numerical value (0-1) giving the target FDR.
#' @return A list containing the following elements:
#' * `pvalue` A numerical value giving the p-value corresponding to the target FDR.
#' * `loess_fit_fdr_all` An object of class `loess`.
#' * `loess_fit_fdr_excl` An object of class `loess`.
#' * `loess_fit_sens_all` An object of class `loess`.
#' * `loess_fit_sens_excl` An object of class `loess`.
#' * `target_FDR` A numerical value (0-1) giving the target FDR.
#' @md
#' @export
#' @importFrom dplyr %>% filter
#' @importFrom stats loess predict
find_pval_target_fdr <- function(data, target_FDR){
  tryCatch({
    
    loess_fdr_all <- stats::loess()
  },
  error = function(e) stop(paste("unable to find p-value corresponding to a target FDR:",e))
  )
}
