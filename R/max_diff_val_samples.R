#' max_diff_val_samples
#'
#' Extracts a maximally distinct set of treatment samples for performance validation.
#'
#' @param all_comb A matrix, as produced by `delboy::all_combinations_treat_samples`.
#' @param max.iter An integer value indicating the maximum number of validation sample combinations.
#'
#' @return A matrix of validation 'treatment' sample combinations.
#' @export
#' @importFrom stats median
max_diff_val_samples <- function(all_comb, max.iter){
  tryCatch({
    tot_combs <- ncol(all_comb)
    if(!is.null(max.iter)){
      num_combs <- min(max.iter, tot_combs)
      # Try to use maximally different sets of samples: 
      # always take first and last as these guaranteed to have a max edit distance.
      val_inds <- c(1, seq(round(stats::median(c(1, tot_combs))), tot_combs-1, 
                           by = ceiling((tot_combs-2)/(num_combs-2))), tot_combs)
      all_comb <- all_comb[,val_inds]
    }
  },
  error = function(e) stop(paste("unable to extract maximally distinct validation samples:",e))
  )
  return(all_comb)
}
