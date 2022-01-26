#' all_combinations_treat_samples
#' 
#' Returns all combinations of treatment samples for validation purposes.
#' 
#' @param samples A character vector of sample names.
#' @param num_treat An integer specifying the number of treatment samples.
#' 
#' @return A matrix of combinations (each column corresponds to a separate combination). There will be `choose(length(samples), length(num_treat))` combinations (columns).
#' @export
#' @importFrom combinat combn
all_combinations_treat_samples <- function(samples, num_treat){
  tryCatch({
    if(num_treat >= length(samples))
      stop(paste("'length(num_treat)' must be less than length(samples):",length(num_treat),"vs",length(samples)))
    comb <- combinat::combn(samples, num_treat)
    return(comb)
  },
  error = function(e) stop(paste("unable to create all combinations of treatment samples for validation:",e))
  )
}
