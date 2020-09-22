#' get_performance_stats
#'
#' Extract performance stats from an object of class `delboy`.
#'
#' @param delboy_res An object of class `delboy`.
#'
#' @return A data frame.
#' @export
get_performance_stats <- function(delboy_res){
  if(!inherits(delboy_res,"delboy")) stop(paste("expecting an object of class 'delboy', got:",class(delboy_res)))

  return(delboy_res$performance_stats_corr_FP)
}
