#' hits
#'
#' Extract hits from an object of class `delboy`.
#'
#' @param delboy_res An object of class `delboy` as produced by `delboy::run_delboy`.
#' @param remove_predicted_false_positives Logical indicating whether predicted false positives should be removed. Defaults to `FALSE`.
#'
#' @return A data frame.
#' @export
#' @importFrom dplyr filter
#' @importFrom magrittr %<>%
hits <- function(delboy_res, remove_predicted_false_positives = FALSE){
  if(!inherits(delboy_res,"delboy"))
    stop(paste("expecting an object of class 'delboy', got",class(delboy_res)))

  hits <- delboy_res$elnet_hits

  if(remove_predicted_false_positives){
    hits %<>%
      dplyr::filter(Predicted_False_Positive == 0)
  }
  return(hits)
}
