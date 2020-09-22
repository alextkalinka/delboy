#' get_deseq2_results
#'
#' Retrieve a data frame of `DESeq2` results from an object of class `delboy`.
#'
#' @param delboy_res An object of class `delboy`.
#'
#' @return A data frame.
#' @export
get_deseq2_results <- function(delboy_res){
  if(!inherits(delboy_res,"delboy")) stop(paste("expecting an object of class 'delboy', got:",class(delboy_res)))

  return(delboy_res$deseq2_results)
}
