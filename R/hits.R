#' hits
#'
#' Extract hits from an object of class `delboy` or `delboy_crispr`.
#'
#' @param delboy_res An object of class `delboy` or `delboy_crispr` as produced by `delboy::run_delboy`.
#' @param all_res Logical indicating whether to return all genes, or just significant hits. Default = `TRUE`.
#' @param remove_predicted_false_positives Logical indicating whether predicted false positives should be removed. Defaults to `TRUE`.
#' @param dir A character string naming the direction of hits required: `pos`, `neg`, if `NULL`, then return both.
#'
#' @return A data frame.
#' @export
#' @importFrom dplyr filter arrange
#' @importFrom magrittr %<>%
hits <- function(delboy_res, all_res = FALSE, remove_predicted_false_positives = TRUE, dir = NULL){
  if(!inherits(delboy_res,"delboy") && !inherits(delboy_res,"delboy_crispr"))
    stop(paste("expecting an object of class 'delboy' or 'delboy_crispr', got",class(delboy_res)))

  if(inherits(delboy_res,"delboy")){
    hits <- delboy_res$elnet_hits
    if(remove_predicted_false_positives){
      hits %<>%
        dplyr::filter(Predicted_False_Positive == 0)
    }
  }else{
    hits <- delboy_res$results_gene
    if(!all_res){
      if(is.null(dir)){
        hits %<>%
          dplyr::filter((significant_hit.pos & !Predicted_False_Positive.pos) | (significant_hit.neg & !Predicted_False_Positive.neg))
      }else if(dir == "pos"){
        hits %<>%
          dplyr::filter((significant_hit.pos & !Predicted_False_Positive.pos))
      }else if(dir == "neg"){
        hits %<>%
          dplyr::filter((significant_hit.neg & !Predicted_False_Positive.neg)) %>%
          dplyr::arrange(pvalue.harmonic_mean.neg)
      }else{
        stop(paste("'dir' should be 'pos', 'neg', or NULL, got:",dir))
      }
    }
  }
  return(hits)
}
