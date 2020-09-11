#' make_delboy_hit_comparison_table
#'
#' Make a TP, FN, FP data frame to aid analysis of `delboy` hits.
#'
#' @param elnet_lr_res An object of class `delboy_elnet`, the output from running `delboy::run_elnet_logistic_reg`.
#' @param deseq2_res The output from `delboy::run_deseq2`.
#' @param lfc_samp A named vector of logFC values where the names are gene names. These are the true positives with their associated logFC values. All other genes are true negatives.
#'
#' @return A data frame.
#' @export
#' @importFrom dplyr filter %>% mutate case_when
make_delboy_hit_comparison_table <- function(elnet_lr_res, deseq2_res, lfc_samp){
  tryCatch({
    elnet_hits <- c(elnet_lr_res$genes.up, elnet_lr_res$genes.down)
    all_hits <- unique(c(names(lfc_samp), elnet_hits))
    TP <- names(lfc_samp)[names(lfc_samp) %in% elnet_hits]
    FN <- names(lfc_samp)[!names(lfc_samp) %in% elnet_hits]
    FP <- elnet_hits[!elnet_hits %in% names(lfc_samp)]
    ret <- deseq2_res %>%
      dplyr::filter(id %in% all_hits) %>%
      dplyr::mutate(hit_type = dplyr::case_when(id %in% TP ~ "True_Positive",
                                                id %in% FN ~ "False_Negative",
                                                id %in% FP ~ "False_Positive"),
                    hit_type = factor(hit_type,levels=c("True_Positive","False_Negative",
                                                        "False_Positive")))
  },
  error = function(e) stop(paste("unable to build delboy-DESeq2 comparison data frame:",e))
  )
  return(ret)
}
