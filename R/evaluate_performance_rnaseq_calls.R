#' evaluate_performance_rnaseq_calls
#'
#' Evaluates the performance of `delboy` on RNAseq data by comparison with `DESeq2` output on the original input data for an experiment, controlling for real signal, and adding known signal (using `seqgendiff`'s binomial-thinning approach) for a sampled number of genes from a logFC distribution matching as closely as possible the real one.
#'
#' @param data A data frame of normalized (and batch-corrected, if necessary) counts for a set of samples.
#' @param group_1 A character vector naming the columns that belong to group 1.
#' @param group_2 A character vector naming the columns that belong to group 2.
#' @param gene_column A character string naming the column containing gene names.
#' @param num_non_null An integer value indicating the number of genes to add signal to.
#' @param lfc A vector of logFC values for non-null cases.
#' @param lfc_dist A vector of density estimates for the logFC values given in `lfc`.
#'
#' @return An object of class `delboy_performance`
#' @export
evaluate_performance_rnaseq_calls <- function(data, group_1, group_2, gene_column,
                                              num_non_null, lfc, lfc_dist){
  # 1. Batch-correct real signal to create true negative dataset.
  data.bc <- delboy::batch_correct(data, group_1, group_2, gene_column)

  # 2.
}
