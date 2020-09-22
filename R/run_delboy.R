#' run_delboy
#'
#' Performs a differential-representation analysis using an elastic-net logistic regression approach for (normalized) count data that is split into two groups.
#'
#' @param data A data frame containing normalized count data for two different groups and their replicates. Can be a path to a file. Ideally the counts will already have been normalized to Transcripts Per Million (TPM), using, for example, the bias-aware quantification methods employed by `salmon` (Patro et al. 2017).
#' @param group_1 A character string naming the columns that belong to group 1.
#' @param group_2 A character string naming the columns that belong to group 2.
#' @param filter_cutoff A numerical value indicating the cutoff below which (summed across all replicates) a gene (or gRNA) will be removed from the data. For example, to keep only genes with more than 1 TPM on average across both groups, set the cutoff to 10 if there are 10 replicates in total.
#' @param gene_column A character string naming the column containing gene names.
#' @param batches A character vector identifying the batch structure of the experiment. The length must equal `length(group_1) + length(group_2)`. If `NULL`, there are no batches, or batches have already been corrected. Defaults to `NULL`.
#' @param batch_corr_method A character string naming the batch correction method. Can be one of `combat_np` (non-parametric) or `combat_p` (parametric). Defaults to `combat_np`. Ignored if `batches` is `NULL`.
#' @param bcorr_data_validation `NULL` if no batch corrected data is already available. Otherwise, a data frame of treatment-corrected data should be supplied (to speed up validation, if already available).
#'
#' @return An object of class `delboy`. Access this object using `delboy::hits`, and `delboy::plot.delboy`.
#' @seealso \code{\link{hits}}, \code{\link{plot.delboy}}
#' @export
#' @importFrom dplyr left_join filter select arrange
#' @md
#' @references
#' * Kalinka, A. T. 2020. Improving the sensitivity of differential-expression analyses for under-powered RNA-seq experiments. bioRxiv.
#' * Patro, R. et al. 2017. Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods 14: 417-419.
run_delboy <- function(data, group_1, group_2, filter_cutoff, gene_column,
                       batches = NULL, batch_corr_method = "combat_np",
                       bcorr_data_validation = NULL){
  # Random samples taken.
  set.seed(1)

  ### 1. Read data.
  if(is.character(data)){
    tryCatch(
      data <- read.delim(data, stringsAsFactors = F),
      error = function(e) stop(paste("unable to read data file",data))
    )
  }else if(!is.data.frame(data)){
    stop(paste("expected 'data' to be a data frame, got instead an object of class",class(data)))
  }

  ### 2A. Sanity checks.
  cat("Checking data...\n")
  if(any(!group_1 %in% colnames(data)))
    stop(paste("unable to find the following group 1 columns:",setdiff(group_1,colnames(data))))
  if(any(!group_2 %in% colnames(data)))
    stop(paste("unable to find the following group 2 columns:",setdiff(group_2,colnames(data))))
  if(!gene_column %in% colnames(data))
    stop(paste("unable to find the gene column",gene_column))
  if(!is.null(bcorr_data_validation)){
    if(!is.data.frame(bcorr_data_validation))
      stop(paste("'bcorr_data_validation' should be a data frame, instead got",
                 class(bcorr_data_validation)))
  }

  ### 2B. Remove any irrelevant columns.
  data <- data[,c(gene_column, group_1, group_2)]

  ### 3. Normalize counts.
  if(!is.null(normalize)){

  }

  ### 4. Filter low count data prior to batch correction.
  cat("Filtering low-count cases...\n")
  tryCatch(data <- data[rowSums(data[,c(group_1, group_2)]) > filter_cutoff,],
           error = function(e) stop(paste("unable to filter low-count data:",e))
  )
  if(nrow(data) < 200)
    stop(paste("only",nrow(data),"rows remain after filtering: consider using a less stringent 'filter_cutoff' value. Cutoff used:",filter_cutoff))

  ### 5. Batch correction.
  if(!is.null(batches)){
    cat("Batch correcting data...\n")
    data <- delboy::batch_correct(data, group_1, group_2, gene_column)
  }

  ### 6. Estimate parameters for performance evaluation.
  ## 6A. Run DESeq2 on the original dataset.
  deseq2_res <- delboy::run_deseq2(data, group_1, group_2, gene_column) %>%
    dplyr::left_join(data, by = c(id = gene_column))

  ## 6B. Number of non-null cases.
  cat("Estimating the number of non-null cases...")
  non.null <- suppressWarnings(delboy::estimate_number_non_nulls(deseq2_res$pvalue))
  cat(non.null$num.non_null,"\n")
  # Sanity-check non-null estimate.
  if(non.null$num.non_null < 40){
    cat("Low estimate of number of non-null cases, setting to 40\n")
    non.null$num.non_null <- 40
  }

  ## 6C. Estimate non-null logFC distribution.
  cat("Estimating non-null logFC distribution...\n")
  lfdr.lfc <- suppressWarnings(delboy::estimate_nonnull_logfc_distr(deseq2_res$log2FoldChange))

  ### 8. Estimate delboy performance relative to DESeq2.
  ## 8A. Batch-correct real signal to create true-negative dataset.
  if(is.null(bcorr_data_validation)){
    cat("Batch correction to create signal-corrected data for validation...\n")
    data.bc <- delboy::batch_correct(data, group_1, group_2, gene_column)
  }else{
    data.bc <- bcorr_data_validation
  }

  ## 8B. Performance evaluation.
  cat("Performance evaluation to validate results...\n")
  perf_eval <- suppressWarnings(
    delboy::evaluate_performance_rnaseq_calls(data.bc, group_1, group_2, gene_column,
                                              non.null$num.non_null,
                                              lfdr.lfc$non_null.lfc,
                                              lfdr.lfc$non_null.dens)
    )

  ### 9. Prep data for Elastic-net analysis.
  cat("Elastic-net logistic regression for hit detection...\n")
  data.elnet <- delboy::prep_elnet_data(data, group_1, group_2, gene_column)

  ### 10. Elastic-net logistic regression to identify differentially-represented genes or gRNAs.
  elnet.lr <- suppressWarnings(
    delboy::run_elnet_logistic_reg(as.matrix(data.elnet[,3:ncol(data.elnet)]),
                                   factor(data.elnet$treat),
                                   alpha = 0.5)
  )

  ### 11. Combine hits with validation hit table to aid analysis of false positives.
  cat("Finishing up...\n")
  hits_orig_val <- delboy::combine_validation_original_hits(elnet.lr, deseq2_res,
                                                            perf_eval$delboy_hit_table)

  ### 12. Create final Elnet hit table.
  elnet_hits <- delboy::assemble_elnet_hits(hits_orig_val, deseq2_res,
                                            perf_eval$svm_validation$svm_validation_fit)

  ### 13. Update performance stats after excluding predicted False Positives.
  pstats_excl_pred_fp <- delboy::exclude_predicted_FP_perf(perf_eval$svm_validation$data_svm,
                                                           perf_eval$performance_stats,
                                                           non.null$num.non_null)

  ### 14. Build object of class 'delboy'.
  ret <- list(non_null = list(nonnull_number = non.null,
                              nonnull_lfc = lfdr.lfc),
              performance_eval = perf_eval,
              data_elnet = data.elnet,
              elnet_hits = elnet_hits,
              elnet_results = elnet.lr,
              deseq2_results = deseq2_res,
              hits_original_validation = hits_orig_val,
              performance_stats_corr_FP = pstats_excl_pred_fp)
  class(ret) <- "delboy"
  return(ret)
}
