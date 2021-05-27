# Helper functions.
.db_message <- function(msg, color){
  msg <- paste(msg,"\n",sep="")
  switch(color,
         blue = cat(crayon::blue(msg)),
         red = cat(crayon::red(msg)),
         green = cat(crayon::green(msg)),
         magenta = cat(crayon::magenta(msg))
         )
}


#' run_delboy
#'
#' Performs a differential-representation analysis using an elastic-net logistic regression approach for normalized count data that is split into two groups.
#'
#' @param data A data frame containing normalized count data for two different groups and their replicates. Can be a path to a file. Ideally the counts will already have been normalized to Transcripts Per Million (TPM), using, for example, the bias-aware quantification methods employed by `salmon` (Patro et al. 2017).
#' @param group_1 A character string naming the columns that belong to group 1.
#' @param group_2 A character string naming the columns that belong to group 2.
#' @param filter_cutoff A numerical value indicating the cutoff below which (summed across all replicates) a gene will be removed from the data. For example, to keep only genes with more than 1 TPM on average across both groups, set the cutoff to 10 if there are 10 replicates in total.
#' @param gene_column A character string naming the column containing gene names.
#' @param batches A named character vector identifying the batch structure with names identifying sample columns in the data input. The length must equal `length(group_1) + length(group_2)`. If `NULL`, there are no batches, or batches have already been corrected. Defaults to `NULL`. Batch correction will be conducted using `sva::ComBat` using non-parametric priors.
#' @param max.iter An integer value indicating the maximum number of validation samples (default = 10). `NULL` indicates all sample combinations should be taken.
#' @param bcorr_data_validation `NULL` if no batch (signal) corrected data is already available for validation. Otherwise, a data frame of treatment-corrected data should be supplied (to speed up validation, if already available). Defaults to `NULL`. Batch correction will be conducted using `sva::ComBat` using non-parametric priors.
#' @param alpha The elastic-net regression penalty, between 0 and 1 (default = 0.5). If `NULL`, alpha is chosen automatically.
#'
#' @return An object of class `delboy`. Access this object using `hits`, `plot.delboy`, `get_performance_stats`, and `get_deseq2_results`.
#' @seealso \code{\link{hits}}, \code{\link{plot.delboy}}, \code{\link{get_performance_stats}}, \code{\link{get_deseq2_results}}
#' @export
#' @importFrom dplyr left_join filter select arrange
#' @importFrom utils read.delim
#' @importFrom crayon blue red green magenta
#' @md
#' @author Alex T. Kalinka
#' @references
#' * Kalinka, A. T. 2020. Improving the sensitivity of differential-expression analyses for under-powered RNA-seq experiments. bioRxiv [10.1101/2020.10.15.340737](https://doi.org/10.1101/2020.10.15.340737).
#' * Patro, R. et al. 2017. Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods 14: 417-419.
run_delboy <- function(data, group_1, group_2, filter_cutoff, gene_column,
                       batches = NULL, max.iter = 10,
                       bcorr_data_validation = NULL, alpha = 0.5){
  # Random samples taken.
  set.seed(1)

  ### 1. Read data.
  if(!is.null(max.iter)){
    if(max.iter < 3) stop("'max.iter' < 3; at least 3 validation samples are required")
  }
  
  if(is.character(data)){
    tryCatch(
      data <- utils::read.delim(data, stringsAsFactors = F),
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
  if(!is.null(batches)){
    if(length(batches) != length(group_1) + length(group_2))
      stop(paste("the number of batched samples must equal the number of total samples across the two groups:\n",
                 "expected",length(group_1) + length(group_2),"but got",length(batches)))
    if(is.null(names(batches))) 
      stop("batches must be a named character vector with names referring to sample columns in the input data")
    if(length(setdiff(names(batches),colnames(data))))
      stop(paste("batch samples missing in data input:",setdiff(names(batches),colnames(data))))
  }

  ### 2B. Remove any irrelevant columns and re-order sample columns by groups.
  data <- data[,c(gene_column, group_1, group_2)]

  ### 3. Filter low count data prior to batch correction.
  cat("Filtering low-count cases...\n")
  tryCatch(data <- data[rowSums(data[,c(group_1, group_2)]) > filter_cutoff,],
           error = function(e) stop(paste("unable to filter low-count data:",e))
  )
  if(nrow(data) < 200)
    stop(paste("only",nrow(data),"rows remain after filtering: consider using a less stringent 'filter_cutoff' value. Cutoff used:",filter_cutoff))

  ### 4. Batch correction.
  if(!is.null(batches)){
    cat("Batch correcting data...\n")
    # Re-order batches to match re-ordered sample columns in data input.
    batches <- batches[match(colnames(data)[2:ncol(data)], names(batches))]
    data <- delboy::batch_correct(data, batches, gene_column)
  }else{
    batches <- NA
  }

  ### 5. Estimate parameters for performance evaluation.
  ## 5A. Run DESeq2 on the original dataset.
  deseq2_res <- delboy::run_deseq2(data, group_1, group_2, gene_column) %>%
    dplyr::left_join(data, by = c(id = gene_column))

  ## 5B. Number of non-null cases.
  cat("Estimating the number of non-null cases...")
  non.null <- suppressWarnings(delboy::estimate_number_non_nulls(deseq2_res$pvalue))
  cat(non.null$num.non_null,"\n")
  # Sanity-check non-null estimate.
  if(non.null$num.non_null == 0) stop("there are zero non-null cases estimated for this dataset")
  
  if(non.null$num.non_null < 40){
    .db_message("low estimate of number of non-null cases, setting to 40", "blue")
    non.null$num.non_null <- 40
  }

  ## 5C. Estimate non-null logFC distribution.
  cat("Estimating non-null logFC distribution...\n")
  lfdr.lfc <- suppressWarnings(delboy::estimate_nonnull_logfc_distr(deseq2_res$log2FoldChange))

  ### 6. Estimate delboy performance relative to DESeq2.
  ## 6A. Batch-correct real signal to create true-negative dataset.
  if(is.null(bcorr_data_validation)){
    cat("Batch correction to create signal-corrected data for validation...\n")
    batches.v <- c(rep("group_1",length(group_1)),rep("group_2",length(group_2)))
    data.bc <- delboy::batch_correct(data, batches.v, gene_column)
  }else{
    batches.v <- NA
    data.bc <- bcorr_data_validation
  }

  ## 6B. Performance evaluation.
  cat("Performance evaluation to validate results...\n")
  perf_eval <- suppressWarnings(
    delboy::evaluate_performance_deg_calls(data.bc, group_1, group_2, gene_column,
                                              max.iter,
                                              non.null$num.non_null,
                                              lfdr.lfc$non_null.lfc,
                                              lfdr.lfc$non_null.dens,
                                              alpha)
    )

  ### 7. Prep data for Elastic-net analysis.
  cat("Elastic-net logistic regression for hit detection...\n")
  data.elnet <- delboy::prep_elnet_data(data, group_1, group_2, gene_column)

  ### 8. Elastic-net logistic regression to identify differential representation.
  elnet.lr <- suppressWarnings(
    delboy::run_elnet_logistic_reg(as.matrix(data.elnet[,3:ncol(data.elnet)]),
                                   factor(data.elnet$treat),
                                   alpha = alpha)
  )

  ### 9. Combine hits with validation hit table to aid analysis of false positives.
  cat("Finishing up...\n")
  hits_orig_val <- delboy::combine_validation_original_hits(elnet.lr, deseq2_res,
                                                            perf_eval$delboy_hit_table)

  # Do we have any hits?
  if(sum(hits_orig_val$hit_type == "Positive") == 0)
    stop("* No hits found in elastic-net regression *")

  ### 10. Create final Elnet hit table.
  elnet_hits <- delboy::assemble_elnet_hits(hits_orig_val, deseq2_res,
                                            perf_eval$svm_validation$decision_boundary)

  ### 11. Update performance stats after excluding predicted False Positives.
  pstats_excl_pred_fp <- delboy::exclude_predicted_FP_perf(perf_eval$svm_validation$data_svm,
                                                           perf_eval$performance_stats,
                                                           # Total num TPs: num non-null * num val combinations.
                                                           non.null$num.non_null * perf_eval$num_val_combinations)

  ### 12. Build return object of class 'delboy'.
  ret <- list(data_input = data,
              data_val_bcorr = data.bc,
              sample_info = data.frame(SampleName = colnames(data)[2:ncol(data)],
                                       Group = c(rep(1,length(group_1)),rep(2,length(group_2))),
                                       Batch = batches,
                                       Batch.validation = batches.v),
              non_null = list(nonnull_number = non.null,
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
