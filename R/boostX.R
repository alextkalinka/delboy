#' boostX
#' 
#' Takes the output from a given differential expression/representation algorithm and tries to boost sensitivity while controlling the FDR.
#' 
#' @param data A data frame containing normalized count data for two different groups and their replicates. Can be a path to a file. Ideally the counts will already have been normalized to Transcripts Per Million (TPM), using, for example, the bias-aware quantification methods employed by `salmon` (Patro et al. 2017).
#' @param diff_exp Differential expression/representation results from a particular algorithm.
#' @param group_1 A character string naming the columns that belong to group 1.
#' @param group_2 A character string naming the columns that belong to group 2.
#' @param gene_column A character string naming the column containing gene names.
#' @param pvalue_column A character string naming the pvalue column in `diff_exp`.
#' @param max.iter An integer value indicating the maximum number of validation samples (default = 10). `NULL` indicates all sample combinations should be taken.
#' @param target_fdr A numerical value (0-1) indicating the target FDR. Defaults to 0.1.
#' @param bcorr_data_validation `NULL` if no batch (signal) corrected data is already available for validation. Otherwise, a data frame of treatment-corrected data should be supplied (to speed up validation, if already available). Defaults to `NULL`. Batch correction will be conducted using `sva::ComBat` using non-parametric priors.
#' @return A list containing the following elements:
#' @md
#' @author Alex T. Kalinka
#' @export
boostX <- function(data, diff_exp, group_1, group_2, gene_column, pvalue_column, max.iter = 10, target_fdr = 0.1,
                   bcorr_data_validation = NULL){
  # Random samples taken.
  set.seed(1)
  
  ### 1. Remove any irrelevant columns and re-order sample columns by groups.
  data <- data[,c(gene_column, group_1, group_2)]
  
  ### 2. Number of non-null cases.
  cat("Estimating the number of non-null cases...")
  non.null <- suppressWarnings(delboy::estimate_number_non_nulls(unlist(diff_exp[,pvalue_column])))
  cat(non.null$num.non_null,"\n")
  # Sanity-check non-null estimate.
  if(non.null$num.non_null == 0) stop("there are zero non-null cases estimated for this dataset")
  
  if(non.null$num.non_null < 40){
    .db_message("low estimate of number of non-null cases, setting to 40", "blue")
    non.null$num.non_null <- 40
  }
  
  ### 3. Estimate non-null logFC distribution.
  cat("Estimating non-null logFC distribution...\n")
  lfdr.lfc <- suppressWarnings(delboy::estimate_nonnull_logfc_distr(deseq2_res$log2FoldChange))
  
  ## 4. Batch-correct real signal to create true-negative dataset.
  if(is.null(bcorr_data_validation)){
    cat("Batch correction to create signal-corrected data for validation...\n")
    batches.v <- c(rep("group_1",length(group_1)),rep("group_2",length(group_2)))
    data.bc <- delboy::batch_correct(data, batches.v, gene_column)
  }else{
    batches.v <- NA
    data.bc <- bcorr_data_validation
  }
  
  ### 5. Performance evaluation.
  cat("Performance evaluation to validate results...\n")
  
  
}
