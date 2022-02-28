#' run_delboy_crispr
#'
#' Calculates a threshold value for a given metric to achieve a specified FDR.
#'
#' @param data A data frame containing count data for two different groups and their replicates.
#' @param controls A character string naming the control sample columns.
#' @param treatments A character string naming the treatment sample columns.
#' @param filter_cutoff A numerical value indicating the cutoff below which (summed across all replicates) a gene will be removed from the data. For example, to keep only genes with more than 1 TPM on average across both groups, set the cutoff to 10 if there are 10 replicates in total.
#' @param grna_column A character string naming the column containing gene names.
#' @param gene_column A character string naming the column containing gene names.
#' @param target_fdr Numeric value (0-1) giving the target FDR. Defaults to 0.1.
#' @param normalize_method A character string naming the read depth normalization method: `relative` (default) or `median_ratio`.
#' @param max.iter An integer value indicating the maximum number of validation sample combinations (default = 3). `NULL` indicates all sample combinations should be taken. `NA` indicates that the validation step should be skipped.
#'
#' @return A data frame.
#' @importFrom dplyr %>% arrange desc filter summarise n mutate
#' @importFrom rlang sym !!
#' @importFrom utils read.delim
#' @export
run_delboy_crispr <- function(data, controls, treatments, filter_cutoff, grna_column, gene_column, 
                              target_fdr = 0.1, normalize_method = "relative", max.iter = 3){
  tryCatch({
    # Random samples taken.
    set.seed(1)
    
    ### 1. Read data.
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
    .check_data_inputs(data, controls, treatments, gene_column, grna_column)
    
    ### 2B. Remove any irrelevant columns and re-order sample columns by groups.
    data.n <- data[,c(grna_column, gene_column, controls, treatments)]
    
    ### 3. Normalize read depth.
    if(normalize_method == "relative"){
      data <- delboy::normalize_library_depth_relative(data)
    }else if(normalize_method == "median_ratio"){
      data <- delboy::normalize_library_depth_median_ratio(data)
    }else{
      stop(paste("'normalize_method' must be one of 'relative' or 'median_ratio', got:",normalize_method))
    }
    
    ### 4. Run DESeq2 and harmonic mean p-value combination.
    res <- delboy::get_crispr_gene_level_hits(data.n, grna_column, gene_column, controls, treatments, target_fdr)
    
    ### 5. Estimate performance and learn FDR threshold.
    if(!is.na(max.iter)){
      
      ## Estimate non-null logFC distribution and adjust to ensure alignment with empirical distribution (from DESeq2).
      el <- delboy::adjust_nonnull_lfc_estimates(res, filter_prop = 0.05)
      
    }else{
      
    }
    
    
  },
  error = function(e) stop(paste("unable to run delboy crispr pipeline:",e))
  )
}
