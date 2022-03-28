# Helper functions.
.print_progress <- function(msg){
  msg <- paste(date(), "*** delboy:",msg,"\n")
  cat(msg)
}


#' run_delboy_crispr
#'
#' Calculates a threshold value for a given metric to achieve a specified FDR.
#'
#' @param data A data frame containing count data for two different groups and their replicates.
#' @param controls A character string naming the control sample columns.
#' @param treatments A character string naming the treatment sample columns.
#' @param grna_column A character string naming the column containing gene names.
#' @param gene_column A character string naming the column containing gene names.
#' @param target_fdr Numeric value (0-1) giving the target FDR. Defaults to 0.1.
#' @param normalize_method A character string naming the read depth normalization method: `relative` (default) or `median_ratio`.
#' @param estim_lfc_distr Logical indicating whether to estimate the non-null log-fold change distribution in the data (default = `FALSE`).
#'
#' @return An object of class `delboy_crispr`.
#' @importFrom utils read.delim packageVersion
#' @export
run_delboy_crispr <- function(data, controls, treatments, grna_column, gene_column, 
                              target_fdr = 0.1, normalize_method = "relative", estim_lfc_distr = FALSE){
  cat(paste("*** delboy version:",utils::packageVersion("delboy"),"\n"))
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
    .print_progress("Checking data...")
    .check_data_inputs(data, controls, treatments, gene_column, grna_column)
    
    ### 2B. Remove any irrelevant columns and re-order sample columns by groups.
    data <- data[,c(grna_column, gene_column, controls, treatments)]
    
    ### 3. Normalize read depth.
    .print_progress("Normalizing read depth...")
    if(normalize_method == "relative"){
      data <- delboy::normalize_library_depth_relative(data, NULL)
    }else if(normalize_method == "median_ratio"){
      data <- delboy::normalize_library_depth_median_ratio(data)
    }else{
      stop(paste("'normalize_method' must be one of 'relative' or 'median_ratio', got:",normalize_method))
    }
    
    .print_progress("Running DESeq2 and combining p-values...")
    ### 4. Run DESeq2 and harmonic mean p-value combination.
    res <- delboy::get_crispr_gene_level_hits(data, grna_column, gene_column, controls, treatments, target_fdr)
    
    ### 5. Flag up predicted False Positives.
    res$hmp_gene_pos <- delboy::mark_up_FPs(res$hmp_gene_pos, "pos")
    
    res$hmp_gene_neg <- delboy::mark_up_FPs(res$hmp_gene_neg, "neg")
    
    .print_progress("Finishing up...")

    #### 6. Build single results data frame for both pos and neg.
    res_both <- delboy::combine_pos_neg_res(res$hmp_gene_pos, res$hmp_gene_neg)
    
    if(estim_lfc_distr){
      ### 7. Estimate non-null logFC distribution and adjust to ensure alignment with empirical distribution (from DESeq2).
      el <- delboy::try_filter_props_nnull_lfc(res, filter_prop = 0.05)
      if(el$fit_ok){
        lfc_distr <- data.frame(lfc = c(el$non_null.pos.lfc, el$non_null.neg.lfc),
                                dens = c(el$non_null.pos.dens/sum(el$non_null.pos.dens),
                                         el$non_null.neg.dens/sum(el$non_null.neg.dens)),
                                type = c(rep("pos",length(el$non_null.pos.lfc)), rep("neg",length(el$non_null.neg.lfc))),
                                stringsAsFactors = F)
      }else{
        lfc_distr <- NA
      }
    }else{
      el <- lfc_distr <- NA
    }
    
    ### 8. Build return object.
    ret <- list(data.norm = data,
                results_gene = res_both,
                results_all = res,
                norm.method = normalize_method,
                target_fdr = target_fdr,
                logfc_nonnull_distr = lfc_distr)
    class(ret) <- "delboy_crispr"
    
  },
  error = function(e) stop(paste("unable to run delboy crispr pipeline:",e))
  )
  return(ret)
}
