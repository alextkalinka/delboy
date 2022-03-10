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
#' @param max.iter An integer value indicating the maximum number of validation sample combinations (default = 3). `NULL` indicates all sample combinations should be taken. `NULL` indicates that the validation step should be skipped.
#' @param filter_prop A numerical value (0-1) to filter the lowest abundance genes when estimating the non-null logFC distribution for performance estimation. Defaults to 0.05.
#'
#' @return An object of class `delboy_crispr`.
#' @importFrom dplyr %>% arrange desc filter summarise n mutate full_join
#' @importFrom rlang sym !!
#' @importFrom utils read.delim packageVersion
#' @export
run_delboy_crispr <- function(data, controls, treatments, grna_column, gene_column, 
                              target_fdr = 0.1, normalize_method = "relative", max.iter = 3, filter_prop = 0.05){
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
    
    ### 5. Estimate performance and learn a summary logfc metric threshold to reduce FDR below the target FDR.
    if(!is.null(max.iter)){
      .print_progress("Estimating algorithm performance...")
      # Correct signal associated with known sample groupings.
      data.sc <- delboy::prep_val_data(data, controls, treatments, "sgRNA", "gene", NULL)
      # Estimate non-null logFC distribution and adjust to ensure alignment with empirical distribution (from DESeq2).
      el <- delboy::adjust_nonnull_lfc_estimates(res, filter_prop = filter_prop)
      # Positive end.
      perf_eval.pos <- delboy::evaluate_performance_crispr_calls(data.sc, res$deseq2_pos, controls, treatments,
                                                                 "gene", "sgRNA", "log2FoldChange",
                                                                 max.iter, 200, el$non_null.pos.lfc, el$non_null.pos.dens,
                                                                 "greater")
      # Negative end.
      perf_eval.neg <- delboy::evaluate_performance_crispr_calls(data.sc, res$deseq2_neg, controls, treatments,
                                                                 "gene", "sgRNA", "log2FoldChange",
                                                                 max.iter, 200, el$non_null.neg.lfc, el$non_null.neg.dens,
                                                                 "less")
      # Combine original and validation data.
      orig_val_hits.pos <- delboy::combine_hits_orig_val_data(res$hmp_gene_pos, perf_eval.pos$hits)
      orig_val_hits.neg <- delboy::combine_hits_orig_val_data(res$hmp_gene_neg, perf_eval.neg$hits)
    }else{
      el <- perf_eval.pos <- perf_eval.neg <- orig_val_hits.pos <- orig_val_hits.neg <- NA
    }
    
    ### 6. Mark up any predicted FPs using any logfc FDR thresholds from the validation data.
    .print_progress("Finishing up...")
    # Positive.
    if(is.list(perf_eval.pos)){
      if(!is.na(perf_eval.pos$lfc_fdr_threshold)){
        res$hmp_gene_pos <- delboy::mark_up_FPs(res$hmp_gene_pos, perf_eval.pos$metr_fdr_thr$metric,
                                                perf_eval.pos$lfc_fdr_threshold, "greater")
      }
    }
    # Negative.
    if(is.list(perf_eval.neg)){
      if(!is.na(perf_eval.neg$lfc_fdr_threshold)){
        res$hmp_gene_neg <- delboy::mark_up_FPs(res$hmp_gene_neg, perf_eval.neg$metr_fdr_thr$metric,
                                                perf_eval.neg$lfc_fdr_threshold, "less")
      }
    }
    # Build single results data frame for both pos and neg.
    res_both <- delboy::combine_pos_neg_res(res$hmp_gene_pos, res$hmp_gene_neg)
    
    # Combine pos and neg non-null logFC distribution estimates for plotting.
    if(is.list(el)){
      lfc_distr <- data.frame(lfc = c(el$non_null.pos.lfc, el$non_null.neg.lfc),
                              dens = c(el$non_null.pos.dens/sum(el$non_null.pos.dens),
                                       el$non_null.neg.dens/sum(el$non_null.neg.dens)),
                              type = c(rep("pos",length(el$non_null.pos.lfc)), rep("neg",length(el$non_null.neg.lfc))),
                              stringsAsFactors = F)
    }else{
      lfc_distr <- NA
    }
    
    ### 7. Build return object.
    ret <- list(data.norm = data,
                norm.method = normalize_method,
                results_gene = res_both,
                results_all = res,
                perf_esitmation = !is.na(max.iter),
                logfc_nonnull_distr = lfc_distr,
                target_fdr = target_fdr,
                max.iter = max.iter,
                perf_estimate_pos = perf_eval.pos,
                perf_estimate_neg = perf_eval.neg,
                comb_orig_valid_hits.pos = orig_val_hits.pos,
                comb_orig_valid_hits.neg = orig_val_hits.neg)
    class(ret) <- "delboy_crispr"
    
  },
  error = function(e) stop(paste("unable to run delboy crispr pipeline:",e))
  )
  return(ret)
}
