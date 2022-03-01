#' evaluate_performance_crispr_calls
#'
#' Evaluates the performance of `delboy` on CRISPR pooled data by controlling for real signal, and adding known signal (using `seqgendiff`'s binomial-thinning approach) for a sampled number of genes from a logFC distribution with both the number and distribution chosen to match as closely as possible the signal in the real data.
#'
#' @param data A data frame of normalized counts (at the gRNA level) for a set of samples with the true signal controlled.
#' @param data_lfc A data frame of logFC estimates for `data`.
#' @param group_1 A character vector naming the columns that belong to group 1.
#' @param group_2 A character vector naming the columns that belong to group 2.
#' @param gene_column A character string naming the column containing gene names.
#' @param grna_column A character string naming the column containing gRNA IDs.
#' @param lfc_column A character string naming the logFC column in `data_lfc`.
#' @param max.iter An integer value indicating the maximum number of validation sample combinations..
#' @param num_non_null An integer value indicating the number of genes to add signal to.
#' @param lfc A vector of logFC values for non-null cases.
#' @param lfc_dens A vector of density estimates for the logFC values given in `lfc`.
#' @param alt_hyp A character string specifying the alternative hypothesis. Can be one of: `greater`, or `less` (see `DESeq2::results` documentation).
#' @param target_fdr Numeric value (0-1) giving the target FDR. Defaults to 0.1.
#'
#' @return An object of class `delboy_performance_crispr`.
#' @export
#' @importFrom seqgendiff thin_diff
#' @importFrom dplyr %>% select filter mutate group_by ungroup n
#' @importFrom magrittr %<>%
#' @importFrom rlang sym !!
#' @importFrom progress progress_bar
evaluate_performance_crispr_calls <- function(data, data_lfc, group_1, group_2, gene_column, grna_column, lfc_column,
                                              max.iter, num_non_null, lfc, lfc_dens, 
                                              alt_hyp, target_fdr = 0.1){
  if(!alt_hyp %in% c("less","greater")) stop("'alt_hyp' must be one of 'less' or 'greater'")
  tryCatch({
    # 1. Prep for seqgendiff.
    data.m <- delboy::prep_count_matrix(data, group_1, group_2, grna_column)

    # 2. All combinations of treatment samples (preserving the number of treatment samples in original data [length(group_2)]).
    all_treat_comb <- delboy::all_combinations_treat_samples(c(group_1, group_2), length(group_2))
    all_treat_comb <- delboy::max_diff_val_samples(all_treat_comb, max.iter)
    num_val_combs <- ncol(all_treat_comb)

    all_val_hits <- NULL
    all_val_perf <- NULL
    deseq <- list()
    counts.bthin <- list()
    # Set up progress bar.
    pb <- progress::progress_bar$new(
      format = "  validating [:bar] :percent time left: :eta",
      total = num_val_combs, clear = FALSE, width = 60)
    pb$tick(0)
    
    # Multiple val samples to improve performance estimate.
    for(i in 1:num_val_combs){
      pb$tick()
      # 3. Sample logFC values for num_non_null cases.
      lfc_ls <- delboy::sample_lfc_genes_guides(data_lfc, num_non_null, gene_column, lfc_column,
                                                     lfc, lfc_dens)
      lfc_samp <- lfc_ls$lfc_samp
      lfc_samp_grna <- lfc_ls$lfc_samp_grna
      
      # 4. Add logFC signal to signal-corrected data.
      treat_samps <- all_treat_comb[,i]
      signal_ls <- delboy::add_signal_val_data(data.m, group_1, group_2, treat_samps, 
                                               lfc_samp_grna, grna_column)

      # 5. Run DESeq2 on bthin data.
      deseq2_res <- delboy::run_deseq2(data.bthin, group_1.v, group_2.v, grna_column, alt_hyp = alt_hyp) %>%
        # Some test results can be NA.
        dplyr::filter(!is.na(pvalue)) %>%
        # Must add a gene column.
        dplyr::mutate(gene = data[match(id, data[,grna_column]),gene_column])

      # 6. Combine gRNA p-values using harmonic meanp approach.
      comb_pvals <- delboy::combine_harmonic_mean_pvals(deseq2_res, "pvalue", "gene", target_fdr = 0.1)

      # 7. Collate TP, FN, and FP into a data frame to aid comparisons.
      delboy_hit_df <- delboy::make_delboy_crispr_hit_comparison_table(comb_pvals, lfc_samp) %>%
        dplyr::mutate(replicate = i)
      
      all_val_hits <- rbind(all_val_hits, delboy_hit_df)
      deseq[[i]] <- deseq2_res
      counts.bthin[[i]] <- data.bthin
    }

    # 8. Estimate aggregate FDR across samples.
    summ_stats <- all_val_hits %>%
      dplyr::summarise(tot_tp = sum(all_val_hits$hit_type == "True_Positive"),
                       tot_fp = sum(all_val_hits$hit_type == "False_Positive"),
                       tot_fn = sum(all_val_hits$hit_type == "False_Negative"),
                       Sensitivity = 100*tot_tp/(tot_tp + tot_fn),
                       FDR = 100*tot_fp/(tot_tp + tot_fp),
                       FDR = ifelse(is.nan(FDR),0,FDR))
    
    # 9. Is FDR greater than 'target_fdr'? If so, attempt to find threshold value of a summary lfc metric to reduce FDR below target.
    metr_ls <- NA
    if(summ_stats$FDR > target_fdr && summ_stats$tot_tp > 5){
      metr_ls <- delboy::infer_metric_for_thr(all_val_hits)
    }
    fdr_thr <- ifelse(!is.na(metr_ls), 
                      delboy::calc_fdr_threshold(all_val_hits, metr_ls$metric, "hit_type", target_fdr), NA)
    
    # 10. Re-calculate summary stats after applying threshold.
    
    
    # 11. Calculate AUPrRc.
    

    # 12. Build return object of class 'delboy_performance'.
    ret <- list(lfc_samp = lfc_ls,
                data.bthin = counts.bthin,
                deseq2_res = deseq,
                hits = all_val_hits,
                all_treat_combinations = all_treat_comb,
                alt_hypothesis = alt_hyp)
    class(ret) <- "delboy_perf_crispr"
  },
  error = function(e) stop(paste("unable to evaluate performance of delboy crispr:",e))
  )
  return(ret)
}
