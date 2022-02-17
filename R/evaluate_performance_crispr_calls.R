#' evaluate_performance_crispr_calls
#'
#' Evaluates the performance of `delboy` on CRISPR pooled data by controlling for real signal, and adding known signal (using `seqgendiff`'s binomial-thinning approach) for a sampled number of genes from a logFC distribution with both the number and distribution chosen to match as closely as possible the signal in the real data.
#'
#' @param data A data frame of normalized (and batch-corrected, if necessary) counts (at the gRNA level) for a set of samples with the true signal controlled (by, for example, batch-correcting it).
#' @param data_lfc A data frame of logFC estimates for `data`.
#' @param group_1 A character vector naming the columns that belong to group 1.
#' @param group_2 A character vector naming the columns that belong to group 2.
#' @param gene_column A character string naming the column containing gene names.
#' @param grna_column A character string naming the column containing gRNA IDs.
#' @param lfc_column A character string naming the logFC column in `data_lfc`.
#' @param max.iter An integer value indicating the maximum number of validation samples.
#' @param num_non_null An integer value indicating the number of genes to add signal to.
#' @param lfc A vector of logFC values for non-null cases.
#' @param lfc_dens A vector of density estimates for the logFC values given in `lfc`.
#' @param alt_hyp A character string specifying the alternative hypothesis. Can be one of: `greater`, or `less` (see `DESeq2::results` documentation).
#' @param target_fdr Numeric value (0-1) giving the target FDR. Defaults to 0.1.
#'
#' @return An object of class `delboy_performance_crispr`.
#' @export
#' @importFrom seqgendiff thin_diff
#' @importFrom dplyr select filter mutate group_by ungroup n
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
    tot_val_combs <- ncol(all_treat_comb)
    if(!is.null(max.iter)){
      num_val_combs <- min(max.iter, tot_val_combs)
      # Try to use maximally different sets of samples.
      val_inds <- seq(1,tot_val_combs, by = ceiling(tot_val_combs/num_val_combs))
      all_treat_comb <- all_treat_comb[,val_inds]
    }else{
      num_val_combs <- ncol(all_treat_comb)
    }
    
    all_val_hits <- NULL
    all_val_perf <- NULL
    deseq <- list()
    counts.bthin <- list()
    # Set up progress bar.
    pb <- progress::progress_bar$new(
      format = "  validating [:bar] :percent time left: :eta",
      total = num_val_combs, clear = FALSE, width = 60)
    pb$tick(0)
    
    # Multiple val samples to improve SVM estimate.
    for(i in 1:num_val_combs){
      pb$tick()
      # 3. Sample logFC values for num_non_null cases.
      lfc_samp <- sample(lfc, num_non_null, prob = lfc_dens/sum(lfc_dens), replace = T)
      # We want the variance of logFC values to reflect the sampled logFC values.
      lfc_samp_df <- delboy::expand_logfc_guides(data_lfc, gene_column, lfc_column, lfc_samp)

      # 4. Sample genes and guides to add signal to.
      genes_signal <- sample(unique(unlist(data[,gene_column],use.names = F)), 
                             length(unique(lfc_samp_df$Gene)), replace = F)
      names(lfc_samp) <- genes_signal
      names(genes_signal) <- 1:length(genes_signal)
      # Map gRNA IDs to genes and filter duplicates (sampled genes could have different numbers of guides).
      lfc_samp_df %<>%
        dplyr::mutate(gene = genes_signal[match(num, names(genes_signal))]) %>%
        dplyr::group_by(gene) %>%
        dplyr::mutate(sgRNA = data_lfc$id[data_lfc$gene == gene[1]][1:dplyr::n()]) %>%
        dplyr::ungroup() %>%
        dplyr::filter(!duplicated(sgRNA))

      lfc_samp_grna <- lfc_samp_df$logFC
      names(lfc_samp_grna) <- lfc_samp_df$sgRNA

      # 5. Create coefficient matrix for seqgendiff.
      coef_mat <- delboy::make_coef_matrix(data, lfc_samp_grna, grna_column)

      # 6. Create design matrix for seqgendiff.
      treat_samps <- all_treat_comb[,i]
      design_mat <- delboy::make_design_matrix(group_1, group_2, treat_samps)

      # 7. Add signal using seqgendiff's binomial-thinning approach.
      thout <- seqgendiff::thin_diff(mat = data.m,
                                     design_fixed = design_mat,
                                     coef_fixed = coef_mat)

      # 8. Prep bthin matrix for use in DiffExp analyses.
      data.bthin <- delboy::prep_bthin_matrix_diffrep(data, thout$mat,
                                                      colnames(data.m),
                                                      as.logical(c(design_mat)),
                                                      grna_column)

      # 9. Run DESeq2 on bthin data.
      group_1.v <- colnames(data.bthin %>%
                              dplyr::select(- !!rlang::sym(grna_column)))[!as.logical(c(design_mat))]
      group_2.v <- colnames(data.bthin %>%
                              dplyr::select(- !!rlang::sym(grna_column)))[as.logical(c(design_mat))]
      
      deseq2_res <- delboy::run_deseq2(data.bthin, group_1.v, group_2.v, grna_column, alt_hyp = alt_hyp) %>%
        # Some test results can be NA.
        dplyr::filter(!is.na(pvalue)) %>%
        # Must add a gene column.
        dplyr::mutate(gene = data[match(id, data[,grna_column]),gene_column])

      # 10. Combine gRNA p-values using harmonic meanp approach.
      comb_pvals <- delboy::combine_harmonic_mean_pvals(deseq2_res, "pvalue", "gene", target_fdr = 0.1)

      # 11. Collate TP, FN, and FP into a data frame to aid comparisons.
      delboy_hit_df <- delboy::make_delboy_crispr_hit_comparison_table(comb_pvals, lfc_samp)
      
      all_val_hits <- rbind(all_val_hits, delboy_hit_df)
      deseq[[i]] <- deseq2_res
      counts.bthin[[i]] <- data.bthin
    }
    return(list(hits = all_val_hits,
                deseq2 = deseq,
                data = counts.bthin))

    # 12. Estimate aggregate FDR across samples.
    tot_tp <- sum(all_val_hits$hit_type == "True_Positive")
    tot_fp <- sum(all_val_hits$hit_type == "False_Positive")
    fdr_est <- 100*tot_fp/(tot_tp + tot_fp)
    if(is.nan(fdr_est)) fdr_est <- 0
    
    # 13. Is FDR greater than 'target_fdr'? If so, attempt to find threshold value of a summary lfc metric to reduce FDR below target.
    if(fdr_est > target_fdr){
      if(tot_tp > 0){
        
      }
    }else{
      fdr_correction <- NA
    }
    


    # 14. Build return object of class 'delboy_performance'.
    ret <- list(lfc_samp = lfc_samp,
                lfc_samp_df = lfc_samp_df,
                data.bthin = data.bthin,
                deseq2_res = deseq2_res,
                delboy_hit_table = all_val_hits,
                svm_validation = svm_validation,
                num_val_combinations = num_val_combs,
                all_treat_combinations = all_treat_comb)
    class(ret) <- "delboy_performance_crispr"
  },
  error = function(e) stop(paste("unable to evaluate performance of delboy crispr:",e))
  )
  return(ret)
}
