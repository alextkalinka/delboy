#' evaluate_performance_boostx
#'
#' Evaluates the performance of `delboy` on RNAseq data by comparison with `DESeq2` output on the original input data for an experiment, controlling for real signal, and adding known signal (using `seqgendiff`'s binomial-thinning approach) for a sampled number of genes from a logFC distribution with both the number and distribution chosen to match as closely as possible the signal in the real data.
#'
#' @param data A data frame of normalized (and batch-corrected, if necessary) counts for a set of samples with the true signal controlled (by, for example, batch-correcting it).
#' @param group_1 A character vector naming the columns that belong to group 1.
#' @param group_2 A character vector naming the columns that belong to group 2.
#' @param gene_column A character string naming the column containing gene names.
#' @param max.iter An integer value indicating the maximum number of validation samples.
#' @param num_non_null An integer value indicating the number of genes to add signal to.
#' @param lfc A vector of logFC values for non-null cases.
#' @param lfc_dens A vector of density estimates for the logFC values given in `lfc`.
#' @param target_fdr A numerical value (0-1) indicating the target FDR.
#' @param algorithm A chracter string naming the algorithm to use (lower-case only): `deseq2`, `mageck`.
#'
#' @return An object of class `boostx_performance`.
#' @export
#' @importFrom seqgendiff thin_diff
#' @importFrom dplyr select
#' @importFrom rlang sym !!
#' @importFrom progress progress_bar
evaluate_performance_boostx <- function(data, group_1, group_2, gene_column, max.iter,
                                           num_non_null, lfc, lfc_dens, target_fdr, algorithm){
  tryCatch({
    # 1. Prep for seqgendiff.
    data.m <- delboy::prep_count_matrix(data, group_1, group_2, gene_column)
    
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
    
    perf <- deg_res <- NULL
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
      
      # 4. Sample genes to add signal to.
      genes_signal <- sample(data[,gene_column], num_non_null, replace = F)
      names(lfc_samp) <- genes_signal
      
      # 5. Create coefficient matrix for seqgendiff.
      coef_mat <- delboy::make_coef_matrix(data, lfc_samp, gene_column)
      
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
                                                      gene_column)
      
      # 9. Run DEG algorithm.
      group_1.v <- colnames(data.bthin %>%
                              dplyr::select(- !!rlang::sym(gene_column)))[!as.logical(c(design_mat))]
      group_2.v <- colnames(data.bthin %>%
                              dplyr::select(- !!rlang::sym(gene_column)))[as.logical(c(design_mat))]
      
      switch(
        type,
        deseq2 = tdeg_res <- delboy::run_deseq2(data.bthin, group_1.v, group_2.v, gene_column) %>%
          dplyr::mutate(abs_log2FoldChange = abs(log2FoldChange)),
        mageck = tdeg_res <- delboy::run_mageck()
      )
      
      # Prelim mark genes with/without added signal prior to knowing the final p-value and abs(LFC) thresholds.
      deg_res <- rbind(deg_res,
                       tdeg_res %>%
                         dplyr::mutate(signal = id %in% genes_signal,
                                       val_repl = i))
      
      # 10. Calculate performance.
      perf <- rbind(perf, delboy::calc_perf_pval_windows(tdeg_res, "pvalue", "id", "abs_log2FoldChange", genes_signal))
    }
    
    # 11. Determine the p-value and abs(LFC) thresholds to apply to original results.
    pval_lfc_thresh <- delboy::find_pval_target_fdr(perf, 100*target_fdr)
    
    # 12. Apply thresholds to validation data.
    deg_res <- delboy::apply_thresholds_val(deg_res, pval_lfc_thresh)
    
    # 13. Performance stats.

    
    
    # 14. Build return object of class 'delboy_performance'.
    ret <- list(lfc_samp = lfc_samp,
                data.bthin = data.bthin,
                deg_res = deg_res,
                pval_lfc_thresholds = pval_lfc_thresh,
                performance = perf,
                num_val_combinations = num_val_combs,
                all_treat_combinations = all_treat_comb)
    class(ret) <- "boostx_performance"
  },
  error = function(e) stop(paste("unable to evaluate performance of boostX:",e))
  )
  return(ret)
}
