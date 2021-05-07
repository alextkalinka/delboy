#' evaluate_performance_rnaseq_calls
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
#' @param alpha The elastic-net regression penalty, between 0 and 1.
#'
#' @return An object of class `delboy_performance`.
#' @export
#' @importFrom seqgendiff thin_diff
#' @importFrom dplyr select
#' @importFrom rlang sym !!
evaluate_performance_rnaseq_calls <- function(data, group_1, group_2, gene_column, max.iter,
                                              num_non_null, lfc, lfc_dens, alpha){
  tryCatch({
    # 1. Prep for seqgendiff.
    data.m <- delboy::prep_count_matrix(data, group_1, group_2, gene_column)

    # 2. All combinations of treatment samples (preserving the number of treatment samples in original data [length(group_2)]).
    all_treat_comb <- delboy::all_combinations_treat_samples(c(group_1, group_2), length(group_2))
    if(is.null(max.iter)){
      num_treat_comb <- ncol(all_treat_comb)
    }else{
      all_treat_comb <- all_treat_comb[,1:max.iter]
    }
    
    all_val_hits <- NULL
    all_val_perf <- NULL
    # Multiple val samples to improve SVM estimate.
    for(i in 1:ncol(all_treat_comb)){
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

      # 9. Run DESeq2 on bthin data.
      group_1.v <- colnames(data.bthin %>%
                              dplyr::select(- !!rlang::sym(gene_column)))[!as.logical(c(design_mat))]
      group_2.v <- colnames(data.bthin %>%
                              dplyr::select(- !!rlang::sym(gene_column)))[as.logical(c(design_mat))]
      deseq2_res <- delboy::run_deseq2(data.bthin, group_1.v, group_2.v, gene_column)

      # 10. Prep data for Elastic-net logistic regression.
      data.elnet <- delboy::prep_elnet_data(data.bthin, group_1.v, group_2.v, gene_column)

      # 11. Run Elastic-net logistic regression on bthin data.
      elnet.lr <- delboy::run_elnet_logistic_reg(as.matrix(data.elnet[,3:ncol(data.elnet)]),
                                               factor(data.elnet$treat),
                                               alpha = alpha)

      # 12. Extract performance statistics.
      perf_stats <- delboy::perf_stats_rnaseq(elnet.lr, deseq2_res, lfc_samp)

      # 13. Collate TP, FN, and FP into a data frame to aid comparisons.
      delboy_hit_df <- delboy::make_delboy_hit_comparison_table(elnet.lr,
                                                                deseq2_res,
                                                                lfc_samp)
      all_val_hits <- rbind(all_val_hits, delboy_hit_df)
      all_val_perf <- rbind(all_val_perf, perf_stats)
    }
    # 14. Performance stats.
    pstats_summ <- delboy::calculate_perf_stats(all_val_perf)
    
    # 15. Is the validation data sufficient for finding SVM decision boundary?
    if(sum(all_val_hits$hit_type == "True_Positive") == 0)
      stop("* no true positives in validation data: unable to validate algorithms *")
    
    if(sum(all_val_hits$hit_type == "False_Positive") < 10){
      cat("* less than 10 False Positive cases in validation data: using lowest 75% of False Negatives *\n")
      use_fn <- TRUE
    }else{
      use_fn <- FALSE
    }
 
    # 16. SVM for false positive classification.
    svm_validation <- delboy::svm_false_positive_classification(all_val_hits, use_fn = use_fn)

    # 17. Build return object of class 'delboy_performance'.
    ret <- list(lfc_samp = lfc_samp,
                data.bthin = data.bthin,
                elnet_lr_res = elnet.lr,
                deseq2_res = deseq2_res,
                delboy_hit_table = all_val_hits,
                performance_stats = pstats_summ,
                svm_validation = svm_validation)
    class(ret) <- "delboy_performance"
  },
  error = function(e) stop(paste("unable to evaluate performance of delboy:",e))
  )
  return(ret)
}
