#' make_ground_truth_data
#' 
#' Add signal to specified num of non-null genes (evenly split across pos and neg cases).
#' 
#' @param data A data frame of counts
#' @param ctrl A character string naming the control sample columns.
#' @param treat A character string naming the treatment sample columns.
#' @param grna_column A character string specifying the gRNA ID column.
#' @param gene_column A character string naming the gene column in `data`.
#' @param num_non_null An integer specifying the number of non-null cases to add signal to.
#' @param treat_sig The samples to be treated as treatment samples.
#' @param data_lfc Output from `delboy::get_crispr_gene_level_hits`. If `NULL` (default) this will be produced within the function.
#' @param data.sc Output from `delboy::prep_val_data`. If `NULL` (default) this will be produced within the function.
#' @param el Output from `delboy::adjust_nonnull_lfc_estimates`. If `NULL` (default) this will be produced within the function.
#' @param method A character string naming the batch correction method: `combat_seq` (default) or `combat`.
#' @param normalize Logical, whether to normalize the data or not (default = `FALSE`).
#' 
#' @return A list.
#' @export
#' @importFrom magrittr %<>%
#' @importFrom dplyr mutate select everything
#' @importFrom rlang sym !!
make_ground_truth_data <- function(data, ctrl, treat, grna_column, gene_column, num_non_null, treat_sig, 
                                   data_lfc = NULL, data.sc = NULL, el = NULL, method="combat_seq", normalize = FALSE){
  # 1. Estimate logFC distribution to sample signal.
  if(normalize)
    data <- delboy::normalize_library_depth_relative(data, NULL)
  if(is.null(data_lfc))
    data_lfc <- delboy::get_crispr_gene_level_hits(data, grna_column, gene_column, ctrl, treat, 0.1)
  # 2. Correct signal in original data.
  if(is.null(data.sc))
    data.sc <- delboy::prep_val_data(data, ctrl, treat, grna_column, gene_column, method=bc_method, normalize_method=NULL)
  # 3. Estimate non-null logFC distribution and adjust to ensure alignment with empirical distribution (from DESeq2).
  if(is.null(el))
    el <- delboy::adjust_nonnull_lfc_estimates(data_lfc, filter_prop = 0.05)
  if(!el$fit_ok) return(el)
  # 4. Combine pos and neg logFC values for plotting.
  lfc_distr <- data.frame(lfc = c(el$non_null.pos.lfc, el$non_null.neg.lfc),
                          dens = c(el$non_null.pos.dens/sum(el$non_null.pos.dens),
                                   el$non_null.neg.dens/sum(el$non_null.neg.dens)),
                          type = c(rep("pos",length(el$non_null.pos.lfc)), rep("neg",length(el$non_null.neg.lfc))),
                          stringsAsFactors = F)
  # 5. Prep data for seqgendiff.
  data.m <- delboy::prep_count_matrix(data.sc, ctrl, treat, grna_column)
  # 6. Sample logFC values.
  lfc_samp <- delboy::sample_lfc_genes_guides(data_lfc$deseq2_pos, num_non_null, gene_column, 
                                              "log2FoldChange", lfc_distr$lfc, lfc_distr$dens)
  
  # 7. Add signal to sampled genes and their guides.
  signal_ls <- delboy::add_signal_val_data(data, data.m, ctrl, treat, treat_sig, lfc_samp$lfc_samp_grna, grna_column)
  
  # 8. Add genes back to ground truth data.
  signal_ls$data.bthin %<>%
    dplyr::mutate(gene = unlist(data[match(sgRNA, unlist(data[,grna_column],use.names = F)),gene_column],use.names = F)) %>%
    dplyr::select(!! rlang::sym(grna_column), gene, dplyr::everything())
  
  return(list(res = data_lfc,
              data.sig_corr = data.sc,
              lfc_distr = lfc_distr,
              lfc_samp = lfc_samp,
              signal_ls = signal_ls,
              treat_samps = treat_sig,
              fit_ok = TRUE))
  
}
