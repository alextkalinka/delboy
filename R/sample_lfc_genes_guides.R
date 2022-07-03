# Helper functions.
# Map gRNA IDs to sampled genes.
.map_guide_ids <- function(data_to, data_from, genes_sig){
  data_to %<>%
    dplyr::mutate(gene = genes_sig[match(num, names(genes_sig))]) %>%
    dplyr::group_by(gene) %>%
    dplyr::mutate(sgRNA = data_from$id[data_from$gene == gene[1]][1:dplyr::n()]) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!duplicated(sgRNA) & !is.na(sgRNA))
  return(data_to)
}


#' sample_lfc_genes_guides
#'
#' For performance estimation, this function samples logFC values, as well as genes and guides to add this signal to.
#'
#' @param data_lfc A data frame containing logFC estimates at the gRNA level.
#' @param num_non_null An integer value indicating the number of genes to add signal to.
#' @param gene_column A character string naming the column containing gene names.
#' @param grna_column A character string naming the column containing gRNA IDs.
#' @param lfc_column A character string naming the logFC column in `data_lfc`.
#' @param lfc A vector of logFC values for non-null cases.
#' @param lfc_dens A vector of density estimates for the logFC values given in `lfc`.
#'
#' @return A list containing the following components:
#' * `lfc_samp`: A named vector of logFC values, where names are gene names.
#' * `lfc_samp_grna`: A named vector of logFC values, where names are gRNA IDs.
#' * `lfc_samp_df`: A data frame of sampled logFC values at the gRNA level.
#' @md
#' @export
sample_lfc_genes_guides <- function(data_lfc, num_non_null, gene_column, lfc_column, lfc, lfc_dens){
  tryCatch({
    # 1. Sample logFC values for num_non_null cases.
    lfc_samp <- sample(lfc, num_non_null, prob = lfc_dens/sum(lfc_dens), replace = T)
    # 2. We want the variance of logFC values to reflect the sampled logFC values.
    lfc_samp_df <- delboy::expand_logfc_guides(data_lfc, gene_column, lfc_column, lfc_samp)
    
    # 3. Sample genes and guides to add signal to.
    genes_signal <- sample(unique(unlist(data_lfc[,gene_column],use.names = F)), 
                           length(unique(lfc_samp_df$num)), replace = F)
    names(lfc_samp) <- genes_signal
    names(genes_signal) <- 1:length(genes_signal)
    
    # 4. Map gRNA IDs to genes and filter duplicates (sampled genes could have different numbers of guides).
    lfc_samp_df <- .map_guide_ids(lfc_samp_df, data_lfc, genes_signal)
    
    lfc_samp_grna <- lfc_samp_df$logFC
    names(lfc_samp_grna) <- lfc_samp_df$sgRNA
  },
  error = function(e) stop(paste("unable to sample logFC values, genes, and guides:",e))
  )
  return(list(lfc_samp = lfc_samp,
              lfc_samp_grna = lfc_samp_grna,
              lfc_samp_df = lfc_samp_df))
}
