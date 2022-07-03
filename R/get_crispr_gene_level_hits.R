#' get_crispr_gene_level_hits
#' 
#' Runs `DESeq2` at the gRNA-level and uses a harmonic mean p-value method to combine the p-values to the gene-level - extracts hits at both the positive and negative ends.
#' 
#' @param data A data frame of normalized counts at the gRNA level.
#' @param grna_column A character string specifying the gRNA ID column.
#' @param gene_column A character string specifying the gene name column.
#' @param group_1 A character vector naming the columns that belong to group 1.
#' @param group_2 A character vector naming the columns that belong to group 2.
#' @param target_fdr Numeric value (0-1) giving the target FDR. Defaults to 0.1.
#' @param dir A character string specifying which end of the fold change distribution to test: `both` (default), `pos`, `neg`.
#' 
#' @return An object of class `delboy_crispr_hits` which is a list containing the following components:
#' * `deseq2_pos`: A data frame containing positive `DESeq2` results at the gRNA level.
#' * `hmp_gene_pos`: A data frame containing positive gene-level results (with p-values derived from the harmonic mean p-value method).
#' * `deseq2_neg`: A data frame containing negative `DESeq2` results at the gRNA level.
#' * `hmp_gene_neg`: A data frame containing negative gene-level results (with p-values derived from the harmonic mean p-value method).
#' @md
#' @export
#' @importFrom dplyr %>% filter mutate select everything
#' @importFrom rlang sym !!
#' @references Wilson, D. J. 2019.The harmonic mean p-value for combining dependent tests. PNAS 116 (4): 1195-1200.
get_crispr_gene_level_hits <- function(data, grna_column, gene_column, group_1, group_2, target_fdr = 0.1, dir = "both"){
  if(!dir %in% c("both","pos","neg"))
    stop(paste("expecting 'dir' to be one of: 'both', 'pos', or 'neg'. Got instead:",dir))
  tryCatch({
    gene_sym <- rlang::sym(gene_column)
    data_deseq2 <- data %>%
      dplyr::select(- (!! gene_sym))
    if(dir %in% c("both","pos")){
      # Positive results.
      des_pos <- delboy::run_deseq2(data_deseq2, group_1, group_2, grna_column, alt_hyp = "greater") %>%
        # Some test results can be NA.
        dplyr::filter(!is.na(pvalue)) %>%
        # Must add a gene column.
        dplyr::mutate(gene = data[match(id, data[,grna_column]), gene_column]) %>%
        dplyr::select(id,gene, dplyr::everything())
      hmp_pos <- delboy::combine_harmonic_mean_pvals(des_pos, "pvalue", "gene", target_fdr)
    }else{
      des_pos <- hmp_pos <- NA
    }
    if(dir %in% c("both","neg")){
      # Negative results.
      des_neg <- delboy::run_deseq2(data_deseq2, group_1, group_2, grna_column, alt_hyp = "less") %>%
        # Some test results can be NA.
        dplyr::filter(!is.na(pvalue)) %>%
        # Must add a gene column.
        dplyr::mutate(gene = data[match(id, data[,grna_column]), gene_column]) %>%
        dplyr::select(id,gene, dplyr::everything())
      hmp_neg <- delboy::combine_harmonic_mean_pvals(des_neg, "pvalue", "gene", target_fdr)
    }else{
      des_neg <- hmp_neg <- NA
    }
    ret <- list(deseq2_pos = des_pos,
                hmp_gene_pos = hmp_pos,
                deseq2_neg = des_neg,
                hmp_gene_neg = hmp_neg)
  },
  error = function(e) stop(paste("unable to extract crispr gene-level hits:",e))
  )
  return(ret)
}
