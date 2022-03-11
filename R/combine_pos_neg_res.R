#' combine_pos_neg_res
#' 
#' Combines positive and negative gene-level CRISPR results.
#' 
#' @param data_pos A data frame of positive hits as produced by `delboy::get_crispr_gene_level_hits`.
#' @param data_neg A data frame of negative hits as produced by `delboy::get_crispr_gene_level_hits`.
#' 
#' @return A data frame.
#' @export
#' @importFrom dplyr %>% full_join mutate rename select arrange
combine_pos_neg_res <- function(data_pos, data_neg){
  tryCatch({
    both <- data_pos %>%
      dplyr::full_join(data_neg, by = "gene", suffix = c(".pos",".neg")) %>%
      # If a hit is in both pos and neg data then set to F (likely due to issue such as clonal outgrowth).
      dplyr::mutate(significant_hit.pos = ifelse((significant_hit.pos & significant_hit.neg),F,significant_hit.pos),
                    significant_hit.neg = ifelse((significant_hit.pos & significant_hit.neg),F,significant_hit.neg)) %>%
      dplyr::rename(mean_log2FoldChange = mean_log2FoldChange.pos,
                    median_log2FoldChange = median_log2FoldChange.pos,
                    mean_baseMean = mean_baseMean.pos,
                    sd_log2FoldChange = sd_log2FoldChange.pos) %>%
      dplyr::select(-mean_log2FoldChange.neg, -median_log2FoldChange.neg,
                    -mean_baseMean.neg, -sd_log2FoldChange.neg) %>%
      dplyr::arrange(pvalue.harmonic_mean.pos, gene)
  },
  error = function(e) stop(paste("unable to combine positive and negative CRISPR results:",e))
  )
  return(both)
}
