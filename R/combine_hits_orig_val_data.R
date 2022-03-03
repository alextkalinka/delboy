#' combine_hits_orig_val_data
#'
#' Combines hits from original and validation data for comparison purposes.
#'
#' @param orig_hits A data frame as produced by `delboy::get_crispr_gene_level_hits`.
#' @param val_hits A data frame as produced by `delboy::make_delboy_crispr_hit_comparison_table`.
#'
#' @return A data frame.
#' @export
#' @importFrom dplyr %>% select filter mutate
combine_hits_orig_val_data <- function(orig_hits, val_hits){
  tryCatch({
    comb <- rbind(orig_hits %>%
                    dplyr::filter(significant_hit) %>%
                    dplyr::mutate(data_type = "Original"),
                  val_hits %>%
                    dplyr::select(-replicate,-TP,-log10_pvalue,
                                  -log10_baseExpr,-abs_log2FoldChange,-hit_type,-lfc_true_mean) %>%
                    dplyr::mutate(data_type = "Validation"))
  },
  error = function(e) stop(paste("unable to combine hits from original and validation data:",e))
  )
  return(comb)
}
