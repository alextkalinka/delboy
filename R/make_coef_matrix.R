#' make_coef_matrix
#'
#' Make a coefficient matrix for use by `seqgendiff`.
#'
#' @param data A data frame of normalized (and batch-corrected, if necessary) counts for a set of samples.
#' @param lfc_samp A named vector of logFC values where names are gene names.
#' @param gene_column A character string naming the column containing gene names.
#'
#' @return A coefficient matrix.
#' @export
#' @importFrom dplyr mutate select
make_coef_matrix <- function(data, lfc_samp, gene_column){
  tryCatch({
    coef_mat <- data %>%
      dplyr::mutate(lfc = ifelse(!!sym(gene_column) %in% names(lfc_samp),
                                 lfc_samp[match(!!sym(gene_column), names(lfc_samp))], 0)) %>%
      dplyr::select(lfc) %>%
      as.matrix
  },
  error = function(e) stop(paste("unable to create coefficient matrix:",e))
  )
  return(coef_mat)
}
