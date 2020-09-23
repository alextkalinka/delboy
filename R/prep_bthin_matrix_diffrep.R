#' prep_bthin_matrix_diffrep
#'
#' Prepare the matrix output from `seqgendiff::thin_diff` for use in differential representation analyses by, for example, `DESeq2` or `delboy`.
#'
#' @param data A data frame of normalized (and batch-corrected, if necessary) counts for a set of samples with the true signal controlled (by, for example, batch-correcting it).
#' @param bthin_matrix The matrix output from `seqgendiff::thin_diff`.
#' @param sample_names A character vector giving the sample names from the original input data provided to `seqgendiff::thin_diff` in their original order.
#' @param treat_cols A vector of 1/0 where a 1 indicates a treatment sample in `seqgendiff::thin_diff`.
#' @param gene_column A character string naming the column containing gene names.
#'
#' @return A data frame.
#' @export
#' @importFrom dplyr mutate select everything %>%
#' @importFrom rlang := sym !!
#' @importFrom magrittr %<>%
prep_bthin_matrix_diffrep <- function(data, bthin_matrix, sample_names, treat_cols, gene_column){
  tryCatch({
    sample_names[treat_cols] <- gsub("^(.*?)$","\\1_T",sample_names[treat_cols])
    sample_names[!treat_cols] <- gsub("^(.*?)$","\\1_ctrl",sample_names[!treat_cols])
    colnames(bthin_matrix) <- sample_names
    bthin_matrix %<>%
      as.data.frame() %>%
      dplyr::mutate(!!rlang::sym(gene_column) := data[,gene_column]) %>%
      dplyr::select(!!rlang::sym(gene_column), dplyr::everything())
    bthin_matrix <- bthin_matrix[!is.na(bthin_matrix[,3]),]
  },
  error = function(e) stop(paste("unable to prep bthin matrix for DiffRep analyses:",e))
  )
  return(bthin_matrix)
}
