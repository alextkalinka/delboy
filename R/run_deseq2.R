#' run_deseq2
#'
#' Runs `DESeq2` on normalized, batch-corrected data containing two groups to be contrasted.
#'
#' @param data A data frame containing normalized count data.
#' @param group_1 A character vector naming the columns that belong to group 1.
#' @param group_2 A character vector naming the columns that belong to group 2.
#' @param gene_column A character string naming the column containing gene names.
#'
#' @return A data frame containing `DESeq2` results.
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom dplyr arrange mutate select everything
#' @export
run_deseq2 <- function(data, group_1, group_2, gene_column){
  tryCatch({
    # Prep data for DESeq2.
    data.m <- delboy::prep_count_matrix(data, group_1, group_2, gene_column)
    # Run DESeq2.
    tr <- delboy::make_treat_factor(data.m, group_1, group_2)
    suppressMessages(dds <- DESeq2::DESeqDataSetFromMatrix(data.m, as.data.frame(tr), ~tr))
    suppressMessages(dds <- DESeq2::DESeq(dds))
    res <- as.data.frame(DESeq2::results(dds)) %>%
      dplyr::arrange(padj) %>%
      dplyr::mutate(id = rownames(.)) %>%
      dplyr::select(id, dplyr::everything())
  },
  error = function(e) stop(paste("unable to run DESeq2:",e))
  )
  return(res)
}
