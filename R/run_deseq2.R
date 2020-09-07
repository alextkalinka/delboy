#' run_deseq2
#'
#' Runs `DESeq2` on normalized, batch-corrected data containing two groups to be contrasted.
#'
#' @param data A data frame containing normalized count data.
#' @param group_1 A character string naming the columns that belong to group 1.
#' @param group_2 A character string naming the columns that belong to group 2.
#' @param gene_column A character string naming the column containing gene names.
#'
#' @return A data frame containing `DESeq2` results.
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom dplyr arrange mutate select everything
#' @export
run_deseq2 <- function(data, group_1, group_2, gene_column){
  tryCatch({
    # Prep data for DESeq2.
    data.m <- as.matrix(data[,c(group_1, group_2)])
    rownames(data.m) <- data[,gene_column]
    data.m <- apply(data.m, 2, round)
    data.m <- apply(data.m, 2, function(x) ifelse(x < 0, 0, x))
    # Run DESeq2.
    tr <- rep("group_1",ncol(data.m))
    tr[which(colnames(data.m) %in% group_2)] <- "group_2"
    tr <- factor(tr)
    dds <- DESeq2::DESeqDataSetFromMatrix(data.m, DataFrame(tr), ~tr)
    dds <- DESeq2::DESeq(dds)
    res <- as.data.frame(DESeq2::results(dds)) %>%
      dplyr::arrange(padj) %>%
      dplyr::mutate(id = rownames(.)) %>%
      dplyr::select(id, dplyr::everything())
  },
  error = function(e) stop(paste("unable to run DESeq2:",e))
  )
  return(res)
}
