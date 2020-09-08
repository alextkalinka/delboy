#' run_delboy
#'
#' Performs a differential-representation analysis using an elastic-net logistic regression approach for (normalized) count data that is split into two groups.
#'
#' @param data A data frame containing count data for two different groups and their replicates. Can be a path to a file. If the count data needs to be normalized, indicate this using the `normalize` argument.
#' @param group_1 A character string naming the columns that belong to group 1.
#' @param group_2 A character string naming the columns that belong to group 2.
#' @param normalize Character string naming the count normalization method. Can be one of `median_ratio`, `relative`, or `NULL` (indicates the data has already been normalized). If the input is RNAseq data, ideally the counts will already have been normalized to Transcripts Per Million (TPM), using, for example, the bias-aware quantification methods employed by `salmon` (Patro et al. 2017).
#' @param filter_cutoff A numerical value indicating the cutoff below which (summed across all replicates) a gene (or gRNA) will be removed from the data. For example, to keep only genes with more than 1 TPM on average across both groups, set the cutoff to 10 if there are 10 replicates in total.
#' @param gene_column A character string naming the column containing gene names.
#' @param batches A character vector identifying the batch structure of the experiment. The length must equal `length(group_1) + length(group_2)`. If `NULL`, there are no batches. Defaults to `NULL`.
#' @param batch_corr_method A character string naming the batch correction method. Can be one of `combat_np` (non-parametric) or `combat_p` (parametric). Defaults to `combat_np`. Ignored if `batches` is `NULL`.
#' @param crispr Logical - is the data from a CRISPR pooled screen. Defaults to `FALSE`.
#' @param grna_column A character string naming the column containing sgRNA IDs. Defaults to `sgRNA`. Ignored if `crispr` set to `FALSE`.
#'
#' @return An object of class `delboy`.
#' @importFrom dplyr select filter
#' @importFrom magrittr %<>%
#' @export
#' @references
#' Kalinka, A. T. 2020.
#' Patro, R. et al. 2017. Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods 14: 417-419.
run_delboy <- function(data, group_1, group_2, normalize, filter_cutoff, gene_column,
                       batches = NULL, batch_corr_method = "combat_np",
                       crispr = FALSE, grna_column = "sgRNA"){
  ### 1. Read data.
  if(is.character(data)){
    tryCatch(
      data <- read.delim(data, stringsAsFactors = F),
      error = function(e) stop(paste("unable to read data file",data))
    )
  }else if(!is.data.frame(data)){
    stop(paste("expected 'data' to be a data frame, got instead an object of class",class(data)))
  }

  ### 2A. Sanity checks.
  if(any(!group_1 %in% colnames(data)))
    stop(paste("unable to find the following group 1 columns:",setdiff(group_1,colnames(data))))
  if(any(!group_2 %in% colnames(data)))
    stop(paste("unable to find the following group 2 columns:",setdiff(group_2,colnames(data))))
  if(!gene_column %in% colnames(data))
    stop(paste("unable to find the gene column",gene_column))
  if(crispr){
    if(!grna_column %in% colnames(data))
      stop(paste("unable to find the gRNA column",grna_column))
  }

  ### 2B. Remove any irrelevant columns.
  if(crispr){
    data <- data[,c(gene_column, grna_column, group_1, group_2)]
  }else{
    data <- data[,c(gene_column, group_1, group_2)]
  }

  ### 3. Normalize counts.
  if(!is.null(normalize)){

  }

  ### 4. Filter low count data prior to batch correction.
  tryCatch(data <- data[rowSums(data[,c(group_1, group_2)]) > filter_cutoff,],
           error = function(e) stop(paste("unable to filter low-count data:",e))
  )
  if(nrow(data) < 200)
    stop(paste("only",nrow(data),"rows remain after filtering: consider using a less stringent 'filter_cutoff' value - used",filter_cutoff))

  ### 5. Batch correction.
  if(!is.null(batches)){

  }

  ### 6. Run DESeq2 on the original dataset.
  if(crispr){
    deseq2_res <- delboy::run_deseq2(data, group_1, group_2, grna_column) %>%
      dplyr::left_join(data, by = c(id = grna_column))
  }else{
    deseq2_res <- delboy::run_deseq2(data, group_1, group_2, gene_column) %>%
      dplyr::left_join(data, by = c(id = gene_column))
  }

  ### 7. Estimate parameters for validation.
  ## 7A. Number of non-null cases.
  non.null <- delboy::estimate_number_non_nulls(deseq2_res$pvalue)

  ## 7B. Estimate non-null logFC distribution.
  lfdr.lfc <- locfdr(deseq2_res$log2FoldChange)
  non_null.dens <- lfdr.lfc$mat[1:which(lfdr.lfc$mat[,11]==0)[1],11]
  non_null.lfc <- lfdr.lfc$mat[1:which(lfdr.lfc$mat[,11]==0)[1],1]

  ### . Prep data for DR analysis.


  ### . Elastic-net logistic regression to identify differentially-represented genes or gRNAs.

}
