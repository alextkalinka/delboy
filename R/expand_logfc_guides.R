#' expand_logfc_guides
#' 
#' Expands a set of summary per-gene logFC values into within-gene per-guide level logFC values by sampling from original guide-level logFC data.
#' 
#' @param data A data frame containing gene and logFC columns at the gRNA level.
#' @param gene_col A character string naming the gene column in `data`.
#' @param lfc_col A character string naming the logFC column in `data`.
#' @param med_lfc A vector of logFC values (originally sampled from an estimated non-null logFC distribution derived from `data`).
#' @param lfc_half_window A real number providing the size of the half-window around each logFC when defining the sampling population. Defaults to 0.3.
#' 
#' @return A data frame containing a gene and logFC column with sampled logFC values.
#' @export
#' @importFrom magrittr %<>%
#' @importFrom dplyr filter between %>% mutate group_by summarise ungroup
#' @importFrom rlang sym !!
expand_logfc_guides <- function(data, gene_col, lfc_col, med_lfc, lfc_half_window = 0.3){
  tryCatch({
    lfc_sym <- rlang::sym(lfc_col)
    gene_sym <- rlang::sym(gene_col)
    # Median logFCs at gene level.
    mlfc <- data %>%
      dplyr::group_by(!! gene_sym) %>%
      dplyr::summarise(median_lfc = median(!! lfc_sym, na.rm=T)) %>%
      dplyr::ungroup()
    # For each FC expand to within-gene variance in FC.
    ret <- NULL
    for(i in 1:length(med_lfc)){
      tfc <- med_lfc[i]
      if(all(med_lfc > 0)){
        min_lfc <- max(0, tfc - lfc_half_window)
        max_lfc <- max(0, tfc + lfc_half_window)
      }else if(all(med_lfc < 0)){
        min_lfc <- min(0, tfc - lfc_half_window)
        max_lfc <- min(0, tfc + lfc_half_window)
      }else{
        stop("'med_lfc' values must be all positive or all negative")
      }
      # Find a gene with a median logFC within window of the focal logFC.
      td <- mlfc %>%
        dplyr::filter(dplyr::between(median_lfc, min_lfc, max_lfc))
      if(nrow(td) == 0) next
      tgene <- sample(unlist(td[,gene_col], use.names=F), 1)
      # Extract median logFC deviations for this gene.
      tfv <- data %>%
        dplyr::filter(!! gene_sym == tgene) %>%
        # Mean deviation of logFC values for this gene.
        dplyr::mutate(mean_dev_lfc = !! lfc_sym - mean(!! lfc_sym, na.rm = T)) %>%
        dplyr::filter(!is.na(mean_dev_lfc))
      if(nrow(tfv) < 2) next
      # Use median deviations to generate sample of logFC values for this gene.
      ret <- rbind(ret, data.frame(Gene = tgene, num = i, logFC = tfc + tfv$mean_dev_lfc))
    }
    if(is.null(ret) || nrow(ret) == 0) stop("unable to find any logFCs satisfying the criteria")
  },
  error = function(e) stop(paste("'expand_med_logfc_guides': unable to expand logFC values per guide:",e))
  )
  return(ret)
}
