#' prep_elnet_data
#'
#' Prepare data for use in an Elastic-net logistic regression.
#'
#' @param data A data frame produced by `delboy::prep_bthin_matrix_diffrep`.
#' @param gene_column A character string naming the column containing gene names.
#'
#' @return A data frame.
#' @export
#' @importFrom dplyr rowwise mutate ungroup sym
#' @importFrom tidyr pivot_longer pivot_wider
prep_elnet_data <- function(data, gene_column){
  tryCatch({
    # Extract treatment groupings.
    data.long <- data %>%
      tidyr::pivot_longer(- !!dplyr::sym(gene_column),
                          names_to = "sample", values_to = "abundance") %>%
      dplyr::rowwise() %>%
      dplyr::mutate(treat = tail(unlist(strsplit(sample,"_")),1)) %>%
      dplyr::ungroup()

    # Genes as columns, samples as rows.
    data.wide <- data.long %>%
      tidyr::pivot_wider(names_from = !!dplyr::sym(gene_column),
                         values_from = abundance)
  },
  error = function(e) stop(paste("unable to prep data for Elastic-net logistic regression:",e))
  )
  return(data.wide)
}
