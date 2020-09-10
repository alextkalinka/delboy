#' prep_elnet_data
#'
#' Prepare data for use in an Elastic-net logistic regression.
#'
#' @param data A data frame produced by `delboy::prep_bthin_matrix_diffrep`.
#' @param group_1 A character vector naming the columns that belong to group 1.
#' @param group_2 A character vector naming the columns that belong to group 2.
#' @param gene_column A character string naming the column containing gene names.
#'
#' @return A data frame.
#' @export
#' @importFrom dplyr rowwise mutate ungroup sym %>%
#' @importFrom tidyr pivot_longer pivot_wider
prep_elnet_data <- function(data, group_1, group_2, gene_column){
  tryCatch({
    # Extract treatment groupings.
    data.long <- data %>%
      tidyr::pivot_longer(- !!dplyr::sym(gene_column),
                          names_to = "sample", values_to = "abundance") %>%
      dplyr::rowwise() %>%
      dplyr::mutate(treat = ifelse(sample %in% group_1,"group_1","group_2")) %>%
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
