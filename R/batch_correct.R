#' batch_correct
#'
#' Performs a batch correction on normalized count data using `sva`'s `ComBat` method.
#'
#' @param data A data frame of normalized count data.
#' @param group_1 A character vector naming the columns that belong to group 1.
#' @param group_2 A character vector naming the columns that belong to group 2.
#' @param gene_column A character string naming the column containing gene names.
#' @param method A character string indicating the `ComBat` method. Can be either `np` (non-parametric), or `p` (parametric). Defaults to `np`.
#'
#' @return A data frame of batch-corrected data (can contain negative values).
#' @importFrom sva ComBat
#' @importFrom dplyr mutate select everything
#' @importFrom rlang := !! sym
#' @export
batch_correct <- function(data, group_1, group_2, gene_column, method = "np"){
  tryCatch({
    # 1. Convert to matrix for ComBat.
    data.m <- as.matrix(data[,c(group_1, group_2)])
    rownames(data.m) <- data[,gene_column]

    # 2. Define treatment groups.
    tr <- delboy::make_treat_factor(data.m, group_1, group_2)

    # 3. Run ComBat.
    data.bc <- sva::ComBat(data.m, batch = tr, par.prior = FALSE)

    # 4. Convert back to data frame.
    data.bc <- as.data.frame(data.bc) %>%
      dplyr::mutate(!!rlang::sym(gene_column) := rownames(.)) %>%
      dplyr::select(!!rlang::sym(gene_column), dplyr::everything())
  },
  error = function(e) stop(paste("unable to batch correct data:",e))
  )
  return(data.bc)
}
