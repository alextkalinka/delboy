#' prep_count_matrix
#'
#' Prepares an integer count matrix from a data frame of normalized counts so the data may be used by functions requiring matrix inputs.
#'
#' @param data A data frame of normalized count data (can contain negative values).
#' @param group_1 A character vector naming the columns that belong to group 1.
#' @param group_2 A character vector naming the columns that belong to group 2.
#' @param gene_column A character string naming the column containing gene names.
#'
#' @return A matrix of integer counts with any negative values set to zero - rownames equal gene names.
#' @export
prep_count_matrix <- function(data, group_1, group_2, gene_column){
  tryCatch({
    as.matrix(data[,c(group_1, group_2)])
    rownames(data.m) <- data[,gene_column]
    data.m <- apply(data.m, 2, round)
    data.m <- apply(data.m, 2, function(x) ifelse(x < 0, 0, x))
  },
  error = function(e) stop(paste("unable to prep count matrix:",e))
  )
  return(data.m)
}
