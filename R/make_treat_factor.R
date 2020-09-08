#' make_treat_factor
#'
#' Returns a treatment factor for use in batch-correction and differential-expression analyses.
#'
#' @param data A matrix of counts with named columns.
#' @param group_1 A character string naming the columns that belong to group 1.
#' @param group_2 A character string naming the columns that belong to group 2.
#'
#' @return A factor vector of treatments corresponding to the columns of the input `data`.
#' @export
make_treat_factor <- function(data, group_1, group_2){
  tryCatch({
    tr <- rep("group_1",ncol(data))
    tr[which(colnames(data) %in% group_2)] <- "group_2"
    tr <- factor(tr)
  },
  error = function(e) stop(paste("unable to make treatment factors:",e))
  )
  return(tr)
}
