#' make_design_matrix
#'
#' Makes a design matrix for use by `seqgendiff`.
#'
#' @param group_1 A character vector naming the columns that belong to group 1.
#' @param group_2 A character vector naming the columns that belong to group 2.
#' @param random_treat Logical indicating whether the number of treatment samples should be a random number between 2 and `length(group_1) + length(group_2) - 2`, or not. Defaults to `FALSE` - there will be `length(group_2)` treatment samples.
#'
#' @return A design matrix.
#' @export
make_design_matrix <- function(group_1, group_2, random_treat = FALSE){
  tryCatch({
    num_samps <- length(group_1) + length(group_2)
    treat_cols <- rep(0,num_samps)
    if(!random_treat){
      # Randomly sample the treatment columns.
      treat_cols[sample(1:num_samps, length(group_2), replace = F)] <- 1
    }else{
      num_treat <- sample(2:(num_samps-2),1)
      treat_cols[sample(1:length(treat_cols), num_treat, replace=F)] <- 1
    }
    design_mat <- matrix(treat_cols)
    colnames(design_mat) <- "treatment"
  },
  error = function(e) stop(paste("unable to make design matrix:",e))
  )
  return(design_mat)
}
