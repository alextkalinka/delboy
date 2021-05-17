#' make_design_matrix
#'
#' Makes a design matrix for use by `seqgendiff`.
#'
#' @param group_1 A character vector naming the columns that belong to group 1.
#' @param group_2 A character vector naming the columns that belong to group 2.
#' @param treat_samps A character vector indicating which samples are treatment samples. `NULL` (default) indicates they should be randomly samples (`length(group_2)`).
#' @param random_treat Logical indicating whether the number of treatment samples should be a random number between 2 and `length(group_1) + length(group_2) - 2`, or not. Defaults to `FALSE` - there will be `length(group_2)` treatment samples.
#'
#' @return A design matrix.
#' @export
make_design_matrix <- function(group_1, group_2, treat_samps = NULL, random_treat = FALSE){
  if(!is.null(treat_samps)){
    if(any(!treat_samps %in% c(group_1, group_2)))
      stop(paste("some of all of treat samps not found in sample set.\nTreatment samples:",
                 treat_samps,"\nAll samps:",c(group_1,group_2)))
  }
  tryCatch({
    samples <- c(group_1, group_2)
    num_samps <- length(samples)
    treat_cols <- rep(0,num_samps)
    if(!random_treat){
      if(!is.null(treat_samps)){
        # 'run_delboy' organises columns into group_1:group_2 order.
        treat_cols[samples %in% treat_samps] <- 1
      }else{
        # Randomly sample the treatment columns.
        treat_cols[sample(1:num_samps, length(group_2), replace = F)] <- 1
      }
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
