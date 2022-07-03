#' normalize_library_depth_relative
#'
#' Normalizes library depth from a counts file by dividing the counts of each sample by its total and multiplying all samples by a constant total number of reads.
#'
#' @param counts A data frame of counts for each sample in the study (samples as columns, gRNAs as rows).
#' @param norm_total An integer specifying the fixed total number of reads each sample will be set to (e.g. 2e7). If `NULL`, then the median total across the samples is used.
#'
#' @return A data frame in which the sample columns are library-depth normalized.
#' @export
normalize_library_depth_relative <- function(counts, norm_total){
  if(ncol(counts) < 3)
    stop("expecting 'counts' to have at least 3 columns")

  tryCatch({
    csums <- colSums(counts[,3:ncol(counts)])
    if(is.null(norm_total)){
      tot <- median(csums, na.rm=T)
    }else{
      tot <- norm_total
    }
    counts[,3:ncol(counts)] <- as.data.frame(tot * sweep(counts[,3:ncol(counts)], 2, csums, "/"))
  },
  error = function(e) stop(paste("unable to normalize (relative) library depth:",e))
  )
  return(counts)
}
