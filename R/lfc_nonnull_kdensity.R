# Helper functions.
.breaks_locfdr <- function(logfc){
  # Based on code in 'locfdr::locfdr'.
  lo <- min(logfc)
  up <- max(logfc)
  breaks <- seq(lo, up, length = 120)
  return(breaks)
}


#' lfc_nonnull_kdensity
#' 
#' Returns a kernel density estimate for the logFC of non-null cases.
#' 
#' @param logfc A vector of logFC values.
#' @param lfdr The return object from `locfdr::locfdr`.
#' @return A list.
#' @mad
#' @export
#' @importFrom dplyr between
#' @importFrom kdensity kdensity
lfc_nonnull_kdensity <- function(logfc, lfdr){
  tryCatch({
    breaks <- .breaks_locfdr(logfc)
    if(length(breaks) != (nrow(lfdr$mat) + 1)) stop("breaks and locfdr output mismatch")
    # 1. Sample logFC according to locfdr estimates 1 - fdr.
    lfc_nonnull <- NULL
    for(i in 1:nrow(lfdr$mat)){
      tfc <- abs(logfc[dplyr::between(logfc, breaks[i], breaks[i+1])])
      num_samp <- pmax(round((1 - lfdr$mat[i,2]) * length(tfc)), 0)
      lfc_nonnull <- append(lfc_nonnull, sample(tfc, num_samp, replace = F))
    }
    kdens <- kdensity::kdensity(lfc_nonnull, start = "gamma", kernel = "gaussian")
    return(list(logfc_nonnull = lfc_nonnull,
                lfc_kdens = kdens))
  },
  error = function(e) stop(paste("unable to calculate kernel density estimate for logFC values:",e))
  )
}
