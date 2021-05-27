# Helper functions.
.min_range <- function(data){
  pos <- data[data > 0]
  neg <- abs(data[data < 0])
  return(min(c(max(pos), max(neg))))
}


#' estimate_nonnull_logfc_distr
#'
#' Estimates the full non-null logFC distribution based on `DESeq2` logFC estimates and mixture model fits employed by `locfdr`.
#'
#' @param logfc A vector of logFC values.
#'
#' @return A list with the following elements:
#' * `non_null.dens`: non-null density estimates.
#' * `non_null.lfc`: logFC mid-points corresponding to the densities.
#' * `misfit`: A numerical estimate of the misfit (larger values indicate worse misfits).
#' @md
#' @export
#' @importFrom locfdr locfdr
#' @importFrom dplyr between
estimate_nonnull_logfc_distr <- function(logfc){
  misfit <- FALSE
  tryCatch({
    withCallingHandlers({
      # Strong asymmetry in the tails can lead to very high misfits.
      minr <- .min_range(logfc)
      lfdr.lfc <- locfdr::locfdr(logfc[dplyr::between(logfc,-minr, minr)], plot = 0)
    },
    warning = function(w) misfit <<- gsub("^.*?misfit = (\\S+)\\.\\s+?.*$","\\1",w)
    )
    non_null.dens <- lfdr.lfc$mat[1:which(lfdr.lfc$mat[,11]==0)[1],11]
    non_null.lfc <- lfdr.lfc$mat[1:which(lfdr.lfc$mat[,11]==0)[1],1]

    ret <- list(non_null.dens = non_null.dens,
                non_null.lfc = non_null.lfc,
                misfit = misfit)
  },
  error = function(e) stop(paste("unable to estimate the non-null logFC distribution:",e))
  )
  return(ret)
}
