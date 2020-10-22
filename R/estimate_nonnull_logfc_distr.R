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
#' @param logfc A vector of logFC values derived from `delboy::run_deseq2`.
#'
#' @return An object of class delboy_logfc
#' @export
#' @importFrom locfdr locfdr
#' @importFrom dplyr between
estimate_nonnull_logfc_distr <- function(logfc){
  misfit <- FALSE
  tryCatch({
    withCallingHandlers({
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
