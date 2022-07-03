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
#' @param signed Logical indicating whether signed logFC distributions should be returned (defaults to `FALSE`).
#'
#' @return An object of class delboy_logfc
#' @export
#' @importFrom locfdr locfdr
#' @importFrom dplyr between
estimate_nonnull_logfc_distr <- function(logfc, signed = FALSE){
  tryCatch({
    logfc <- .cleanup_qnorm(logfc)
    if(!signed){
      # To stabilise estimates we exclude extreme logFC asymmetries between negative and positive ends of the distribution.
      minr <- .min_range(logfc)
      lfc <- logfc[dplyr::between(logfc,-minr, minr)]
    }else{
      lfc <- logfc
    }
    tryCatch(
      lfdr.lfc <- locfdr::locfdr(lfc, plot = 0),
      error = function(e) cat("locfdr cannot estimate non-null logFC distribution\n")
    )
    if(!exists("lfdr.lfc"))
      return(list(fit_ok = FALSE))
    if(!signed){
      non_null.dens <- lfdr.lfc$mat[1:which(lfdr.lfc$mat[,11]==0)[1],11]
      non_null.lfc <- lfdr.lfc$mat[1:which(lfdr.lfc$mat[,11]==0)[1],1]
      non_null.pos.dens <- NA
      non_null.pos.lfc <- NA
      non_null.neg.dens <- NA
      non_null.neg.lfc <- NA
    }else{
      non_null.dens <- NA
      non_null.lfc <- NA
      non_null.neg.dens <- lfdr.lfc$mat[1:which(lfdr.lfc$mat[,11]==0)[1],11]
      non_null.neg.lfc <- lfdr.lfc$mat[1:which(lfdr.lfc$mat[,11]==0)[1],1]
      lrev <- lfdr.lfc$mat
      lrev[,1] <- rev(lrev[,1])
      lrev[,11] <- rev(lrev[,11])
      non_null.pos.dens <- lrev[1:which(lrev[,11]==0)[1],11]
      non_null.pos.lfc <- lrev[1:which(lrev[,11]==0)[1],1]
    }

    ret <- list(non_null.dens = non_null.dens,
                non_null.lfc = non_null.lfc,
                non_null.pos.dens = non_null.pos.dens,
                non_null.pos.lfc = non_null.pos.lfc,
                non_null.neg.dens = non_null.neg.dens,
                non_null.neg.lfc = non_null.neg.lfc,
                fit_ok = TRUE)
  },
  error = function(e) stop(paste("unable to estimate the non-null logFC distribution:",e))
  )
  return(ret)
}
