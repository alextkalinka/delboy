# Helper functions.

# locfdr call.
.call_locfdr <- function(qvals){
  misfit <- FALSE
  withCallingHandlers({
    lfdr <- locfdr::locfdr(qvals, plot = 0)
  },
  warning = function(w) misfit <<- gsub("^.*?misfit = (\\S+)\\.\\s+?.*$","\\1",w)
  )
  mf <- .calculate_locfdr_misfit(qvals, lfdr)
  return(list(locfdr = lfdr, misfit = mf, misfit_reported = misfit))
}


# Calculate locfdr misfit (code adapted from 'locfdr' source code).
.calculate_locfdr_misfit <- function(zz, lfdr){
  f <- lfdr$mat[,"f"]
  # Calculate misfit.
  df <- 7
  bre <- 120
  lo <- min(zz)
  up <- max(zz)
  zzz <- pmax(pmin(zz, up), lo)
  breaks <- seq(lo, up, length = bre)
  zh <- hist(zzz, breaks = breaks, plot = F)
  y <- zh$counts
  K <- length(y)
  D <- (y - f)/(f + 1)^0.5
  D <- sum(D[2:(K - 1)]^2)/(K - 2 - df)
  return(D)
}


# Smooth p-values from DESeq2 to enable 'locfdr' estimate of non-null abundance.
# DESeq2 sometimes produces an excess of p-values > 0.9 and this lumpiness in the p-value distribution can hamper the Gaussian mixture models used by locfdr.
.smooth_pvals <- function(pval, start=0.5){
  pick <- dplyr::between(pval,start,1)
  tot <- sum(pick)
  # 28 breaks with 7576 pvals.
  num_breaks <- round(0.0037*tot)
  if(num_breaks <= 1)
    num_breaks <- 7
  pv <- data.frame(pvalue = pval[pick]) %>%
    dplyr::mutate(pv_cuts = cut(pvalue, breaks = num_breaks)) %>%
    dplyr::group_by(pv_cuts) %>%
    dplyr::mutate(midp = stats::median(pvalue,na.rm=T)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(prob = midp/sum(unique(midp))) %>%
    dplyr::group_by(pv_cuts) %>%
    dplyr::sample_n(round(tot * unique(prob)), replace=T) %>%
    dplyr::ungroup()

  ret <- c(pval[!pick], pv$pvalue)
  diff <- length(pval) - length(ret)
  if(diff > 0){
    ret <- c(ret,rep(ret[1],diff))
  }else if(diff < 0){
    ret <- ret[1:length(pval)]
  }
  return(ret)
}


# Lots of small p-values can lead to an asymmetric Gaussian with a long left tail when transformed.
.symmetrise_q_vals <- function(qval){
  qval <- qval[!qval > 0]
  qval <- append(qval, abs(qval))
  return(qval)
}


.cleanup_qnorm <- function(qn){
  qn <- qn[!is.infinite(qn)]
  qn <- qn[!is.na(qn)]
  return(qn)
}

#' estimate_number_non_nulls
#'
#' Estimates the number of non-null cases in a dataset using `locfdr` with raw p-values as inputs.
#'
#' @param pvals A numeric vector of raw (unadjusted) p-values.
#'
#' @return A list containing an integer giving an estimate of the number of non-null cases, and an estimate of the misfit of the `locfdr` mixture model.
#' @importFrom locfdr locfdr
#' @importFrom dplyr between group_by ungroup mutate sample_n %>%
#' @importFrom stats median qnorm
#' @export
estimate_number_non_nulls <- function(pvals){
  pvals <- pvals[!is.na(pvals)]
  qvals <- list()
  fits <- list()
  tryCatch({
    # 1. Unaltered p-vals.
    qvals[[1]] <- .cleanup_qnorm(stats::qnorm(pvals))
    fits[[1]] <- .call_locfdr(qvals[[1]])
    
    # 2. Symmetrise transformed p-vals.
    qvals[[2]] <- .symmetrise_q_vals(.cleanup_qnorm(stats::qnorm(pvals)))
    fits[[2]] <- .call_locfdr(qvals[[2]])
    
    # 3. Smooth p-vals (start=0.5).
    qvals[[3]] <- .cleanup_qnorm(stats::qnorm(.smooth_pvals(pvals, start = 0.5)))
    fits[[3]] <- .call_locfdr(qvals[[3]])
    
    # 4. Smooth p-vals (start=0.9).
    qvals[[4]] <- .cleanup_qnorm(stats::qnorm(.smooth_pvals(pvals, start = 0.9)))
    fits[[4]] <- .call_locfdr(qvals[[4]])
    
    misfits <- unlist(lapply(fits, function(x) x$misfit))
    min_misfit <- which(misfits == min(misfits, na.rm=T))[1]
    best_fit <- fits[[min_misfit]]
    qn <- qvals[[min_misfit]]
    lfdr <- best_fit$locfdr

    num.non_null <- round(sum(lfdr$mat[1:which(lfdr$mat[,11]==0)[1], 11]))
  },
  error = function(e) stop(paste("unable to estimate number of non-nulls:",e))
  )
  return(list(num.non_null = num.non_null, original.num.non_null = num.non_null,
              misfit = best_fit$misfit, qvals = qn, locfdr = lfdr, fit_type = min_misfit,
              all_misfits = misfits))
}
