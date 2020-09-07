# Smooth p-values > 0.5 from DESeq2 to enable 'locfdr' estimate of non-null abundance.
# DESeq2 sometimes produces an excess of p-values > 0.9 and this can hamper the Gaussian mixture models used by locfdr.
.smooth_pvals <- function(pval, start=0.5){
  pick <- between(pval,start,1)
  tot <- sum(pick)
  pv <- data.frame(pvalue = pval[pick]) %>%
    mutate(pv_cuts = cut(pvalue, breaks=28)) %>%
    group_by(pv_cuts) %>%
    mutate(midp = median(pvalue,na.rm=T)) %>%
    ungroup %>%
    mutate(prob = midp/sum(unique(midp))) %>%
    group_by(pv_cuts) %>%
    sample_n(round(tot * unique(prob)), replace=T) %>%
    ungroup

  ret <- c(pval[!pick], pv$pvalue)
  diff <- length(pval) - length(ret)
  if(diff > 0){
    ret <- c(ret,rep(ret[1],diff))
  }else if(diff < 0){
    ret <- ret[1:length(pval)]
  }
  return(ret)
}


#' estimate_number_non_nulls
#'
#' Estimates the number of non-null cases in a dataset using `locfdr` with raw p-values as inputs.
#'
#' @param pvals A numeric vector of raw (unadjusted) p-values.
#'
#' @return A list containing an integer giving an estimate of the number of non-null cases, and an estimate of the misfit of the mixture model.
#' @importFrom locfdr locfdr
#' @export
estimate_number_non_nulls <- function(pvals){
  misfit <- FALSE
  tryCatch({
    withCallingHandlers({
      qn <- qnorm(.smooth_pvals(pvals))
      lfdr <- locfdr::locfdr(qn)
    },
    warning = function(w) misfit <<- gsub("^.*?misfit = (\\d+?).*$","\\1",w)
    )
    num.non_null <- round(sum(lfdr$mat[1:which(lfdr$mat[,11]==0)[1], 11]))
  },
  error = function(e) stop(paste("unable to estimate number of non-nulls:",e))
  )
  return(list(num.non_null = num.non_null, misfit = misfit))
}
