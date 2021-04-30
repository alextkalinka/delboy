# Helper functions.
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
  misfit <- FALSE
  tryCatch({
    # 1. Unaltered p-vals.
    misfit_1 <- FALSE
    withCallingHandlers({
      qn <- .cleanup_qnorm(stats::qnorm(pvals))
      lfdr_1 <- locfdr::locfdr(qn, plot = 0)
    },
    warning = function(w) misfit_1 <<- gsub("^.*?misfit = (\\S+)\\.\\s+?.*$","\\1",w)
    )
    
    if(is.character(misfit_1) & !grepl(" ",misfit_1)){
      # 2. Symmetrise transformed p-vals.
      misfit_2 <- FALSE
      withCallingHandlers({
        qn <- .symmetrise_q_vals(.cleanup_qnorm(stats::qnorm(pvals)))
        lfdr_2 <- locfdr::locfdr(qn, plot = 0)
      },
      warning = function(w) misfit_2 <<- gsub("^.*?misfit = (\\S+)\\.\\s+?.*$","\\1",w)
      )

      if(is.character(misfit_2) & !grepl(" ",misfit_2)){
        # 3. Smoothes p-vals (start=0.5).
        misfit_3 <- FALSE
        withCallingHandlers({
          qn <- .cleanup_qnorm(stats::qnorm(.smooth_pvals(pvals)))
          lfdr_3 <- locfdr::locfdr(qn, plot = 0)
        },
        warning = function(w) misfit_3 <<- gsub("^.*?misfit = (\\S+)\\.\\s+?.*$","\\1",w)
        )

        if(is.character(misfit_3) & !grepl(" ",misfit_3)){
          # 4. Smoothes p-vals (start=0.9).
          misfit_4 <- FALSE
          withCallingHandlers({
            qn <- .cleanup_qnorm(stats::qnorm(.smooth_pvals(pvals, start=0.9)))
            lfdr_4 <- locfdr::locfdr(qn, plot = 0)
          },
          warning = function(w) misfit_4 <<- gsub("^.*?misfit = (\\S+)\\.\\s+?.*$","\\1",w)
          )

          if(is.character(misfit_4) & !grepl(" ",misfit_4)){
            if(as.numeric(misfit_3) <= as.numeric(misfit_4)){
              lfdr <- lfdr_3
              misfit <- misfit_3
              fit_type <- 3
            }else{
              lfdr <- lfdr_4
              misfit <- misfit_4
              fit_type <- 4
            }
          }else{
            lfdr <- lfdr_4
            misfit <- misfit_4
            fit_type <- 4
          }
        }else{
          lfdr <- lfdr_3
          misfit <- misfit_3
          fit_type <- 3
        }
      }else{
        lfdr <- lfdr_2
        misfit <- misfit_2
        fit_type <- 2
      }
    }else{
      lfdr <- lfdr_1
      misfit <- misfit_1
      fit_type <- 1
    }

    num.non_null <- round(sum(lfdr$mat[1:which(lfdr$mat[,11]==0)[1], 11]))
  },
  error = function(e) stop(paste("unable to estimate number of non-nulls:",e))
  )
  return(list(num.non_null = num.non_null, original.num.non_null = num.non_null,
              misfit = misfit, qvals = qn, locfdr = lfdr, fit_type = fit_type))
}
