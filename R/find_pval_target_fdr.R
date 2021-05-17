#' find_pval_target_fdr
#' 
#' Returns the unadjusted p-value that corresponds to a target FDR using a LOESS fit to the FDR estimates.
#' 
#' @param data A data frame containing FDR estimates and the corresponding unadjusted p-values.
#' @param target_FDR A numerical value (0-1) giving the target FDR.
#' @return A list containing the following elements:
#' * `pvalue_target_FDR`: A numerical value giving the p-value corresponding to the target FDR.
#' * `data`: A modified data frame containing LOESS FDR predictions.
#' * `loess_fit_fdr_all` An object of class `loess`.
#' * `loess_fit_fdr_excl` An object of class `loess`.
#' * `loess_fit_sens_all` An object of class `loess`.
#' * `loess_fit_sens_excl` An object of class `loess`.
#' * `target_FDR`: A numerical value (0-1) giving the target FDR.
#' @md
#' @export
#' @importFrom dplyr %>% filter
#' @importFrom stats loess predict
find_pval_target_fdr <- function(data, target_FDR){
  tryCatch({
    # Densify data.
    data_excl <- as.data.frame(smoothr::smooth_ksmooth(as.matrix(data %>%
                                                              dplyr::filter(type == "Excl_Pred_FP") %>%
                                                              dplyr::select(pvalue,FDR.percent)))) %>%
      dplyr::mutate(type = "Excl_Pred_FP")
    colnames(data_excl) <- c("pvalue","FDR.percent","type")
    data_all <- as.data.frame(smoothr::smooth_ksmooth(as.matrix(data %>%
                                                                   dplyr::filter(type == "All") %>%
                                                                   dplyr::select(pvalue,FDR.percent)))) %>%
      dplyr::mutate(type = "All")
    colnames(data_all) <- c("pvalue","FDR.percent","type")
    data <- rbind(data_all,data_excl)
    return(data)
    # LOESS fits.
    loess_fdr_all <- stats::loess(pvalue ~ FDR.percent, data %>%
                                    dplyr::filter(type == "All"),
                                  control = loess.control(surface = "direct"))
    #loess_sens_all <- stats::loess(Sensitivity.percent ~ pvalue, data %>%
    #                                  dplyr::filter(type == "All"),
    #                               control = loess.control(surface = "direct"))
    loess_fdr_excl <- stats::loess(pvalue ~ FDR.percent, data %>%
                                    dplyr::filter(type == "Excl_Pred_FP"),
                                   control = loess.control(surface = "direct"))
    #loess_sens_excl <- stats::loess(pvalue ~ Sensitivity.percent, data %>%
    #                                 dplyr::filter(type == "Excl_Pred_FP"),
    #                                control = loess.control(surface = "direct"))
    # pvalue corresponding to target FDR after exclusion of FPs.
    pv_fdr_targ <- stats::predict(loess_fdr_excl, newdata = data.frame(FDR.percent = target_FDR))
    fdr_excl <- data.frame(FDR.percent = seq(0,100,by=1), type = "Excl_Pred_FP") %>%
      dplyr::mutate(pvalue = stats::predict(loess_fdr_excl, newdata = .))
    fdr_all <- data.frame(FDR.percent = seq(0,100,by=1), type = "All") %>%
      dplyr::mutate(pvalue = stats::predict(loess_fdr_all, newdata = .))
    return(list(pvalue_target_FDR = pv_fdr_targ,
                FDR_data = rbind(fdr_all,fdr_excl),
                target_FDR = target_FDR))
  },
  error = function(e) stop(paste("unable to find p-value corresponding to a target FDR:",e))
  )
}
