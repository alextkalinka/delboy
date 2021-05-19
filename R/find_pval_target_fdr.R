# Helper functions.
.chull_df <- function(df){
  chull_inds <- grDevices::chull(as.matrix(df))
  df <- df[chull_inds,]
  return(df)
}


#' find_pval_target_fdr
#' 
#' Returns the unadjusted p-value that corresponds to a target FDR.
#' 
#' @param data A data frame containing FDR estimates and the corresponding unadjusted p-values.
#' @param target_FDR A numerical value (0-1) giving the target FDR.
#' @return A list containing the following elements:
#' * `pvalue_target_FDR`: A numerical value giving the p-value corresponding to the target FDR.
#' * `data`: A modified data frame containing LOESS FDR predictions.
#' * `target_FDR`: A numerical value (0-1) giving the target FDR.
#' @md
#' @export
#' @importFrom dplyr %>% filter select mutate
#' @importFrom zoo rollmean
find_pval_target_fdr <- function(data, target_FDR){
  tryCatch({
    data %<>%
      dplyr::filter(!is.na(pvalue))
    # Rolling mean used to smooth FDR curve.
    data_fdr_excl <- data %>%
      dplyr::filter(type == "Excl_Pred_FP") %>%
      dplyr::select(FDR.percent, pvalue)
    # Drop down to 0,0.
    fdr_excl_zero <- as.data.frame(smoothr::smooth_ksmooth(as.matrix(rbind(data.frame(FDR.percent = 0, pvalue = 0),
                                                             data_fdr_excl %>% 
                                                               dplyr::filter(FDR.percent == min(FDR.percent,na.rm = T)))))) %>%
      dplyr::rename(FDR.percent = V1, pvalue = V2)
    data_fdr_excl <- rbind(data_fdr_excl, fdr_excl_zero) %>%
      dplyr::arrange(pvalue)
    data_fdr_excl <- data.frame(FDR.percent = zoo::rollmean(data_fdr_excl$FDR.percent, 5),
                                pvalue = zoo::rollmean(data_fdr_excl$pvalue, 5))
    data_fdr_excl <- as.data.frame(smoothr::smooth_ksmooth(as.matrix(data_fdr_excl))) %>%
      dplyr::mutate(type = "Excl_Pred_FP")
    colnames(data_fdr_excl) <- c("FDR.percent","pvalue","type")
    
    data_fdr_all <- data %>%
      dplyr::filter(type == "All") %>%
      dplyr::select(FDR.percent, pvalue)
    # Drop down to 0,0.
    all_excl_zero <- as.data.frame(smoothr::smooth_ksmooth(as.matrix(rbind(data.frame(FDR.percent = 0, pvalue = 0),
                                                                           data_fdr_all %>% 
                                                                             dplyr::filter(FDR.percent == min(FDR.percent,na.rm = T)))))) %>%
      dplyr::rename(FDR.percent = V1, pvalue = V2)
    data_fdr_all <- rbind(data_fdr_all, all_excl_zero) %>%
      dplyr::arrange(pvalue)
    data_fdr_all <- data.frame(FDR.percent = zoo::rollmean(data_fdr_all$FDR.percent, 5),
                               pvalue = zoo::rollmean(data_fdr_all$pvalue, 5))
    data_fdr_all <- as.data.frame(smoothr::smooth_ksmooth(as.matrix(data_fdr_all))) %>%
      dplyr::mutate(type = "All")
    colnames(data_fdr_all) <- c("FDR.percent","pvalue","type")
    
    return(list(#pvalue_target_FDR = pv_fdr_targ,
                FDR_data = rbind(data_fdr_all,data_fdr_excl),
                target_FDR = target_FDR))
  },
  error = function(e) stop(paste("unable to find p-value corresponding to a target FDR:",e))
  )
}
