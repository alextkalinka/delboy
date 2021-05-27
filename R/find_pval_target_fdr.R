# Helper functions.
.smooth_fdr_pvalue <- function(data, stat, dtype){
  col <- rlang::sym(stat)
  # Drop vals down to 0,0.
  zeros <- as.data.frame(smoothr::smooth_ksmooth(as.matrix(rbind(tidyr::tibble(!!col := 0, pvalue = 0),
                                                                         data %>% 
                                                                           dplyr::filter(!!col == min(!!col, na.rm = T)))))) %>%
    dplyr::rename(!!col := V1, pvalue = V2)
  data <- rbind(data, zeros) %>%
    dplyr::arrange(pvalue)
  # Rolling means.
  data <- tidyr::tibble(!!col := zoo::rollmean(data[,stat], 5),
                        pvalue = zoo::rollmean(data$pvalue, 5))
  # Run twice to help smooth.
  data <- tidyr::tibble(!!col := zoo::rollmean(data[,stat], 5),
                        pvalue = zoo::rollmean(data$pvalue, 5))
  # Add zero point back (can be lost when averaging).
  data <- rbind(data, tidyr::tibble(!!col := 0, pvalue = 0)) %>%
    dplyr::arrange(pvalue)
  
  # Kernel smooth.
  data <- as.data.frame(smoothr::smooth_ksmooth(as.matrix(data))) %>%
    dplyr::mutate(type = dtype)
  colnames(data) <- c(stat,"pvalue","type")
  return(data)
}


#' find_pval_target_fdr
#' 
#' Returns the unadjusted p-value that corresponds to a target FDR.
#' 
#' @param data A data frame containing FDR estimates and the corresponding unadjusted p-values.
#' @param target_FDR A numerical value (0-100) giving the target FDR percent (%).
#' @return A list containing the following elements:
#' * `pvalue_target_FDR`: The p-value corresponding to the target FDR.
#' * `abs_Log2FoldChange.argmax_KS_dist`: The absolute log fold change FP threshold corresponding to the target FDR.
#' * `data_FDR`: A modified data frame containing smoothed FDR estimates.
#' * `data_Sensitivity`: A modified data frame containing smoothed Sensitivity estimates.
#' * `target_FDR`: A numerical value (0-100) giving the target FDR %.
#' @md
#' @export
#' @importFrom dplyr %>% filter select mutate arrange rename
#' @importFrom zoo rollmean
#' @importFrom smoothr smooth_ksmooth
#' @importFrom utils tail
#' @importFrom tidyr tibble
#' @importFrom rlang !! sym :=
find_pval_target_fdr <- function(data, target_FDR){
  tryCatch({
    data %<>%
      dplyr::filter(!is.na(pvalue))
    # 1. FDR.
    # 1A. Excluding predicted FPs.
    data_fdr_excl <- data %>%
      dplyr::filter(type == "Excl_Pred_FP") %>%
      dplyr::select(FDR.percent, pvalue)

    data_fdr_excl <- .smooth_fdr_pvalue(data_fdr_excl, "FDR.percent", "Excl_Pred_FP")
    
    # 1B. All.
    data_fdr_all <- data %>%
      dplyr::filter(type == "All") %>%
      dplyr::select(FDR.percent, pvalue)
    
    data_fdr_all <- .smooth_fdr_pvalue(data_fdr_all, "FDR.percent", "All")
    
    # 2. Sensitivity.
    # 2A. Excluding predicted FPs.
    data_sens_excl <- data %>%
      dplyr::filter(type == "Excl_Pred_FP") %>%
      dplyr::select(Sensitivity.percent, pvalue)
    
    data_sens_excl <- .smooth_fdr_pvalue(data_sens_excl, "Sensitivity.percent", "Excl_Pred_FP")
    
    # 2B. All.
    data_sens_all <- data %>%
      dplyr::filter(type == "All") %>%
      dplyr::select(Sensitivity.percent, pvalue)
    
    data_sens_all <- .smooth_fdr_pvalue(data_sens_all, "Sensitivity.percent", "All")
    
    # 3. Unadjusted p-value that corresponds to target FDR.
    targfdr <- data_fdr_excl %>%
      dplyr::filter(type == "Excl_Pred_FP") %>%
      dplyr::filter(abs(FDR.percent - target_FDR) == min(abs(FDR.percent - target_FDR), na.rm = T))
    pval_targfdr <- utils::tail(targfdr$pvalue,1)
    lfc_thresh_targfdr <- utils::tail((data %>%
                                        dplyr::filter(abs(pvalue - pval_targfdr) == min(abs(pvalue - pval_targfdr),na.rm=T)))$abs_Log2FoldChange.argmax_KS_dist,1)
  
    return(list(pvalue_target_FDR = pval_targfdr,
                abs_Log2FoldChange.argmax_KS_dist = lfc_thresh_targfdr,
                data_FDR = rbind(data_fdr_all, data_fdr_excl),
                data_Sensitivity = rbind(data_sens_all, data_sens_excl),
                target_FDR = target_FDR))
  },
  error = function(e) stop(paste("unable to find p-value corresponding to a target FDR:",e))
  )
}
