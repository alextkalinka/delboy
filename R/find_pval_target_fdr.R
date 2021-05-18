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
#' @importFrom grDevices chull
find_pval_target_fdr <- function(data, target_FDR){
  tryCatch({
    # Convex hull of points -> densified.
    data_fdr_excl <- data %>%
      dplyr::filter(type == "Excl_Pred_FP") %>%
      dplyr::select(FDR.percent, pvalue) %>%
      .chull_df()
    data_fdr_excl <- as.data.frame(smoothr::smooth_ksmooth(as.matrix(data_fdr_excl))) %>%
      dplyr::mutate(type = "Excl_Pred_FP")
    colnames(data_fdr_excl) <- c("FDR.percent","pvalue","type")
    
    data_fdr_all <- data %>%
      dplyr::filter(type == "All") %>%
      dplyr::select(FDR.percent, pvalue) %>%
      .chull_df()
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
