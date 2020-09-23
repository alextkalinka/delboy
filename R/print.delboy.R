#' print.delboy
#'
#' Prints summary of a `delboy` object to the console.
#'
#' @param x Output from `delboy::run_delboy`.
#' @param ... Other arguments to be passed to `print`.
#'
#' @return Prints to the console.
#' @export
print.delboy <- function(x, ...){
  cat(paste("*** delboy summary ***\n\nPerformance estimates:\n",
            "Sensitivity (%):\n  delboy: ",round(x$performance_stats_corr_FP$Sensitivity.percent[1],3),
            " (",x$performance_stats_corr_FP$Num_true_calls[1]," genes)",
            "\n  DESeq2: ",round(x$performance_stats_corr_FP$Sensitivity.percent[2],3),
            " (",x$performance_stats_corr_FP$Num_true_calls[2]," genes)",
            "\nFDR (%):\n  delboy: ",round(x$performance_stats_corr_FP$FDR.percent[1],3),
            " (",x$performance_stats_corr_FP$Num_false_calls[1]," genes)",
            "\n  DESeq2: ",round(x$performance_stats_corr_FP$FDR.percent[2],3),
            " (",x$performance_stats_corr_FP$Num_false_calls[2]," genes)",
            sep=""))
}
