# Helper functions.
.plotFCExpr <- function(delboy, xlim, ylim){
  if(is.null(xlim)){
    xlim_l <- 0.5
    xlim_u <- 4
  }else{
    xlim_l <- xlim[1]
    xlim_u <- xlim[2]
  }
  if(is.null(ylim)){
    ylim_l <- 0
    ylim_u <- 1.5
  }else{
    ylim_l <- ylim[1]
    ylim_u <- ylim[2]
  }
  db <- delboy$performance_eval$svm_validation$decision_boundary
  pl <- delboy$hits_original_validation %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(data = delboy$hits_original_validation %>%
                          filter(hit_type!="False_Negative"),
                 ggplot2::aes(log10_baseExpr, abs_log2FoldChange, color=hit_type)) +
    ggplot2::geom_line(data = db,
                       ggplot2::aes(log10_baseExpr, abs_log2FoldChange), size = 0.8) +
    ggplot2::facet_grid(~data_type) +
    ggplot2::coord_cartesian(xlim = c(xlim_l, xlim_u),
                             ylim = c(ylim_l, ylim_u))
  print(pl)
}


.plotFCExprFN <- function(delboy, xlim, ylim){
  if(is.null(xlim)){
    xlim_l <- 0.5
    xlim_u <- 4
  }else{
    xlim_l <- xlim[1]
    xlim_u <- xlim[2]
  }
  if(is.null(ylim)){
    ylim_l <- 0
    ylim_u <- 1.5
  }else{
    ylim_l <- ylim[1]
    ylim_u <- ylim[2]
  }
  data <- delboy$hits_original_validation %>%
    dplyr::filter(hit_type != "Positive")
  db <- delboy$performance_eval$svm_validation$decision_boundary
  pl <- data %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(data = data,
                        ggplot2::aes(log10_baseExpr, abs_log2FoldChange, color=hit_type)) +
    ggplot2::geom_line(data = db,
                       ggplot2::aes(log10_baseExpr, abs_log2FoldChange), size = 0.8) +
    ggplot2::coord_cartesian(xlim = c(xlim_l, xlim_u),
                             ylim = c(ylim_l, ylim_u))
  print(pl)
}


.plotBinDev <- function(delboy){
  if(!is.na(delboy$elnet_results$cvfit)){
    plot(delboy$elnet_results$cvfit)
  }else{
    cat("too few replicates to perform cross-validation")
  }
}


.plotMisClass <- function(delboy){
  if(!is.na(delboy$elnet_results$cvfit.class)){
    plot(delboy$elnet_results$cvfit.class)
  }else{
    cat("too few replicates to perform cross-validation")
  }
}


.plotLFCNonNull <- function(delboy){
  data <- data.frame(log2FoldChange = abs(delboy$non_null$nonnull_lfc$non_null.lfc),
                     Density = delboy$non_null$nonnull_lfc$non_null.dens/length(delboy$non_null$nonnull_lfc$non_null.lfc))
  ggplot2::ggplot(data, ggplot2::aes(log2FoldChange, Density)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = 0, linetype="dashed", color="red") +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::ggtitle("Estimated abs(logFC) distribution for non-null cases")
}


#' plot.delboy
#'
#' Plotting for `delboy` objects.
#'
#' @param x Output from `delboy::run_delboy`.
#' @param type A character string naming the plot type: `fc_expr`, `fc_expr_FN`, `lfc_nonnull`, `deviance`, or `misclass`. Defaults to `fc_expr`.
#' @param xlim xlim values for x-axis. Defaults to `NULL` for `c(0.5,4)`.
#' @param ylim xlim values for y-axis. Defaults to `NULL` for `c(0,1.5)`.
#' @param ... Other arguments to be passed to `plot`.
#'
#' @return Used for side-effect of plotting.
#' @export
#' @importFrom ggplot2 ggplot aes facet_grid geom_point geom_line geom_vline geom_hline ggtitle coord_cartesian
plot.delboy <- function(x, type = "fc_expr", xlim = NULL, ylim = NULL, ...){
  if(!inherits(x,"delboy")) stop(paste("expecting an object of class 'delboy', got:",class(x)))
  switch(type,
         fc_expr = .plotFCExpr(x, xlim, ylim),
         fc_expr_FN = .plotFCExprFN(x, xlim, ylim),
         lfc_nonnull = .plotLFCNonNull(x),
         deviance = .plotBinDev(x),
         misclass = .plotMisClass(x)
         )
}
