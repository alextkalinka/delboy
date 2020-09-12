# Helper functions.
.plotFCExpr <- function(delboy){
  grid <- delboy$performance_eval$svm_validation$grid
  pl <- delboy$hits_original_validation %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(data = delboy$hits_original_validation %>%
                 filter(hit_type!="False_Negative"),
               ggplot2::aes(log10_baseExpr, abs_log2FoldChange, color=hit_type)) +
    ggplot2::geom_contour(data = grid,
                          ggplot2::aes(log10_baseExpr, abs_log2FoldChange, z=Predicted_FP)) +
    ggplot2::facet_grid(~data_type) +
    ggplot2::xlim(0.5,4) +
    ggplot2::ylim(0,1.5)
  print(pl)
}


# Helper functions.
.plotFCExprFN <- function(delboy){
  data <- delboy$hits_original_validation %>%
    dplyr::filter(hit_type != "Positive")
  grid <- delboy$performance_eval$svm_validation$grid
  pl <- data %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(data = data,
                        ggplot2::aes(log10_baseExpr, abs_log2FoldChange, color=hit_type)) +
    ggplot2::geom_contour(data = grid,
                          ggplot2::aes(log10_baseExpr, abs_log2FoldChange, z=Predicted_FP)) +
    ggplot2::xlim(0.5,4) +
    ggplot2::ylim(0,1.5)
  print(pl)
}


#' plot.delboy
#'
#' Plotting for `delboy` objects.
#'
#' @param delboy_res Output from `delboy::run_delboy`.
#' @param type A character string naming the plot type: `fc_expr`.
#' @param ... Other arguments to be passed to `plot`.
#'
#' @return Used for side-effect of plotting.
#' @export
#' @importFrom ggplot2 ggplot aes facet_grid geom_point geom_contour
plot.delboy <- function(delboy_res, type, ...){
  switch(type,
         fc_expr = .plotFCExpr(delboy_res),
         fc_expr_FN = .plotFCExprFN(delboy_res))
}
