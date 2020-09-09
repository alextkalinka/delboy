#' run_elnet_logist_reg
#'
#' Run an Elastic-net logistic regression using `glmnet`.
#'
#' @param data A matrix of normalized, batch-corrected counts - samples as rows, genes as columns.
#' @param treat A factor vector grouping samples into one of at most two groups.
#' @param alpha The elastic-net resgression penalty, between 0 and 1.
#'
#' @return An object of class `delboy_elnet.
#' @export
#' @importFrom glmnet glmnet cv.glmnet
run_elnet_logist_reg <- function(data, treat, alpha){
  tryCatch({
    withCallingHandlers({
      fit.elnet <- glmnet::glmnet(data, factor(treat), family = "binomial", alpha = alpha)
      fit.cv_dev <- glmnet::cv.glmnet(data, factor(treat), family = "binomial", alpha = alpha,
                                  type.measure = "deviance")
      fit.cv_class <- glmnet::cv.glmnet(data, factor(treat), family = "binomial", alpha = alpha,
                                      type.measure = "class")
    },
    warning = function(w)


  },
  error = function(e) stop(paste("unable to run Elastic-net logistic regression:",e))
  )

}
