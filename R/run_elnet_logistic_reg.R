#' run_elnet_logist_reg
#'
#' Run an Elastic-net logistic regression using `glmnet`.
#'
#' @param data A matrix of normalized, batch-corrected counts - samples as rows, genes as columns.
#' @param treat A factor vector grouping samples into one of at most two groups.
#' @param alpha The elastic-net resgression penalty, between 0 and 1.
#'
#' @return An object of class `delboy_elnet`.
#' @export
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats coef
#' @importFrom utils tail
run_elnet_logistic_reg <- function(data, treat, alpha){
  tryCatch({
    max.num_reps <- sort(table(treat),decreasing = T)[1]
    warn <- NULL
    withCallingHandlers({
      fit.elnet <- glmnet::glmnet(data, factor(treat), family = "binomial", alpha = alpha)
      if(max.num_reps > 3){
        fit.cv_dev <- glmnet::cv.glmnet(data, factor(treat), family = "binomial", alpha = alpha,
                                    type.measure = "deviance")
        fit.cv_class <- glmnet::cv.glmnet(data, factor(treat), family = "binomial", alpha = alpha,
                                        type.measure = "class")
      }else{
        fit.cv_dev <- NA
        fit.cv_class <- NA
      }
    },
    warning = function(w) warn <<- append(warn,w)
    )
    # Extract non-zero coefficients at point where fit is best (lambda min).
    if(max.num_reps > 3){
      genes.elnet <- stats::coef(fit.elnet, s = fit.cv_dev$lambda.min)
    }else{
      genes.elnet <- stats::coef(fit.elnet, s = utils::tail(fit.elnet$lambda,1))
    }

    # Over-represented in group 2.
    genes.up <- genes.elnet[genes.elnet[,1] > 0,]
    genes.up <- names(genes.up)[2:length(genes.up)]

    # Under-represented in group 2.
    genes.down <- genes.elnet[genes.elnet[,1] < 0,]
    genes.down <- names(genes.down)[2:length(genes.down)]

    # Build object of 'delboy_elnet' class.
    ret <- list(fit = fit.elnet,
                cvfit = fit.cv_dev,
                cvfit.class = fit.cv_class,
                genes.elnet = genes.elnet,
                genes.up = genes.up,
                genes.down = genes.down,
                warnings = warn)
    class(ret) <- "delboy_elnet"
  },
  error = function(e) stop(paste("unable to run Elastic-net logistic regression:",e))
  )
  return(ret)
}
