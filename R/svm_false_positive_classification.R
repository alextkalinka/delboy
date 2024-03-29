#' svm_false_positive_classification
#'
#' Utilises an SVM to learn the decision boundary between false and true positives for validation data.
#'
#' @param data A data frame as produced by `delboy::make_delboy_hit_comparison_table`.
#' @param use_fn Logical indicating whether false negatives should be used in conjunction with any false positives for the decision boundary. Defaults to `FALSE`.
#' @param grid_fc_dens An integer giving the density of grid sampling for the decision boundary on the fold change axis. Defaults to 600.
#' @param grid_ex_dens An integer giving the density of grid sampling for the decision boundary on the expression axis. Defaults to 800.
#'
#' @return A list containing the SVM fit, a data frame containing predicted false positive status, and a data frame produced by `expand.grid`.
#' @export
#' @importFrom e1071 svm
#' @importFrom svmpath svmpath
#' @importFrom dplyr %>% filter select mutate
#' @importFrom magrittr %<>%
#' @importFrom stats predict quantile
svm_false_positive_classification <- function(data, use_fn = FALSE, grid_fc_dens = 600, grid_ex_dens = 800){
  tryCatch({
    if(use_fn){
      if(sum(data$hit_type == "False_Negative") < 10)
        stop("* insufficient False Positive and False Negative data to validate algorithms *")
      # FNs up to 75% quantile used as faux FPs.
      fn_75 <- as.numeric(stats::quantile((data %>%
                                 dplyr::filter(hit_type == "False_Negative"))$abs_log2FoldChange,
                               probs = 0.75))
      data %<>%
        dplyr::mutate(hit_type = as.character(hit_type),
                      hit_type = ifelse((hit_type == "False_Negative" & abs_log2FoldChange < fn_75),
                                         "False_Positive", hit_type))
    }

    # Prep data for SVM.
    data_svm <- data %>%
      dplyr::filter(hit_type != "False_Negative" & !is.na(hit_type) &
                      !is.na(abs_log2FoldChange) &
                      !is.na(log10_baseExpr) & !is.infinite(log10_baseExpr) &
                      !is.infinite(abs_log2FoldChange)) %>%
      dplyr::mutate(False_Positive = as.factor(ifelse(hit_type=="True_Positive",0,1))) %>%
      dplyr::select(False_Positive, abs_log2FoldChange, log10_baseExpr)

    # 1. Try 'svmpath::svmpath' first - due to fully optimized regularization, there is minimal over-fitting.
    #    If there is a singularity or other error, then try 'e1071::svm'.
    e1071 <- FALSE
    tryCatch({
      svm_val <- svmpath::svmpath(as.matrix(data_svm %>% select(-False_Positive)),
                                  ifelse(data_svm$False_Positive == 0,-1,1),
                                  param.kernel = 2,
                                  trace = F,
                                  ridge = 1e-12, lambda.min = 1e-3, eps = 1e-9)

      data.grid <- expand.grid(abs_log2FoldChange = seq(from=0,to=4, length = grid_fc_dens),
                               log10_baseExpr = seq(from=0,to=6, length = grid_ex_dens)) %>%
        dplyr::mutate(Predicted_FP = as.numeric(stats::predict(svm_val, as.matrix(.),
                                                               lambda = tail(svm_val$lambda,1),
                                                               type = "class")),
                      Predicted_FP = ifelse(Predicted_FP == -1,0,1))
    },
    error = function(e) e1071 <<- TRUE
    )

    if(e1071){
      svm_val <- e1071::svm(False_Positive ~ ., data_svm, 
                            kernel = "polynomial", degree = 2, scale = F)

      data.grid <- expand.grid(abs_log2FoldChange = seq(from=0,to=4, length = grid_fc_dens),
                               log10_baseExpr = seq(from=0,to=6, length = grid_ex_dens)) %>%
        dplyr::mutate(Predicted_FP = as.numeric(as.character(stats::predict(svm_val, as.matrix(.)))))
    }

    # 2. Extract and smooth the decision boundary using convex hull points.
    db_sm <- delboy::process_decision_boundary(data.grid)
    
    if(min(db_sm$abs_log2FoldChange, na.rm = T) < 0.05)
      .db_message("low fold-change limit for the false positive decision boundary:\n consider using a higher 'max.iter' value to minimise the risk of false positives", 
                  "blue")
    
    # 3. Predict FP-TP status for validation data.
    pred_fp <- delboy::predict_FP_delboy(data_svm, db_sm)$pred_FP
    data_svm %<>%
      dplyr::mutate(Predicted_FP = pred_fp)

    svm_method <- ifelse(e1071,"e1071","svmpath")
  },
  error = function(e) stop(paste("unable to perform SVM on validation data:",e))
  )
  return(list(svm_validation_fit = svm_val, data_svm = data_svm, grid = data.grid, 
              decision_boundary = db_sm, svm.method = svm_method))
}
