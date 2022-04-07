#' batch_correct
#'
#' Performs a batch correction on normalized count data using `sva`'s `ComBat` method.
#'
#' @param data A data frame of normalized count data.
#' @param batches A character vector naming the batches that the samples belong to.
#' @param gene_column A character string naming the column containing gene names.
#' @param method A character string naming the batch correction method: `combat_seq` or `combat`.
#' @param parametric Logical, whether to perform a parametric correction or not (defaults to `FALSE`).
#'
#' @return A data frame of batch-corrected data (can contain negative values).
#' @importFrom sva ComBat
#' @importFrom dplyr %>% mutate select everything
#' @importFrom rlang := !! sym
#' @export
batch_correct <- function(data, batches, gene_column, method, parametric = FALSE){
  tryCatch({
    # 1. Convert to matrix for ComBat.
    data.m <- as.matrix(data %>%
                          dplyr::select(- !!rlang::sym(gene_column)))
    rownames(data.m) <- data[,gene_column]

    # 2. Batch correct using ComBat.
    if(method == "combat" & !parametric){
      data.bc <- sva::ComBat(data.m, batch = factor(batches), par.prior = FALSE)
    }else if(method == "combat" & parametric){
      data.bc <- sva::ComBat(data.m, batch = factor(batches), par.prior = TRUE)
    }else if(method == "combat_seq"){
      data.bc <- sva::ComBat_seq(data.m, batch = factor(batches))
    }else{
      stop(paste("unrecognised batch correction method:",method))
    }

    # 3. Convert back to data frame.
    data.bc <- as.data.frame(data.bc) %>%
      dplyr::mutate(!!rlang::sym(gene_column) := rownames(.)) %>%
      dplyr::select(!!rlang::sym(gene_column), dplyr::everything())
    
    rownames(data.bc) <- NULL
  },
  error = function(e) stop(paste("unable to batch correct data:",e))
  )
  return(data.bc)
}
