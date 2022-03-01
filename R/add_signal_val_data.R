#' add_signal_val_data
#'
#' For performance estimation, this function adds logFC signal to a known set of genes and guides in signal-controlled count data.
#'
#' @param data A data frame of normalized counts (at the gRNA level) for a set of samples with the true signal controlled.
#' @param group_1 A character vector naming the columns that belong to group 1.
#' @param group_2 A character vector naming the columns that belong to group 2.
#' @param treat_samps A character vector naming the treatment samples in `data`.
#' @param lfc_samp_grna A named vector of logFC samples, with gRNA ID names.
#' @param grna_column A character string naming the column containing gRNA IDs.
#' 
#' @return A list containing the following components:
#' * `data.bthin`: A data frame containing ground truth signal.
#' * `coef_mat`: A coefficient matrix to be used by `seqgendiff`.
#' * `design_mat`: A design matrix to be used by `seqgendiff`.
#' * `thout`: output from `seqgendiff`.
#' * `group_1`: A character vector of sample names in group_1.
#' * `group_2`: A character vector of sample names in group_2.
#' * 
#' @md
#' @export
#' @importFrom dplyr %>% select
#' @importFrom rlang sym !!
add_signal_val_data <- function(data, group_1, group_2, treat_samps, lfc_samp_grna, grna_column){
  tryCatch({
    # 1. Create coefficient matrix for seqgendiff.
    coef_mat <- delboy::make_coef_matrix(data, lfc_samp_grna, grna_column)
    
    # 2. Create design matrix for seqgendiff.
    design_mat <- delboy::make_design_matrix(group_1, group_2, treat_samps)
    
    # 3. Add signal using seqgendiff's binomial-thinning approach.
    thout <- seqgendiff::thin_diff(mat = data.m,
                                   design_fixed = design_mat,
                                   coef_fixed = coef_mat)
    
    # 4. Prep bthin matrix for use in DiffExp analyses.
    data.bthin <- delboy::prep_bthin_matrix_diffrep(data, thout$mat,
                                                    colnames(data.m),
                                                    as.logical(c(design_mat)),
                                                    grna_column)
    
    # 5. New treatment groups.
    group_1.v <- colnames(data.bthin %>%
                            dplyr::select(- !!rlang::sym(grna_column)))[!as.logical(c(design_mat))]
    group_2.v <- colnames(data.bthin %>%
                            dplyr::select(- !!rlang::sym(grna_column)))[as.logical(c(design_mat))]
  },
  error = function(e) stop(paste("unable to add signal to validation data:",e))
  )
  return(list(data.bthin = data.bthin,
              coef_mat = coef_mat,
              design_mat = design_mat,
              thout = thout,
              group_1 = group_1.v,
              group_2 = group_2.v))
}
