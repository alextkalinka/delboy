#' prep_val_data
#'
#' Normalizes and batch-corrects signal in count data in preparation for use in performance evaluation.
#'
#' @param data A data frame of counts (at the gRNA level) for a set of samples.
#' @param group_1 A character vector naming the columns that belong to group 1.
#' @param group_2 A character vector naming the columns that belong to group 2.
#' @param grna_column A character string naming the column containing gene names.
#' @param gene_column A character string naming the column containing gene names.
#' @param bc_method A character string naming the batch correction method: `combat_seq` (default) or `combat`.
#' @param normalize_method A character string naming the read depth normalization method: `relative` (default) or `median_ratio`. `NULL` indicates no normalization.
#'
#' @return A data frame of normalized, signal-corrected data.
#' @export
#' @importFrom dplyr %>% select everything mutate across
#' @importFrom rlang sym !! :=
#' @importFrom magrittr %<>%
prep_val_data <- function(data, group_1, group_2, grna_column, gene_column, bc_method = "combat_seq", normalize_method = "relative"){
  tryCatch({
    data <- data[,c(grna_column, gene_column, group_1, group_2)]
    if(!is.null(normalize_method)){
      if(normalize_method == "relative"){
        data <- delboy::normalize_library_depth_relative(data, NULL)
      }else{
        data <- delboy::normalize_library_depth_median_ratio(data)
      }
    }
    # Need integer counts for combat_seq.
    if(bc_method == "combat_seq"){
      data %<>%
        dplyr::mutate(dplyr::across(c(3:ncol(data)), round))
    }
    # Correct signal associated with known sample groupings.
    gene_sym <- rlang::sym(gene_column)
    grna_sym <- rlang::sym(grna_column)
    batches <- c(rep("a", length(group_1)), rep("b", length(group_2)))
    data <- delboy::batch_correct(data %>% dplyr::select(- (!!gene_sym)), 
                                  batches = batches, grna_column, method = bc_method, parametric = TRUE) %>%
      dplyr::mutate(!!gene_sym := unlist(data[match(!!grna_sym, unlist(data[,grna_column], use.names = F)),gene_column], use.names = F)) %>%
      dplyr::select(!!grna_sym, !!gene_sym, dplyr::everything())
  },
  error = function(e) stop(paste("unable to prep validation data:",e))
  )
  return(data)
}
