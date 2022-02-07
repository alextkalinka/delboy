# Helper function.
.geom_mean <- function(...){
  x <- c(unlist(list(...)))
  return(exp(sum(log(x[x > 0]), na.rm=T) / length(x)))
}


#' calc_median_normalization_size_factors
#'
#' Calculates the size factors used in median normalization for count data: median of the geometric mean of counts across all experiments.
#'
#' @param counts A data frame of counts for each sample in the study (samples as columns, gRNAs as rows).
#'
#' @return A data frame with 1 row and as many columns as `counts` with size factors for each experiment.
#' @author Alex T. Kalinka, \email{alex.kalinka@@cancer.org.uk}
#' @importFrom dplyr mutate summarise select vars
#' @importFrom purrr pmap_dbl
#' @export
calc_median_normalization_size_factors <- function(counts){
  size_factors <- counts %>%
    dplyr::mutate(geom_mean = purrr::pmap_dbl(
      .l = dplyr::select(., -sgRNA, -gene),
      .f = .geom_mean
    )) %>%
    dplyr::mutate_at(dplyr::vars(-sgRNA, -gene, -geom_mean),
                     list(~ . / geom_mean)) %>%
    dplyr::summarise_at(dplyr::vars(-sgRNA, -gene, -geom_mean),
                        median)
  return(size_factors)
}


#' normalize_library_depth_median_ratio
#'
#' Normalizes library depth from a counts file using the median-ratio method [1].
#'
#' @param counts A data frame of counts for each sample in the study (samples as columns, gRNAs as rows).
#'
#' @return A data frame in which the sample columns are library-depth normalized.
#' @author Alex T. Kalinka, \email{alex.kalinka@@cancer.org.uk}
#' @importFrom dplyr mutate
#' @references Wang T, Wei JJ, Sabatini DM, Lander ES. Genetic screens in human cells using the CRISPR-Cas9 system. Science. 2014, 343: 80-84.
#' Anders S, Huber W. Differential expression analysis for sequence count data. Genome Biol. 2010, 11: R106
#' @export
normalize_library_depth_median_ratio <- function(counts){
  tryCatch({
    size_factors <- counts %>%
      fgcQC::calc_median_normalization_size_factors()
    # Make sure none of the size factors are 0.
    size_factors <- apply(size_factors, 2, function(x) ifelse(x==0,1,x))
    counts[,3:ncol(counts)] <- as.data.frame(t(apply(counts[,3:ncol(counts)],
                                                     1, function(x) x/size_factors)))
  },
  error = function(e) stop(paste("unable to calculate median ratio normalization:",e))
  )
  return(counts)
}
