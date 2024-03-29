% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Tday35}
\alias{Tday35}
\title{Tday35}
\description{
A sample of signal-corrected expression data from a bulk RNA organoid experiment with signal added to a known set of genes.
}
\details{
A list with the following elements:
\itemize{
\item \code{data_bcorr}: A data frame of batch-corrected counts.
\item \code{data_val_bcorr}: A data frame of signal-corrected counts.
\item \code{lfc_samp}: A named vector of log-fold changes added to \code{data_bcorr} in which the names are gene names.
\item \code{group_1}: The samples that belong to group 1.
\item \code{group_2}: The samples that belong to group 2.
}
}
\references{
Kelava et al. 2021. Male sex hormones increase production of excitatory neurons in human cerebral organoids. bioRxiv. \href{https://doi.org/10.1101/2020.10.24.353359}{10.1101/2020.10.24.353359}
}
\keyword{data}
