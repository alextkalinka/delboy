% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svm_false_positive_classification.R
\name{svm_false_positive_classification}
\alias{svm_false_positive_classification}
\title{svm_false_positive_classification}
\usage{
svm_false_positive_classification(data, kernel)
}
\arguments{
\item{data}{A data frame as produced by \code{delboy::make_delboy_hit_comparison_table}.}

\item{kernel}{A character string naming the kernel type: \code{linear}, \code{polynomial}, \code{radial}.}
}
\value{
A list containing the SVM fit, a data frame containing predicted false positive status, and a data frame produced by \code{expand.grid} for use in \code{ggplot2::geom_contour}.
}
\description{
Performs an SVM to learn the decision boundary between false and true positives for validation data.
}
