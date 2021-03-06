% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hits.R
\name{hits}
\alias{hits}
\title{hits}
\usage{
hits(delboy_res, remove_predicted_false_positives = FALSE)
}
\arguments{
\item{delboy_res}{An object of class \code{delboy} as produced by \code{delboy::run_delboy}.}

\item{remove_predicted_false_positives}{Logical indicating whether predicted false positives should be removed. Defaults to \code{FALSE}.}
}
\value{
A data frame.
}
\description{
Extract hits from an object of class \code{delboy}.
}
