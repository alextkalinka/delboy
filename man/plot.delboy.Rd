% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.delboy.R
\name{plot.delboy}
\alias{plot.delboy}
\title{plot.delboy}
\usage{
\method{plot}{delboy}(delboy_res, type, xlim = NULL, ylim = NULL, ...)
}
\arguments{
\item{delboy_res}{Output from \code{delboy::run_delboy}.}

\item{type}{A character string naming the plot type: \code{fc_expr}, \code{fc_expr_FN}, \code{lfc_nonnull} \code{deviance}, \code{misclass}.}

\item{xlim}{xlim values for x-axis. Defaults to \code{NULL} for \code{c(0.5,4)}.}

\item{ylim}{xlim values for y-axis. Defaults to \code{NULL} for \code{c(0,1.5)}.}

\item{...}{Other arguments to be passed to \code{plot}.}
}
\value{
Used for side-effect of plotting.
}
\description{
Plotting for \code{delboy} objects.
}