% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_bayes_factors_2dkcde.R
\name{calculate_bayes_factors_2dkcde}
\alias{calculate_bayes_factors_2dkcde}
\title{calculate_bayes_factors_2dkcde}
\usage{
calculate_bayes_factors_2dkcde(
  data,
  data_bs,
  fc_column,
  expr_column,
  group_column
)
}
\arguments{
\item{data}{A data frame containing data for which Bayes Factors will be calculated.}

\item{data_bs}{A data frame containing bootstrap-sampled data for which a 2D kernel cumulative distribution estimate will be calculated for the data points in \code{data}.}

\item{fc_column}{A character string naming the log fold change column.}

\item{expr_column}{A character string naming the expression level column.}

\item{group_column}{A character string naming the grouping column containing only \code{True_Positive} and \code{False_Positive} entries.}
}
\value{
A list containing the following elements:
\itemize{
\item \code{kcde_TP} - an object of class \code{kcde} for \code{True_Positives}.
\item \code{kcde_FP} - an object of class \code{kcde} for \code{False_Positives}.
\item \code{cum_prob_BF} - a data frame containing cumulative probability estimates and Bayes Factors for the data in \code{data}.
}
}
\description{
Calculates log Bayes Factors using a 2d kernel cumulative distribution estimate for True and False Positives along the dimensions of log fold change (LFC) and expression level (EL): Pr(LFC, EL | TP)/Pr(LFC, EL | FP)).
}
