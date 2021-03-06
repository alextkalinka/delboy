% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep_count_matrix.R
\name{prep_count_matrix}
\alias{prep_count_matrix}
\title{prep_count_matrix}
\usage{
prep_count_matrix(data, group_1, group_2, gene_column)
}
\arguments{
\item{data}{A data frame of normalized count data (can contain negative values).}

\item{group_1}{A character vector naming the columns that belong to group 1.}

\item{group_2}{A character vector naming the columns that belong to group 2.}

\item{gene_column}{A character string naming the column containing gene names.}
}
\value{
A matrix of integer counts with any negative values set to zero - rownames equal gene names.
}
\description{
Prepares an integer count matrix from a data frame of normalized counts so the data may be used by functions requiring matrix inputs.
}
