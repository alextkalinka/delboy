# `delboy`

## Summary

`delboy` is an `R` package for conducting differential-expression analyses on RNA-seq data in which there are exactly **two groups** to be contrasted. The method is designed to improve sensitivity for under-powered data-sets, in which the effect sizes are small and there are few replicates, while controlling the False Discovery Rate (FDR).

`delboy` - **D**ifferential-representation analysis by **E**lastic-net **L**ogistic regression with **B**in**O**mial-thinning validit**Y** tests.

You can read about the method in the companion [manuscript](https://www.biorxiv.org/content/10.1101/2020.10.15.340737v1.full).

## Installation

```r
# install.packages("devtools")
devtools::install_github("alextkalinka/delboy")
```

## Usage

Input data should be a data frame of normalized counts in which there is a gene column with the remaining columns being sample columns.

```r
db <- delboy::run_delboy(
		data = expr_data_frame,
		group_1 = c("ctrl-1","ctrl-2","ctrl-3"),
		group_2 = c("treat-1","treat-2","treat-3"),
		filter_cutoff = 40,
		gene_column = "gene_id",
		batches = NULL
)

# To print a summary report to the console:
db

# To extract a data frame of hits
# (includes a 'Predicted_False_Positive' column):
my_hits <- delboy::hits(db)

```

To plot validation performance next to original data showing the false-positive decision boundary (axes limits can be controlled using the `xlim` and `ylim` arguments):

```r
plot(db, type = "fc_expr")
```

To visualize false negatives in the validation data relative to the false-positive decision boundary:

```r
plot(db, type = "fc_expr_FN")
```

To plot the distrubution of log-fold changes used for the validation data:

```r
plot(db, type = "lfc_nonnull")
```

## References

Kalinka, A. T. (2020). Improving the sensitivity of differential-expression analyses for under-powered RNA-seq experiments. bioRxiv [10.1101/2020.10.15.340737](https://www.biorxiv.org/content/10.1101/2020.10.15.340737v1.full).

## Bugs, Issues, or Requests

Please contact [Alex Kalinka](mailto:alex.t.kalinka@gmail.com).
