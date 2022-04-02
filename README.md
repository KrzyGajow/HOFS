
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HOFS

<!-- badges: start -->

<!-- badges: end -->

Higher Order Mutual Information Approximation for Feature Selection.

## Installation

You can install the released version of HOFS from
[GitHub](https://github.com/) with:

``` r
library(remotes)
install_github("KrzyGajow/HOFS")
```

## Example

This is a basic example which shows you how to solve a common problem:


```r
library("HOFS")
data(iris)
Results <- HOFS( data = iris[ , -5 ], label = iris[ , 5, drop = F ], Nfeatures = 4 )

# Run Shiny web application 
runShinyHOFS()
```
