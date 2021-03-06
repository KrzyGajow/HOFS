% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HOFS.R
\name{HOFS}
\alias{HOFS}
\title{Higher Order Mutual Information Approximation for Feature Selection}
\usage{
HOFS( data, label, C = 0.2, Nfeatures = 10, samplePct = 0.1, verbose = F, seed = 666 )
}
\arguments{
\item{data}{Data.frame with the invesigated features.}
\item{label}{Column with the target feature.}
\item{C}{Predefined constant C, assessing whether create new subset.}
\item{Nfeatures}{Nunber of features to be selected.}
\item{samplePct}{Percentage of observations that are sampled to speed up the computation.}
\item{verbose}{Switch for tracing selection.}
\item{seed}{Seed for Pseudo-Random Number Generator.}
}


\value{
List with the selected features along with the coresponding rank.
}

\examples{
library("HOFS")
data(iris)
Results <- HOFS( data = iris[ , -5 ], label = iris[ , 5, drop = F ], Nfeatures = 4 )
}
