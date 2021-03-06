% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/knn_impute.R
\name{knn_impute}
\alias{knn_impute}
\title{Impute metabolomic data by k-nearest neighbour imputation}
\usage{
knn_impute(met, cutoff)
}
\arguments{
\item{met}{MultiAssayExperiment object with slot "raw"}

\item{cutoff}{Cutoff at which fraction of missing values, a metabolite should be removed from the data. The recommended cutoff is 0.4, e.g. 40 percent.}
}
\value{
MultiAssayExperiment object with slots "raw" and "imputed"
}
\description{
Impute metabolomic data by k-nearest neighbour imputation
}
\examples{
knn_impute(create_mae(assay,rowData,colData), cutoff=0.4)
}
