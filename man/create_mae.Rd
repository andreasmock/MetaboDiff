% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_mae.R
\name{create_mae}
\alias{create_mae}
\title{Create MultiAssayExperiment object for MetaboDiff}
\usage{
create_mae(assay, rowData, colData)
}
\arguments{
\item{assay}{a matrix containing the relative metabolic measurements}

\item{rowData}{a dataframe containing the available metabolite annotation}

\item{colData}{a dataframe containing sample metadata}
}
\value{
MultiAssayExperiment object with slots "raw"
}
\description{
Create MultiAssayExperiment object for MetaboDiff
}
\examples{
create_mae(assay,rowData,colData)
}
