% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/name_modules.R
\name{name_modules}
\alias{name_modules}
\title{Name metabolic correlation modules}
\usage{
name_modules(met, pathway_annotation)
}
\arguments{
\item{met}{fully preprocessed MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed"}

\item{pathway_annotation}{Pathway annotation from rowData}
}
\value{
names of metabolic correlation modules stored in metadata(met)$METree$labels
}
\description{
Name metabolic correlation modules
}
\examples{
name_modules(met_example, pathway_annotation="SMPDB|Pathway.Name")
}
