% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/volcano_plot.R
\name{volcano_plot}
\alias{volcano_plot}
\title{Volcano plot of differential analysis}
\usage{
volcano_plot(met, group_factor, label_colors, p_adjust = TRUE, ...)
}
\arguments{
\item{met}{fully preprocessed MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed"}

\item{group_factor}{name of column in colData for grouping}

\item{label_colors}{vector of colors for levels of sample_label}
}
\value{
Volcano plot of differential analysis
}
\description{
Volcano plot of differential analysis
}
\examples{
volcano_plot(met_example, group_factor="tumor_normal", label_colors=c("darkseagreen","dodgerblue"),p_adjust=TRUE)
}
