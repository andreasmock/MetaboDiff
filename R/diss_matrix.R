#' Construct dissimilarity matrix for metabolic correlation network
#'
#' @param met MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed" and metadata slot containing the results of the differential analysis
#' @return MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed" and metadata slot containing the results of the differential analysis and dissimilarity matrix
#' @examples
#' diss_matrix(met_example)
#' @export
diss_matrix <- function(met) {
    mat = t(assay(met[["norm_imputed"]]))
    s = abs(bicor(mat))
    beta = 3
    a = s^beta
    w = 1-a
    metadata(met)[["diss_matrix"]] = w
    met
}

