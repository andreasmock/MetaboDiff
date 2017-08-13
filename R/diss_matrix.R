#' Construct dissimilarity matrix for metabolic correlation network
#'
#' @param met fully preprocessed MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed"
#' @return adds dissimilarity matrix to metadata slot
#' @examples
#' diss_matrix(met_example)
#' @export
diss_matrix <- function(met) {
    mat = t(assay(met[["norm_imputed"]]))
    s = abs(WGCNA::bicor(mat))
    beta = 3
    a = s^beta
    w = 1-a
    metadata(met)[["diss_matrix"]] = w
    met
}

