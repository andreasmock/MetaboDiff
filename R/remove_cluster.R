#' Remove cluster from metabolomic data
#'
#' @param met MultiAssayExperiment object with slots "raw" and "impute"
#' @param cluster numeric vector of cluster to remove
#' @return MultiAssayExperiment object where cluster x is removed from slots "raw" and "impute"
#' @examples
#' remove_cluster(met_example,cluster=1)
#' @export
remove_cluster = function(met, cluster){
    mat = log2(assays(met)[["imputed"]] / apply(assays(met)[["imputed"]],1,median))
    clustering =  cutree(hclust(dist(t(mat))),k=2)
    met[["imputed"]] = met[["imputed"]][,!clustering==cluster]
    met[["raw"]] = met[["raw"]][,!clustering==cluster]
    met
}

