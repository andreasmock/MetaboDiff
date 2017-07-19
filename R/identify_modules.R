#' Identification of metabolic correlation modules
#'
#' @param met MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed" and metadata slot containing the results of the differential analysis
#' @return MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed" and metadata slot containing the results of the differential analysis and dissimilarity matrix
#' @examples
#' identify_modules(met_example)
identify_modules <- function(met, min_module_size) {
    w = metadata(met)$diss_matrix
    #create gene tree by average linkage hierarchical clustering
    metadata(met)$tree = hclust(as.dist(w), method = 'average')

    #module identification using dynamic tree cut algorithm
    metadata(met)$modules = cutreeDynamic(dendro = tree, distM = w, deepSplit = 4, pamRespectsDendro = FALSE,
                            minClusterSize = min_module_size)
    #assign module colours
    metadata(met)$module_colors = labels2colors(modules)
    met
}

