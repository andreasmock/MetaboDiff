#' Identification of metabolic correlation modules
#'
#' @param met MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed" and metadata slot containing the results of the differential analysis
#' @return MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed" and metadata slot containing the results of the differential analysis and dissimilarity matrix
#' @examples
#' identify_modules(met_example)
#' @export
identify_modules <- function(met, min_module_size) {
    w = metadata(met)$diss_matrix
    #create gene tree by average linkage hierarchical clustering
    metadata(met)$tree = hclust(as.dist(w), method = 'average')

    #module identification using dynamic tree cut algorithm
    metadata(met)$modules = as.vector(dynamicTreeCut::cutreeDynamic(dendro =  metadata(met)$tree, distM = w, deepSplit = 4, pamRespectsDendro = FALSE,
                            minClusterSize = min_module_size))
    #assign module colour vector
    metadata(met)$module_color_vector = WGCNA::labels2colors(metadata(met)$modules)

    #calculate eigengenes
    metadata(met)$MEs = WGCNA::moduleEigengenes(t(assay(met[["norm_imputed"]])),
                           colors = metadata(met)$modules,
                           excludeGrey = FALSE)$eigengenes
    colnames(metadata(met)$MEs) = 0:(ncol(metadata(met)$MEs)-1)

    #calculate dissimilarity of module eigengenes
    metadata(met)$MEDiss = 1-cor( metadata(met)$MEs);

    #cluster module eigengenes
    metadata(met)$METree = hclust(as.dist(metadata(met)$MEDiss),
                                  method = 'average');

    #assign module colours
    range_modules = range(metadata(met)$modules)
    metadata(met)$module_colors = WGCNA::labels2colors(c(range_modules[1]:range_modules[2]))[metadata(met)$METree$order]

    #calculate module membership
    metadata(met)$MM = abs(WGCNA::bicor(t(assay(met[["norm_imputed"]])),
                                        metadata(met)$MEs))
    met
}

