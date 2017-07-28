#' Module significance plot
#'
#' @param met MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed" and metadata slot containing the results of the differential analysis
#' @return MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed" and metadata slot containing the results of the differential analysis and dissimilarity matrix
#' @examples
#' MS_plot(met_example, group_factor = "tumor_normal")
#' @export
MS_plot = function(met, group_factor){
    tree = ape::as.phylo(metadata(met)$METree)
    x = metadata(met)[[paste0("MS_",group_factor)]]$av_adj_pval
    names(x) = tree$tip.label
    obj=phytools::contMap(tree=tree,
                      x=x,
                      res=200,
                      plot=FALSE,
                      lims=c(0,1))
    obj = phytools::setMap(obj, colors=c("red",rep("lightgreen",3)))
    plot(obj)
}
