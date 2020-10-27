#' Module significance plot
#'
#' @param met fully preprocessed MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed"
#' @param group_factor name of column in colData for grouping
#' @return Module significance plot
#' @examples
#' MS_plot(met_example, group_factor = "tumor_normal",p_value_cutoff=0.05, p_adjust=FALSE)
#' @export
MS_plot = function(met, group_factor,p_value_cutoff,p_adjust=FALSE){
    id = grep(group_factor,names(metadata(met)))[2]
    tree = ape::as.phylo(metadata(met)$METree)
    if(p_adjust==TRUE){
        x = -log10(metadata(met)[[id]]$av_adj_pval)
    } else {
        x = -log10(metadata(met)[[id]]$av_pval)
    }

    names(x) = tree$tip.label
    obj=phytools::contMap(tree=tree,
                      x=x,
                      res=400,
                      plot=FALSE,
                      lims=c(0,-log10(p_value_cutoff))
                      )
    obj = phytools::setMap(obj, colors=c(rep("white",9),"brown"))
    plot(obj)
}


