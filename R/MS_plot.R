#' Module significance plot
#'
#' @param met fully preprocessed MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed"
#' @param group_factor name of column in colData for grouping
#' @return Module significance plot
#' @examples
#' MS_plot(met_example, group_factor = "tumor_normal",label_colors=c("darkseagreen","dodgerblue"),p_value_cutoff=0.05)
#' @export
MS_plot = function(met, group_factor,label_colors,p_value_cutoff,p_adjust=TRUE){
    id = grep(group_factor,names(metadata(met)))[2]
    tree = ape::as.phylo(metadata(met)$METree)
    if(p_adjust==TRUE){
        x = -log10(metadata(met)[[id]]$av_adj_pval)
    } else {
        x = -log10(metadata(met)[[id]]$av_pval)
    }
    x[metadata(met)[[id]]$av_fold_change<0] = x[metadata(met)[[id]]$av_fold_change<0]*(-1)
    names(x) = tree$tip.label
    obj=phytools::contMap(tree=tree,
                      x=x,
                      res=200,
                      plot=FALSE,
                      lims=c(round(log10(p_value_cutoff),digits = 2),round(-log10(p_value_cutoff),digits = 2)))
    obj = phytools::setMap(obj, colors=c(label_colors[1],rep("white",3), label_colors[2]))
    plot(obj)
}


