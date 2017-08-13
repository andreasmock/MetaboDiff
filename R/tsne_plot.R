#' tSNE plot
#'
#' @param met fully preprocessed MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed"
#' @param group_factor name of column in colData for grouping
#' @param label_colors vector of colors for levels of sample_label
#' @return tSNE plot
#' @examples
#' tsne_plot(met_example, group_factor="tumor_normal", label_colors=c("darkseagreen","dodgerblue"))
#' @export
tsne_plot <- function(met, group_factor, label_colors) {
    df = data.frame(tsne::tsne(X = t(assay(met[["norm_imputed"]]))))
    colnames(df) = c("tSNE1","tSNE2")
    df$grouping = as.vector(colData(met)[[group_factor]])
    p = ggplot(df, aes(x=tSNE1,y=tSNE2,colour=grouping)) + geom_point() + theme_classic() + scale_colour_manual(values=label_colors)
    return(p)
}

