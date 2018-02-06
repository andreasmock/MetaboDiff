#' pca plot
#'
#' @param met fully preprocessed MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed"
#' @param group_factor name of column in colData for grouping
#' @param label_colors vector of colors for levels of sample_label
#' @return pca plot
#' @examples
#' pca(met_example, group_factor="tumor_normal", label_colors=c("darkseagreen","dodgerblue"))
#' @export
pca_plot <- function(met, group_factor, label_colors) {
    pca = prcomp(t(assay(met[["norm_imputed"]])))
    df = data.frame(pca$x)[,1:2]
    eigs <- pca$sdev^2
    vars = round(eigs[1:2] / sum(eigs),digits=2)*100
    df$grouping = as.vector(colData(met)[[group_factor]])
    p = ggplot(df, aes(x=PC1,y=PC2,colour=grouping)) +
        geom_point() + theme_classic() +
        scale_colour_manual(values=label_colors) +
        xlab(paste0("PC1 ","(",vars[1],"% of variance)")) +
        ylab(paste0("PC2 ","(",vars[2],"% of variance)"))
    return(p)
}

