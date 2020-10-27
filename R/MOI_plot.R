#' Module of interest plot
#'
#' @param met fully preprocessed MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed"
#' @param group_factor name of column in colData for grouping
#' @param MOI module of interest
#' @return plot exploring the module of interest
#' @examples
#' MOI_plot(met_example, group_factor = "tumor_normal", MOI=2, label_colors=c("darkseagreen","dodgerblue"), p_adjust=TRUE)
#' @export
MOI_plot = function(met, group_factor, MOI, label_colors, p_adjust=TRUE, ...){
    id = grep(group_factor,names(metadata(met)))[1]
    mets = rowData(met[["norm_imputed"]])$BIOCHEMICAL[metadata(met)$modules==MOI]
    if(p_adjust==TRUE){
        x = -log10(metadata(met)[[id]]$adj_pval[metadata(met)$modules==MOI])
        y =  metadata(met)$MM[,MOI+1][metadata(met)$modules==MOI]
        fc = metadata(met)[[id]]$dm[metadata(met)$modules==MOI]
        df = data.frame(mets=mets,x=x,y=y,fc=fc)
        ggplot(df, aes(x=x,y=y,color=fc)) +
            geom_point() +
            scale_color_gradient2(low=label_colors[1],mid="grey",high=label_colors[2],midpoint=0) +
            geom_text(aes(label=mets), vjust=1.5) +
            theme_classic() +
            xlab("adjusted p-value (-log10)") +
            ylab("module membership") +
            xlim(range(x) + c(-0.5,+1)) +
            geom_vline(aes(xintercept = 1.30103),
                       colour="darkorange3",
                       linetype="dashed")
    } else {
        x = -log10(metadata(met)[[id]]$pval[metadata(met)$modules==MOI])
        y =  metadata(met)$MM[,MOI+1][metadata(met)$modules==MOI]
        fc = metadata(met)[[id]]$dm[metadata(met)$modules==MOI]
        df = data.frame(mets=mets,x=x,y=y,fc=fc)
        ggplot(df, aes(x=x,y=y,color=fc)) +
            geom_point() +
            scale_color_gradient2(low=label_colors[1],mid="black",high=label_colors[2],midpoint=0) +
            geom_text(aes(label=mets), vjust=1.5) +
            theme_classic() +
            xlab("p-value (-log10)") +
            ylab("module membership") +
            xlim(range(x) + c(-0.5,+1)) +
            geom_vline(aes(xintercept = 1.30103),
                       colour="darkorange3",
                       linetype="dashed")
    }
}
