#' Module of interest plot
#'
#' @param met fully preprocessed MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed"
#' @param group_factor name of column in colData for grouping
#' @param MOI module of interest
#' @return plot exploring the module of interest
#' @examples
#' MOI_plot(met_example, group_factor = "tumor_normal", MOI=2)
#' @export
MOI_plot = function(met, group_factor, MOI){
    mets = rowData(met[["norm_imputed"]])$BIOCHEMICAL[metadata(met)$modules==MOI]
    x = -log10(metadata(met)[[paste0("ttest_",group_factor)]]$adj_pval[metadata(met)$modules==MOI])
    y =  metadata(met)$MM[,MOI+1][metadata(met)$modules==MOI]
    fc = metadata(met)[[paste0("ttest_",group_factor)]]$fold_change[metadata(met)$modules==MOI]
    df = data.frame(mets=mets,x=x,y=y,fc=fc)
    ggplot(df, aes(x=x,y=y)) + geom_point() + geom_text(aes(label=mets), vjust=1.5) + theme_classic() +
        xlab("p-value (-log10)") + ylab("module membership") + geom_vline(aes(xintercept = 1.30103),colour="darkorange3",linetype="dashed")
}
