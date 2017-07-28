#' Module of interest plot
#'
#' @param met MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed" and metadata slot containing the results of the differential analysis
#' @return MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed" and metadata slot containing the results of the differential analysis and dissimilarity matrix
#' @examples
#' MOI_plot(met_example, group_factor = "tumor_normal", MOI=2)
#' @export
MOI_plot = function(met, group_factor, MOI){
    mets = rowData(met[["norm_imputed"]])$BIOCHEMICAL[metadata(met)$modules==MOI]
    x = -log10(metadata(met)[[paste0("ttest_",group_factor)]]$adj_pval[metadata(met)$modules==MOI])
    y =  metadata(met)$MM[,MOI+1][metadata(met)$modules==MOI]
    fc = metadata(met)[[paste0("ttest_",group_factor)]]$fold_change[metadata(met)$modules==MOI]
    df = data.frame(mets=mets,x=x,y=y,fc=fc)
    googleVis::gvisBubbleChart(df, idvar="mets",xvar="x",yvar="y", sizevar="fc",
         options=list(width=1000,
                      height=1000,
                      vAxis="{title:'module membership'}",
                      hAxis="{title:'p-value (-log10)'}"
                      ))
}
