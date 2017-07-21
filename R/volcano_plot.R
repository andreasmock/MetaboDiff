#' Volcano plot of differential analysis
#'
#' @param met MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed" and metadata slot containing the results of the differential analysis
#' @param group_factor name of column in colData for grouping
#' @param label_colors vector of colors for levels of sample_label
#' @return Volcano plot of differential analysis
#' @examples
#' volcano_plot(met_example, group_factor="tumor_normal", label_colors=c("darkseagreen","dodgerblue"))
#' @export
volcano_plot <- function(met, group_factor, label_colors) {
    df=metadata(met)[[paste0("ttest_", group_factor)]]
    lv = levels(as.factor(colData(met)[[group_factor]]))
    xlabel = paste0("log2 fold-change"," [",lv[2]," vs. ",lv[1], "]")
    plot(df$fold_change, -1*log10(df$adj_pval),bty="n",pch=20,
         xlab=xlabel, ylab="adjusted -log10 p-value")
    abline(h=-log10(0.05),lty=2)
    abline(v=log2(1.5),lty=2)
    abline(v=-log2(1.5),lty=2)
    Tsig = df$fold_change<(-log2(1.5))&df$adj_pval<0.05
    points(x = df$fold_change[Tsig],y=-1*log10(df$adj_pval)[Tsig],col=label_colors[1],pch=20)
    Nsig = df$fold_change>(log2(1.5))&df$adj_pval<0.05
    points(x = df$fold_change[Nsig],y=-1*log10(df$adj_pval)[Nsig],col=label_colors[2],pch=20)
}

