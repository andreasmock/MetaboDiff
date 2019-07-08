#' Volcano plot of differential analysis
#'
#' @param met fully preprocessed MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed"
#' @param group_factor name of column in colData for grouping
#' @param label_colors vector of colors for levels of sample_label
#' @return Volcano plot of differential analysis
#' @examples
#' volcano_plot(met_example, group_factor="tumor_normal", label_colors=c("darkseagreen","dodgerblue"),p_adjust=TRUE)
#' @export
volcano_plot <- function(met, group_factor, label_colors,p_adjust=TRUE, ...) {

    id = grep(group_factor,names(metadata(met)))[1]
    df=metadata(met)[[id]]
    name = names(metadata(met))[id]
    lv = levels(as.factor(colData(met)[[group_factor]]))
    xlabel = paste0("log2 fold-change"," [",paste(lv,collapse = "_vs_"),"]")

    if(p_adjust==TRUE){
        plot(df$fold_change, -1*log10(df$adj_pval),bty="n",pch=20,
             xlab=xlabel, ylab="adjusted -log10 p-value", ...)
        abline(h=-log10(0.05),lty=2)
        abline(v=log2(1.5),lty=2)
        abline(v=-log2(1.5),lty=2)
        Tsig = df$fold_change<(-log2(1.5))&df$adj_pval<0.05
        points(x = df$fold_change[Tsig],y=-1*log10(df$adj_pval)[Tsig],col=label_colors[1],pch=20)
        Nsig = df$fold_change>(log2(1.5))&df$adj_pval<0.05
        points(x = df$fold_change[Nsig],y=-1*log10(df$adj_pval)[Nsig],col=label_colors[2],pch=20)
    } else {
        plot(df$fold_change, -1*log10(df$pval),bty="n",pch=20,
             xlab=xlabel, ylab="-log10 p-value", ...)
        abline(h=-log10(0.05),lty=2)
        abline(v=log2(1.5),lty=2)
        abline(v=-log2(1.5),lty=2)
        Tsig = df$fold_change<(-log2(1.5))&df$pval<0.05
        points(x = df$fold_change[Tsig],y=-1*log10(df$pval)[Tsig],col=label_colors[1],pch=20)
        Nsig = df$fold_change>(log2(1.5))&df$pval<0.05
        points(x = df$fold_change[Nsig],y=-1*log10(df$pval)[Nsig],col=label_colors[2],pch=20)
    }
}

