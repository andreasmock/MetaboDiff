#' Volcano plot of differential analysis
#'
#' @param met fully preprocessed MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed"
#' @param group_factor name of column in colData for grouping
#' @param label_colors vector of colors for levels of sample_label
#' @return Volcano plot of differential analysis
#' @examples
#' volcano_plot(met_example, group_factor="tumor_normal", label_colors=c("darkseagreen","dodgerblue"),p_adjust=TRUE)
#' @export
volcano_plot <- function(met, group_factor, label_colors, dm_cutoff=0.5, p_adjust=TRUE, ...) {

    id = grep(group_factor,names(metadata(met)))[1]
    df=metadata(met)[[id]]
    name = names(metadata(met))[id]
    lv = levels(as.factor(colData(met)[[group_factor]]))
    xlabel = paste0("difference in means"," [",paste(lv,collapse = "-"),"]")

    if(p_adjust==TRUE){
        plot(df$dm, -1*log10(df$adj_pval),bty="n",pch=20,
             xlab=xlabel, ylab="adjusted -log10 p-value", ...)
        abline(h=-log10(0.05),lty=2)
        abline(v=dm_cutoff,lty=2)
        abline(v=-dm_cutoff,lty=2)
        Tsig = df$dm<(-dm_cutoff)&df$adj_pval<0.05
        points(x = df$dm[Tsig],y=-1*log10(df$adj_pval)[Tsig],col=label_colors[1],pch=20)
        Nsig = df$dm>(dm_cutoff)&df$adj_pval<0.05
        points(x = df$dm[Nsig],y=-1*log10(df$adj_pval)[Nsig],col=label_colors[2],pch=20)
    } else {
        plot(df$dm, -1*log10(df$pval),bty="n",pch=20,
             xlab=xlabel, ylab="-log10 p-value", ...)
        abline(h=-log10(0.05),lty=2)
        abline(v=dm_cutoff,lty=2)
        abline(v=-dm_cutoff,lty=2)
        Tsig = df$dm<(-dm_cutoff)&df$pval<0.05
        points(x = df$dm[Tsig],y=-1*log10(df$pval)[Tsig],col=label_colors[1],pch=20)
        Nsig = df$dm>(dm_cutoff)&df$pval<0.05
        points(x = df$dm[Nsig],y=-1*log10(df$pval)[Nsig],col=label_colors[2],pch=20)
    }
}

