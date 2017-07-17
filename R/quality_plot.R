#' Quality plot of processing steps
#'
#' @param met MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed"
#' @param group_factor name of column in colData for grouping
#' @param label_colors vector of colors for levels of sample_label
#' @return Quality plot of processing steps
#' @examples
#' quality_plot(met_example, group_factor="tumor_normal", label_colors=c("darkseagreen","dodgerblue"))
quality_plot <- function(met, group_factor, label_colors) {
    mdata = as.data.frame(longFormat(met[,,1],colDataCols=group_factor))
    mdata$value = log2(mdata$value)
    plot1 = ggplot(mdata,
           mapping=aes(x=colname,
                       y=value,
                       fill=tumor_normal)) +
        geom_boxplot() + xlab("") + ggtitle("raw") +
        theme(axis.text.x=element_blank()) + ylab("log2(raw abundance)") + scale_fill_manual(values=label_colors)

    mdata = as.data.frame(longFormat(met[,,2],colDataCols=group_factor))
    mdata$value = log2(mdata$value)
    plot2 = ggplot(mdata,
                   mapping=aes(x=colname,
                               y=value,
                               fill=tumor_normal)) +
        geom_boxplot() + xlab("") + ggtitle("imputed") +
        theme(axis.text.x=element_blank()) + ylab("log2(imputed raw abundance)") + scale_fill_manual(values=label_colors)

    mdata = as.data.frame(longFormat(met[,,3],colDataCols=group_factor))
    plot3 = ggplot(mdata,
                   mapping=aes(x=colname,
                               y=value,
                               fill=tumor_normal)) +
        geom_boxplot() + xlab("") + ggtitle("norm") +
        theme(axis.text.x=element_blank()) + ylab("vsn normalized abundance") + scale_fill_manual(values=label_colors)

    mdata = as.data.frame(longFormat(met[,,4],colDataCols=group_factor))
    mdata$value = log2(mdata$value)
    plot4 = ggplot(mdata,
                   mapping=aes(x=colname,
                               y=value,
                               fill=tumor_normal)) +
        geom_boxplot() + xlab("") + ggtitle("norm_imputed") +
        theme(axis.text.x=element_blank()) + ylab("vsn normalized and imputed abundance") + scale_fill_manual(values=label_colors)
    plot_grid(plot1, plot2, plot3,plot4, align='h', labels=c('A', 'B','C', 'D'),label_size = 18)
}

