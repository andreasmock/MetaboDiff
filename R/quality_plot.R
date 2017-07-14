#' Quality plot of normalization
#'
#' @param met MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed"
#' @return Quality plot of normalization
#' @examples
#' quality_plot(met_example)
quality_plot <- function(met, cutoff) {

    mdata = as.data.frame(longFormat(met[,,1],colDataCols="tumor_normal"))
    mdata$value = log2(mdata$value)
    plot1 = ggplot(mdata,
           mapping=aes(x=colname,
                       y=value,
                       fill=tumor_normal)) +
        geom_boxplot() + xlab("") + ggtitle("raw") +
        theme(axis.text.x=element_blank()) + ylab("log2(raw abundance)") + scale_fill_manual(values=c("darkseagreen","dodgerblue"))

    mdata = as.data.frame(longFormat(met[,,2],colDataCols="tumor_normal"))
    mdata$value = log2(mdata$value)
    plot2 = ggplot(mdata,
                   mapping=aes(x=colname,
                               y=value,
                               fill=tumor_normal)) +
        geom_boxplot() + xlab("") + ggtitle("imputed") +
        theme(axis.text.x=element_blank()) + ylab("log2(imputed raw abundance)") + scale_fill_manual(values=c("darkseagreen","dodgerblue"))

    mdata = as.data.frame(longFormat(met[,,3],colDataCols="tumor_normal"))
    plot3 = ggplot(mdata,
                   mapping=aes(x=colname,
                               y=value,
                               fill=tumor_normal)) +
        geom_boxplot() + xlab("") + ggtitle("norm") +
        theme(axis.text.x=element_blank()) + ylab("vsn normalized abundance") + scale_fill_manual(values=c("darkseagreen","dodgerblue"))

    mdata = as.data.frame(longFormat(met[,,4],colDataCols="tumor_normal"))
    mdata$value = log2(mdata$value)
    plot4 = ggplot(mdata,
                   mapping=aes(x=colname,
                               y=value,
                               fill=tumor_normal)) +
        geom_boxplot() + xlab("") + ggtitle("norm_imputed") +
        theme(axis.text.x=element_blank()) + ylab("vsn normalized and imputed abundance") + scale_fill_manual(values=c("darkseagreen","dodgerblue"))
    plot_grid(plot1, plot2, plot3,plot4, align='h', labels=c('A', 'B','C', 'D'),label_size = 18)
}

