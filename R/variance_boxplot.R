#' Boxplot comparing the distribution of variance in metabolic measurements across the pathways
#'
#' @param met MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed"
#' @param rowAnnotation character name of rowAnnotation to stratify metabolites
#' @return Boxplot comparing the distribution of variance in metabolic measurements across the pathways
#' @examples
#' variance_boxplot(met_example, rowAnnotation="Pathway.Name")
variance_boxplot <- function(met, rowAnnotation) {
    mdata = as.data.frame(longFormat(t(assays(met)[["norm_imputed"]]),
                                     colDataCols=))
    ggplot(mdata,
           mapping=aes(x=colname,
                       y=value,
                       fill=rowData(met[["Pathway.Name"]]))) +
        geom_boxplot() + xlab("") + ggtitle("raw") +
        theme(axis.text.x=element_blank()) + ylab("log2(raw abundance)")
}

