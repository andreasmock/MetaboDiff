#' Boxplot comparing the distribution of variance in metabolic measurements across the pathways
#'
#' @param met MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed"
#' @param rowAnnotation character name of rowAnnotation to stratify metabolites
#' @return Boxplot comparing the distribution of variance in metabolic measurements across the pathways
#' @examples
#' variance_boxplot(met_example, rowAnnotation="SMPDB.Pathway.Name")
#' @export
variance_boxplot <- function(met, rowAnnotation) {
    df = data.frame(value=apply(assays(met)[["norm_imputed"]],1,var),
                     pathway=as.vector(rowData(met[["norm_imputed"]])[[rowAnnotation]]))
    df$pathway_ordered = reorder(df$pathway,df$value,median)
    ggplot(df,
           mapping=aes(x=pathway_ordered,
                       y=value,
                       fill=pathway_ordered)) + coord_flip() + theme_minimal() +
        geom_boxplot() + xlab("") + guides(fill=FALSE) + ylab("vsn normalized abundance")
}
