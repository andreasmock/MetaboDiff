#' Heatmap to visualize outliers in the study set
#'
#' @param met MultiAssayExperiment object with slots "raw" and "imputed"
#' @param sample_labels factor of sample labels
#' @param label_colors vector of colors for levels of sample_label
#' @return heatmap to visualize outliers in the study set
#' @examples
#' outlier_heatmap(met_example, colData(met_example)$tumor_normal, c("darkseagreen","dodgerblue"))
outlier_heatmap = function(met, sample_labels,label_colors) {
    sample_labels = as.factor(sample_labels)
    met_na = is.na(assay(met))*1
    col_na = apply(met_na,2,sum)/nrow(met_na)
    mat = log2(assays(met)[["imputed"]] / apply(assays(met)[["imputed"]],1,median))
    clusters = cutree(hclust(dist(t(mat))),k=2)
    col_list = list(grouping=label_colors,cluster=c("black","grey"))
    names(col_list$grouping)=levels(sample_labels)
    names(col_list$cluster)=c(1,2)
    rowanno = columnAnnotation(df=data.frame(grouping=colData$tumor_normal,cluster=as.vector(clusters)),
                               col=col_list,
                               missing_metabolites=anno_barplot(col_na[order(col_na)],
                                                    gp=gpar(fill="brown",
                                                            col="brown"),
                                                    axis=TRUE,
                                                    border=FALSE),
                               annotation_height=c(1,1,5),
                               show_annotation_name = TRUE)
    Heatmap(mat,
            show_column_names = FALSE,
            show_row_names = FALSE,
            name="log2 rel. abundance",
            top_annotation = rowanno,
            clustering_distance_rows = "euclidean",
            column_title = "samples",
            column_title_side = "bottom",
            row_title= "metabolites")
}