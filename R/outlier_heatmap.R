#' Heatmap to visualize outliers in the study set
#'
#' @param met MultiAssayExperiment object with slots "raw" and "imputed"
#' @param group_factor name of column in colData for grouping
#' @param label_colors vector of colors for levels of sample_label
#' @return heatmap to visualize outliers in the study set
#' @examples
#' outlier_heatmap(met_example, group_factor="tumor_normal", label_colors=c("darkseagreen","dodgerblue"), k=2)
#' @export
outlier_heatmap = function(met, group_factor, label_colors, k) {
    colors = c("black","grey","darkmagenta","blue3","gold3","firebrick2","yellow","brown","orange","green")
    sample_labels = as.factor(colData(met)[[group_factor]])
    met_na = is.na(assay(met))*1
    col_na = apply(met_na,2,sum)/nrow(met_na)
    mat = log2(assays(met)[["imputed"]] / apply(assays(met)[["imputed"]],1,median))
    clusters = cutree(hclust(dist(t(mat))),k=k)
    col_list = list(grouping=label_colors,cluster=colors[1:k])
    names(col_list$grouping)=levels(droplevels(sample_labels))
    names(col_list$cluster)=c(1:k)
    rowanno = columnAnnotation(df=data.frame(grouping=colData(met)[[group_factor]],cluster=as.vector(clusters)),
                               col=col_list,
                               missing_metabolites=anno_barplot(col_na[order(col_na)],
                                                    gp=gpar(fill="brown",
                                                            col="brown"),
                                                    axis=TRUE,
                                                    border=FALSE),
                               annotation_height=c(1,1,5),
                               show_annotation_name = TRUE,
                               annotation_name_gp=gpar(cex=0.7))
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
