#' Heatmap to visualize missing metabolites across the samples
#'
#' @param met MultiAssayExperiment object with slot "raw"
#' @param group_factor name of column in colData for grouping
#' @param label_colors vector of colors for levels of sample_labels
#' @return heatmap visualizing missing metabolites across the samples
#' @examples
#' na_heatmap(met_example, group_factor="tumor_normal", label_colors=c("darkseagreen","dodgerblue"))
#' @export
na_heatmap = function(met, group_factor, label_colors) {
    met_na = is.na(assay(met))*1
    row_na = apply(met_na,1,sum)/ncol(met_na)
    col_na = apply(met_na,2,sum)/nrow(met_na)
    sample_labels = as.factor(colData(met)[[group_factor]])[order(col_na)]
    col_list = list(grouping=label_colors)
    names(col_list$grouping)=levels(droplevels(sample_labels))
    colanno = columnAnnotation(df=data.frame(grouping=sample_labels),
                               col=col_list,
                               barplot=anno_barplot(col_na[order(col_na)],
                                                    gp=gpar(fill="brown",
                                                            col="brown"),
                                                    axis=TRUE,
                                                    border=FALSE),
                               annotation_height=c(1,3))
    rowanno = rowAnnotation(barplot=row_anno_barplot(row_na[order(row_na)],
                                                     gp=gpar(fill="brown",
                                                             col="brown"),
                                                     axis = TRUE,
                                                     border=FALSE),
                            width=unit(1.5,"cm"))

    ha = Heatmap(matrix = met_na[order(row_na),order(col_na)],
                 name="missing",
                 col=c("lightgrey","brown"),
                 cluster_row = FALSE,
                 cluster_columns = FALSE,
                 column_title = "samples",
                 column_title_side = "bottom",
                 show_column_names = FALSE,
                 show_row_names = FALSE,
                 heatmap_legend_param = list(labels=c("no","yes")),
                 row_title= "metabolites",
                 bottom_annotation = colanno,
                 bottom_annotation_height = unit(1.5,"cm"))
    ha + rowanno
}
