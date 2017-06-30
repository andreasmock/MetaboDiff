#' Heatmap to visualize missing metabolites across the samples
#'
#' @param metabolite_ids vector of HMDB, KEGG or ChEBI identifiers.
#' @return dataframe containing SMPDB annotation for metabolite identifiers.
#' @examples
#' get_SMPDBanno(c("HMDB00538", "HMDB00250"))
#' get_SMPDBanno(c("C00002", "C00020"))
#' get_SMPDBanno(c("15422", "16027"))
#'
na_heatmap = function(met){
    met_na = is.na(assay(met))*1
    row_na = apply(met_na,1,sum)/ncol(met_na)
    col_na = apply(met_na,2,sum)/nrow(met_na)

    rowanno = rowAnnotation(barplot=row_anno_barplot(row_na[order(row_na)],
                                                     gp=gpar(fill="brown",
                                                             col="brown"),
                                                     axis = TRUE,
                                                     border=FALSE),
                            width=unit(1.5,"cm"))
    colanno = columnAnnotation(barplot=anno_barplot(col_na[order(col_na)],
                                                    gp=gpar(fill="brown",
                                                            col="brown"),
                                                    axis=TRUE,
                                                    border=FALSE))

    ha = Heatmap(matrix = met_na[order(row_na),order(col_na)],
                 name="missing",
                 col=c("lightgrey","brown"),
                 cluster_row = FALSE,
                 cluster_columns = FALSE,
                 column_title = "samples",
                 column_title_side = "bottom",
                 show_column_names = FALSE,
                 heatmap_legend_param = list(labels=c("no","yes")),
                 row_title= "metabolites",
                 bottom_annotation = colanno,
                 bottom_annotation_height = unit(1.5,"cm"))
    ha + rowanno
}
