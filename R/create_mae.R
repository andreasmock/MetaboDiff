#' Create MultiAssayExperiment object for MetaboDiff
#'
#' @param assay a matrix containing the relative metabolic measurements
#' @param rowData a dataframe containing the available metabolite annotation
#' @param colData a dataframe containing sample metadata
#' @return MultiAssayExperiment object with slots "raw"
#' @examples
#' create_mae(assay,rowData,colData)
#' @export
create_mae = function(assay,rowData,colData){
    rownames(colData) = colnames(assay)
    se = SummarizedExperiment(assays=as.matrix(assay),
                              rowData=rowData)
    experiment_list = list(raw=se)
    sampleMap = data.frame(primary=rownames(colData),
                           colname=colnames(se))
    sampleMap_list = listToMap(list(raw=sampleMap))
    met = MultiAssayExperiment(experiments = experiment_list,
                               colData = colData,
                               sampleMap = sampleMap_list)
    met
}
