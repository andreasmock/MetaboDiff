#' Name metabolic correlation modules
#'
#' @param met MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed" and metadata slot containing the results of the differential analysis
#' @param pathway_annotation Pathway annotation from rowData
#' @return names of metabolic correlation modules stored in metadata(met)$METree$labels
#' @examples
#' name_modules(met_example, pathway_annotation="SMPDB|Pathway.Name")
#' @export
name_modules = function(met, pathway_annotation){
    modules = metadata(met)$modules
    module_names = rep(NA,max(modules)+1)
    res = rep(NA,max(modules)+1)
    for(i in 1:length(module_names)){
        ids = as.numeric(which(modules==(i-1)))
        paths = table(rowData(met[["norm_imputed"]])[[pathway_annotation]][ids])
        names = names(which(paths[paths>0]==max(paths[paths>0])))
        if (length(names)>1){
            module_names[i] = paste0(i-1, " | ", paste(names[1:2],collapse=" | "))
        } else {
            module_names[i] = paste0(i-1, " | ", names)
        }
    }
    metadata(met)$METree$labels = module_names
    met
}

