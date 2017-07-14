#' Impute metabolomic data by k-nearest neighbour imputation
#'
#' @param met MultiAssayExperiment object with slots "raw"
#' @param cutoff Cutoff at which fraction of missing values, a metabolite should be removed from the data. The recommended cutoff is 0.4, e.g. 40%.
#' @return MultiAssayExperiment object with splots "raw" and "imputed
#' @examples
#' knn_impute(met_example,cutoff=0.4)
knn_impute <- function(met, cutoff) {
    met_temp = met[["raw"]]
    met_temp=met_temp[!rowSums(is.na(assay(met_temp)))>(ncol(assay(met_temp))*cutoff),]
    assay(met_temp) = impute.knn(assay(met_temp))$data
    sampleMap_imp = sampleMap(met)
    levels(sampleMap_imp$assay) <- "imputed"
    met2 = c(x=met,
            imputed=met_temp,
            sampleMap=sampleMap_imp)
    met2
}

