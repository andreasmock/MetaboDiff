#' Annotation using Small Molecule Pathway Database (SMPDB)
#'
#' @param metabolite_ids vector of HMDB, KEGG or ChEBI identifiers.
#' @return dataframe containing SMPDB annotation for metabolite identifiers.
#' @examples
#' get_SMPDBanno(c("HMDB00538", "HMDB00250"))
#' get_SMPDBanno(c("C00002", "C00020"))
#' get_SMPDBanno(c("15422", "16027"))
#'
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

