#' Normalize metabolomic data by vsn
#'
#' @param @param met MultiAssayExperiment object with slots "raw" and "imputed"
#' @return MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed"
#' @examples
#' normalize_met(met)
#' @export
normalize_met <- function(met) {
    raw_temp = met[["raw"]]
    imputed_temp = met[["imputed"]]
    assay(raw_temp) = justvsn(assay(raw_temp))
    assay(imputed_temp) = justvsn(assay(imputed_temp))


    sampleMap_raw = sampleMap(met)[1:nrow(colData(met)),]
    sampleMap_imputed = sampleMap(met)[1:nrow(colData(met)),]

    sampleMap_raw$assay = droplevels(sampleMap_raw$assay)
    sampleMap_imputed$assay = droplevels(sampleMap_imputed$assay)

    levels(sampleMap_raw$assay) <- "norm"
    levels(sampleMap_imputed$assay) <- "norm_imputed"

    met2 = c(x=met,
            norm=raw_temp,
            sampleMap=sampleMap_raw)

    met3 = c(x=met2,
             norm_imputed=imputed_temp,
             sampleMap=sampleMap_imputed)

    met3
}

