#' Annotation using Small Molecule Pathway Database (SMPDB)
#'
#' @param metabolite_ids vector of HMDB, KEGG or ChEBI identifiers.
#' @return dataframe containing SMPDB annotation for metabolite identifiers.
#' @examples
#' get_SMPDBanno(c("HMDB00538", "HMDB00250"))
#' get_SMPDBanno(c("C00002", "C00020"))
#' get_SMPDBanno(c("15422", "16027"))
#'
normalize_met <- function(met) {
    raw_temp = met[["raw"]]
    imputed_temp = met[["imputed"]]
    assay(raw_temp) = justvsn(assay(raw_temp))
    assay(imputed_temp) = justvsn(assay(imputed_temp))


    sampleMap_raw = sampleMap(met)[1:86,]
    sampleMap_imputed = sampleMap(met)[1:86,]

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

