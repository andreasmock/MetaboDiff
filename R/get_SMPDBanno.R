#' Annotation using Small Molecule Pathway Database (SMPDB)
#'
#' @param metabolite_ids vector of HMDB, KEGG or ChEBI identifiers.
#' @return dataframe containing SMPDB annotation for metabolite identifiers.
#' @examples
#' get_SMPDBanno(c("HMDB00538", "HMDB00250"))
#' get_SMPDBanno(c("C00002", "C00020"))
#' get_SMPDBanno(c("15422", "16027"))
get_SMPDBanno <- function(metabolite_ids) {
    db = read.csv(system.file("extdata", "metabolites.csv", package = "MetaboDiff"))
    if (length(grep(metabolite_id,pattern = "HMDB"))>0) {
        id_type="HMDB.ID"
    }
    else if (length(grep(metabolite_id,pattern = "C"))>0) {
        id_type="KEGG.ID"
    } else {
        id_type="ChEBI.ID"
    }
    which_col = which(colnames(db)==id_type)
    db[match(metabolite_id,db[,which_col]),]
}
