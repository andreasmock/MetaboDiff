#' Annotation using Small Molecule Pathway Database (SMPDB)
#'
#' @param metabolite_ids vector of HMDB, KEGG or ChEBI identifiers.
#' @return dataframe containing SMPDB annotation for metabolite identifiers.
#' @examples
#' get_SMPDBanno(c("HMDB00538", "HMDB00250"))
#' get_SMPDBanno(c("C00002", "C00020"))
#' get_SMPDBanno(c("15422", "16027"))
#'
get_SMPDBanno <- function(rowData,
                          column_kegg_id,
                          column_hmdb_id,
                          column_chebi_id) {
    db = read.csv(system.file("extdata", "metabolites.csv", package = "MetaboDiff"))
    res = matrix(NA,nrow=nrow(rowData),ncol=ncol(db))
    colnames(res) = colnames(db)

    if(!is.na(column_kegg_id)){
        temp1 = as.matrix(db[match(rowData[,column_kegg_id],db$KEGG.ID),])
    }
    if(!is.na(column_hmdb_id)){
        temp2 = as.matrix(db[match(rowData[,column_hmdb_id],db$HMDB.ID),])
    }
    if(!is.na(column_chebi_id)){
        temp3 = as.matrix(db[match(rowData[,column_chebi_id],db$ChEBI.ID),])
    }
    if(exists("temp1")){
        res[is.na(res[,1]),] = temp1[is.na(res[,1]),]
    }
    if(exists("temp2")){
        res[is.na(res[,1]),] = temp2[is.na(res[,1]),]
    }
    if(exists("temp3")){
        res[is.na(res[,1]),] = temp3[is.na(res[,1]),]
    }
    cbind(rowData,res)
}
