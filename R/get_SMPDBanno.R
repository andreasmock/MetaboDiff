#' Annotation using Small Molecule Pathway Database (SMPDB)
#'
#' @param met MultiAssayExperiment object created with create_mae
#' @param column_kegg_id Column number containing KEGG ids.
#' @param column_hmbd_id Column number containing HMDB ids.
#' @param column_chebi_id Column number containing ChEBI ids.
#' @return MultiAssayExperiment object with SMPDB annotation
#' @examples
#' get_SMPDBanno(met_example,column_kegg_id=6,column_hmdb_id=7,column_chebi_id=NA)
#' @export
get_SMPDBanno <- function(met,
                          column_kegg_id,
                          column_hmdb_id,
                          column_chebi_id) {
    rowData = rowData(met[["raw"]])
    db = read.csv(system.file("extdata", "metabolites.csv", package = "MetaboDiff"))
    res = matrix(NA,nrow=nrow(rowData),ncol=ncol(db))
    colnames(res) = colnames(db)
    colnames(res) = paste0("SMPDB|",colnames(res))

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
    rowData(met[["raw"]]) = data.frame(rowData,res)
    met
}
