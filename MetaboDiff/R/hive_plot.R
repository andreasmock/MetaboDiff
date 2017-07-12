#' Annotation using Small Molecule Pathway Database (SMPDB)
#'
#' @param metabolite_ids vector of HMDB, KEGG or ChEBI identifiers.
#' @return dataframe containing SMPDB annotation for metabolite identifiers.
#' @examples
#' get_SMPDBanno(c("HMDB00538", "HMDB00250"))
#' get_SMPDBanno(c("C00002", "C00020"))
#' get_SMPDBanno(c("15422", "16027"))
#'
#hiveplot <- function(met) {

cd = rowData(met[["norm"]])
df = data.frame(lab1=cd$SUPER_PATHWAY,lab2=cd$SUB_PATHWAY,lab3=cd$BIOCHEMICAL)
lab0 = paste0(levels(df$lab1),"|","Metabolome")
lab12 = paste0(df$lab1,"|",df$lab2)
lab23 = paste0(df$lab2,"|",df$lab3)
df2=strsplit(unique(c(lab0,lab12,lab23)),"[|]")
edge_df = data.frame(node1=sapply(df2,"[",i=1),node2=sapply(df2,"[",i=2))

edge2HPD_mod = function (edge_df= NULL, axis.cols = NULL, type = "2D", desc = NULL,
                         ...){
    if (is.null(edge_df)) {
        stop("No edge data provided")
    }
    if (!is.data.frame(edge_df)) {
        stop("edge_df is not a data frame")
    }
    lab1 <- unlist(edge_df[, 1])
    lab1 <- as.character(lab1)
    lab2 <- unlist(edge_df[, 2])
    lab2 <- as.character(lab2)
    nn <- length(unique(c(lab1, lab2)))
    size <- rep(1, nn)
    id <- 1:nn
    axis <- rep(1, nn)
    color <- as.character(rep("black", nn))
    radius <- rep(1, nn)
    HPD <- list()
    HPD$nodes$id <- id
    HPD$nodes$lab <- unique(c(lab1, lab2))
    HPD$nodes$axis <- axis
    HPD$nodes$radius <- radius
    HPD$nodes$size <- size
    HPD$nodes$color <- color
    ne <- nrow(edge_df)
    edge_df[, 1] <- as.character(edge_df[, 1])
    edge_df[, 2] <- as.character(edge_df[, 2])
    HPD$edges$id1 <- rep(NA, ne)
    HPD$edges$id2 <- rep(NA, ne)
    for (n in 1:ne) {
        pat1 <- paste("", edge_df[n, 1], "", sep = "")
        pat2 <- paste("", edge_df[n, 2], "", sep = "")
        HPD$edges$id1[n] <- which(HPD$nodes$lab==pat1)
        HPD$edges$id2[n] <- which(HPD$nodes$lab==pat2)
    }
    if (ncol(edge_df) > 2) {
        if (is.numeric(edge_df[, 3]) | is.integer(edge_df[, 3])) {
            edge_weight <- edge_df[, 3]
        }
        else {
            warning("No edge weight column detected. Setting default edge weight to 1")
            edge_weight <- rep(1, ne)
        }
    }
    HPD$edges$weight <- edge_weight
    HPD$edges$color <- rep("gray", ne)
    HPD$nodes <- as.data.frame(HPD$nodes)
    HPD$edges <- as.data.frame(HPD$edges)
    if (is.null(desc)) {
        desc <- "No description provided"
    }
    HPD$desc <- desc
    if (is.null(axis.cols)) {
        axis.cols <- brewer.pal(length(unique(HPD$nodes$axis)),
                                "Set1")
    }
    HPD$axis.cols <- axis.cols
    HPD$nodes$axis <- as.integer(HPD$nodes$axis)
    HPD$nodes$size <- as.numeric(HPD$nodes$size)
    HPD$nodes$color <- as.character(HPD$nodes$color)
    HPD$nodes$lab <- as.character(HPD$nodes$lab)
    HPD$nodes$id <- as.integer(HPD$nodes$id)
    HPD$edges$id1 <- as.integer(HPD$edges$id1)
    HPD$edges$id2 <- as.integer(HPD$edges$id2)
    HPD$edges$weight <- as.numeric(HPD$edges$weight)
    HPD$edges$color <- as.character(HPD$edges$color)
    HPD$type <- type
    class(HPD) <- "HivePlotData"
    chkHPD(HPD)
    return(HPD)
}

hive1 = edge2HPD_mod(edge_df,desc = "Test",axis.cols="red")

hive2 <- mineHPD(hive1, option = "rad <- tot.edge.count")
sumHPD(hive2)

hive3 <- mineHPD(hive2, option = "axis <- source.man.sink")
sumHPD(hive3, chk.all = TRUE)

plotHive(hive3,bkgnd="white")

#}

