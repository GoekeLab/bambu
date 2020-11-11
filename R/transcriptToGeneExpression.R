#' Reduce transcript expression to gene expression
#' @title transcript to gene expression
#' @param se a summarizedExperiment object from \code{\link{bambu}}
#' @return A SummarizedExperiment object
#' @export
#' @examples
#' se <- readRDS(system.file("extdata",
#'     "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",
#'     package = "bambu"
#' ))
#' transcriptToGeneExpression(se)
transcriptToGeneExpression <- function(se) {
    counts <- as.data.table(assays(se)$counts, keep.rownames = TRUE)
    runnames <- colnames(counts)[-1]
    colnames(counts)[-1] <- rename_duplicatedNames(runnames)
    counts <- melt(counts, id.vars = "rn", measure.vars = colnames(counts)[-1])
    setnames(counts, "rn", "TXNAME")
    rowDataSe <- as.data.table(rowData(se))
    counts <- rowDataSe[, .(TXNAME, GENEID)][counts, on = "TXNAME"]
    counts[, valueGene := sum(value), by = list(variable, GENEID)]
    counts[, valueGeneCPM := valueGene / max(sum(value), 1) * 10^6,
        by = list(variable)]
    ## counts
    counts_gene <- dcast(unique(counts[, .(GENEID, variable, valueGene)]),
        GENEID ~ variable, value.var = "valueGene")
    counts_gene_CPM <- dcast(unique(counts[, .(GENEID, variable,
        valueGeneCPM)]), GENEID ~ variable, value.var = "valueGeneCPM")
    ## geneRanges
    exByGene <- txRangesToGeneRanges(rowRanges(se),
        TXNAMEGENEID_Map = rowDataSe[, .(TXNAME, GENEID)])
    if ("newTxClass" %in% colnames(rowDataSe)) {
        rowDataSe <- rowDataSe[, .(TXNAME, GENEID, newTxClass)]
        rowDataSe[, newGeneClass := ifelse(grepl("ENSG", GENEID),
            "annotation", unique(newTxClass)), by = GENEID]
        mcols(exByGene) <- unique(rowDataSe[, .(GENEID,
            newGeneClass)])[match(names(exByGene), GENEID)]
    }
    ## SE
    counts_gene <- setDF(counts_gene)
    RowNames <- counts_gene$GENEID
    rownames(counts_gene) <- RowNames
    counts_gene_CPM <- setDF(counts_gene_CPM)
    rownames(counts_gene_CPM) <- RowNames
    ColNames <- colnames(counts_gene)[-1]
    seOutput <- SummarizedExperiment(
        assays = SimpleList(
            counts = as.matrix(counts_gene[, -1], ncol = length(ColNames),
                dimnames = list(RowNames, ColNames)),
            CPM = as.matrix(counts_gene_CPM[match(RowNames,
                counts_gene_CPM$GENEID), -1], ncol = length(ColNames),
                dimnames = list(RowNames, ColNames))),
        rowRanges = exByGene[RowNames],
        colData = colData(se))
    return(seOutput)
}

#' rename runnames when there are duplicated names
#' @title rename_duplicatedNames
#' @param runnames sample names
#' @noRd
rename_duplicatedNames <- function(runnames){
    ## rename runnames when duplicated names are found
    if (length(which(duplicated(runnames)))) {
        iter <- 1
        while (length(which(duplicated(runnames)))) {
            if (iter == 1) {
                runnames[which(duplicated(runnames))] <-
                    paste0(runnames[which(duplicated(runnames))], "...", iter)
            } else {
                runnames[which(duplicated(runnames))] <-
                    gsub(paste0("...", iter - 1, "$"), paste0("...", iter),
                    runnames[which(duplicated(runnames))])
            }
            iter <- iter + 1
        }
    }
    return(runnames)
}
