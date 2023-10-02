#' Reduce transcript expression to gene expression
#' @title transcript to gene expression
#' @param se a summarizedExperiment object from \code{\link{bambu}}
#' @return A SummarizedExperiment object
#' @import data.table 
#' @export
#' @examples
#' se <- readRDS(system.file("extdata",
#'     "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",
#'     package = "bambu"
#' ))
#' transcriptToGeneExpression(se)
transcriptToGeneExpression <- function(se) {
    counts <- as.data.table(as.matrix(assays(se)$counts), keep.rownames = TRUE)
    runnames <- colnames(counts)[-1]
    colnames(counts)[-1] <- rename_duplicatedNames(runnames)
    colData(se)@rownames <- rename_duplicatedNames(colData(se)@rownames)
    counts <- melt(counts, id.vars = "rn", measure.vars = colnames(counts)[-1])
    setnames(counts, "rn", "TXNAME")
    rowDataSe <- as.data.table(rowData(se))
    counts <- rowDataSe[, .(TXNAME, GENEID)][counts, on = "TXNAME"]
    
    incompatibleCounts <- metadata(se)$incompatibleCounts
    incompatibleCounts[, TXNAME := "incompatible"]
    counts_incompatible <- melt(incompatibleCounts, id.vars = c("GENEID","TXNAME"), 
        measure.vars = setdiff(colnames(incompatibleCounts), c("GENEID","TXNAME")))
    # GENEID, TXNAME, variable, value
    counts <- rbind(counts, counts_incompatible[variable %in% unique(counts$variable)])
    counts[, valueGene := sum(value), by = list(variable, GENEID)]
    counts[, valueGeneCPM := valueGene / max(sum(value), 1) * 10^6,
           by = list(variable)]

    ## counts
    counts_gene <- dcast(unique(counts[, .(GENEID, variable,
        valueGene)]), GENEID ~ variable, value.var = "valueGene")
    counts_gene_CPM <- dcast(unique(counts[, .(GENEID, variable,
        valueGeneCPM)]), GENEID ~ variable, value.var = "valueGeneCPM")
    ## geneRanges
    exByGene <- reducedRangesByGenes(rowRanges(se))
    if ("txClassDescription" %in% colnames(rowDataSe)) {
        rowDataSe <- rowDataSe[, .(TXNAME, GENEID, txClassDescription)]
        rowDataSe[, newGeneClass := ifelse(grepl("ENSG", GENEID),
            "annotation", unique(txClassDescription)), by = GENEID]
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
    ColData <- colData(se)
    ColData@rownames <- ColNames
    ColData@listData$name <- ColNames
    seOutput <- SummarizedExperiment(
    assays = SimpleList(counts = as.matrix(counts_gene[, -1, drop = FALSE],
            ncol = length(ColNames),
            dimnames = list(RowNames, ColNames)),
            CPM = as.matrix(counts_gene_CPM[match(RowNames,
            counts_gene_CPM$GENEID), -1, drop = FALSE], ncol = length(ColNames),
            dimnames = list(RowNames, ColNames))),
        rowRanges = exByGene[RowNames],
        colData = ColData)
    
    if(is(assays(se)$counts, "sparseMatrix")) {
        assays(seOutput)$counts <- as(assays(seOutput)$counts, "sparseMatrix")
        assays(seOutput)$CPM <- as(assays(seOutput)$CPM, "sparseMatrix")
    }
    
    return(seOutput)
}