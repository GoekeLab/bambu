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
    counts <- assays(se)$counts
    runnames <- colnames(counts)[-1]
    rowDataSe <- as.data.table(rowData(se))
    
    counts  = fac2sparse(rowData(se)$GENEID) %*% counts
    incompatibleCounts <- metadata(se)$incompatibleCounts
    incompatibleCounts = incompatibleCounts[match(rownames(counts), rownames(incompatibleCounts)),]
    counts = counts + incompatibleCounts
    counts.total = colSums(counts)
    counts.total[counts.total==0] = 1
    counts.CPM = counts/counts.total * 10^6

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
    RowNames <- rownames(counts)
    ColNames <- colnames(counts)
    ColData <- colData(se)
    ColData@rownames <- ColNames
    ColData@listData$name <- ColNames
    seOutput <- SummarizedExperiment(
    assays = SimpleList(counts = counts,
            CPM = counts.CPM),
        rowRanges = exByGene[RowNames],
        colData = ColData)
    
    return(seOutput)
}