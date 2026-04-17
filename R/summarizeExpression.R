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
    uniqueCounts <- assays(se)$uniqueCounts
    rowDataSe <- as.data.table(rowData(se))

    uniqueCounts = fac2sparse(factor(rowData(se)$GENEID, levels = unique(rowData(se)$GENEID))) %*% uniqueCounts
    incompatibleCounts <- metadata(se)$incompatibleCounts
    nonuniqueCounts <- metadata(se)$nonuniqueCounts
    counts = uniqueCounts + incompatibleCounts + nonuniqueCounts
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

#' Generate a SummarizedExperiment of unique counts from quantData
#' @description This function is intended to be used after the transcript
#'   discovery and \code{assignDist} steps in \code{\link{bambu}}. It builds a
#'   transcript-level SummarizedExperiment containing raw unique counts (reads
#'   uniquely assigned to a single transcript) without EM estimation, which can
#'   be passed directly to \code{\link{transcriptToGeneExpression}} to obtain
#'   gene-level counts as uniqueCounts + nonuniqueCounts + incompatibleCounts.
#' @param quantData a list of quantData objects produced by the assignDist step
#' @param annotations a GRangesList of transcript annotations
#' @return A SummarizedExperiment object with \code{assays$uniqueCounts},
#'   \code{metadata$incompatibleCounts}, and \code{metadata$nonuniqueCounts}
#' @import data.table
#' @noRd
generateUniqueCountsSEFromQuantData <- function(quantData, annotations) {
    uniqueCountsList <- lapply(quantData, function(x) {
        readClassDt <- getReadClassDt(x)
        x_filtered <- readClassDt %>% filter(!multi_align & !is.na(eqClass.match))

        uniqueCounts <- if (nrow(x_filtered) == 0) {
            sparseMatrix(i = 1, j = 1, x = 0, dims = c(length(annotations), nrow(getSampleData(x))))
        } else {
            txids <- mcols(annotations)$txid
            i <- rep(match(x_filtered$txid, txids), lengths(x_filtered$columnIds))
            j <- unlist(x_filtered$columnIds)
            x_vals <- unlist(x_filtered$columnCounts)
            sparseMatrix(i = i, j = j, x = x_vals, dims = c(length(annotations), nrow(getSampleData(x))))
        }
        rownames(uniqueCounts) <- names(annotations)
        colnames(uniqueCounts) <- rownames(getSampleData(x))
        return(uniqueCounts)
    })
    uniqueCounts <- do.call(cbind, uniqueCountsList)

    incompatibleCounts <- do.call(cbind, lapply(quantData, getIncompatibleCounts))
    nonuniqueCounts <- do.call(cbind, lapply(quantData, function(x) {
        generateNonUniqueCountMatrix(getReadClassDt(x), annotations, getSampleData(x)$id)
    }))

    colData <- do.call(rbind, lapply(quantData, getSampleData))

    se <- SummarizedExperiment(assays = SimpleList(uniqueCounts = uniqueCounts))
    rowRanges(se) <- annotations
    colData(se) <- DataFrame(colData)
    metadata(se)$incompatibleCounts <- incompatibleCounts
    metadata(se)$nonuniqueCounts <- nonuniqueCounts

    return(se)
}
