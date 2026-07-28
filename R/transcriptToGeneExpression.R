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
    rowDataSe <- as.data.table(rowData(se))

    counts = fac2sparse(factor(rowData(se)$GENEID, levels = unique(rowData(se)$GENEID))) %*% counts
    if (!is.null(metadata(se)$incompatibleCounts)) {
        incompatibleCounts <- metadata(se)$incompatibleCounts
        if ("nonuniqueCounts" %in% names(metadata(se)))
            incompatibleCounts <- incompatibleCounts + metadata(se)$nonuniqueCounts
        incompatibleCounts <- Matrix(incompatibleCounts[match(rownames(counts), rownames(incompatibleCounts)), ], sparse = TRUE)
        counts <- counts + incompatibleCounts
    }
    counts.total = colSums(counts)
    counts.total[counts.total==0] = 1
    counts.CPM = counts %*% Diagonal(x = 1 / counts.total) * 10^6

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
    metadata(seOutput)$seType <- SE_TYPES[["geneCounts"]]
    return(seOutput)
}

#' @title Unique-count SummarizedExperiment from Bambu read-to-transcript assignments
#' @description This function is intended to be used after the transcript
#' discovery and read-to-transcript assignment steps in \code{\link{bambu}} /
#' \code{\link{bambu.singlecell}}. It generates a transcript-level SummarizedExperiment
#' containing raw unique counts (reads uniquely assigned to a single transcript) without
#' EM estimation. This function is useful for highly multiplexed, sparse data (such as
#' single cell and spatial data) where the EM does not have sufficient information to
#' provide accurate transcript expression estimates.
#' @param quantData A list of \code{quantData} objects, one per sample, produced by the
#' read-to-transcript assignment step of \code{\link{bambu.singlecell}}
#' (\code{output = "quantData"}) or the equivalent \code{assignDist = TRUE} run of
#' \code{\link{bambu}}.
#' @param annotations A \code{GRangesList} of transcript annotations matching the ones used to
#' produce \code{quantData}, typically the extended annotations from transcript discovery.
#' @return A \code{SummarizedExperiment} with one row per transcript and one column per cell
#' (or sample). Pass it to \code{\link{transcriptToGeneExpression}} to collapse the unique
#' counts to the gene level. It contains:
#' \describe{
#'     \item{\code{assays(se)$counts}}{a sparse matrix of unique counts, i.e. reads uniquely
#'     assigned to a single transcript.}
#'     \item{\code{rowRanges(se)}}{the transcript \code{annotations} provided in the argument.}
#'     \item{\code{colData(se)}}{per-cell (or per-sample) metadata carried over from
#'     \code{quantData}, such as \code{id}, \code{sampleName}, and \code{barcode}.}
#'     \item{\code{metadata(se)$incompatibleCounts}}{per-gene counts of reads not compatible
#'     with any annotated transcript. \code{\link{transcriptToGeneExpression}} adds these back
#'     into the gene-level counts to give more accurate gene expression estimates, so reads
#'     that cannot be pinned to one transcript still count toward their gene.}
#'     \item{\code{metadata(se)$nonuniqueCounts}}{per-gene counts of reads compatible with more
#'     than one transcript (ambiguous assignments). Like \code{incompatibleCounts}, these are
#'     added back into the gene-level counts by \code{\link{transcriptToGeneExpression}} for
#'     more accurate gene expression estimates, and also indicate how many reads were
#'     ambiguously assigned.}
#'     \item{\code{metadata(se)$seType}}{a label identifying the \code{SummarizedExperiment} object as \code{"uniqueCounts"} type}
#' }
#' @seealso \code{\link{bambu.singlecell}} and \code{\link{bambu}} for producing \code{quantData};
#' \code{\link{transcriptToGeneExpression}} to collapse the result to gene-level counts.
#' @examples
#' ## This works on both bulk and single-cell quantData;
#' ## here we demonstrate with the single-cell test data.
#' rds.dir <- system.file("extdata", "single_cell", package = "bambu")
#' quantData <- readRDS(file.path(rds.dir,
#'     "quantData_GIS_cellMix_HepG2-A549-H9-HEYA8_5primeSingleCell_multisample_chr9_1_1000000.rds"))
#' extendedAnnotations <- readRDS(file.path(rds.dir,
#'     "extendedAnnotations_GIS_cellMix_HepG2-A549-H9-HEYA8_5primeSingleCell_multisample_chr9_1_1000000.rds"))
#' uniqueCountsSe <- getUniqueCountsSe(quantData, extendedAnnotations)
#' @import data.table
#' @export
getUniqueCountsSe <- function(quantData, annotations) {
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
    uniqueCounts <- do.call(cbind, unname(uniqueCountsList))

    incompatibleCounts <- do.call(cbind, unname(lapply(quantData, getIncompatibleCounts)))
    nonuniqueCounts <- do.call(cbind, unname(lapply(quantData, function(x) {
        generateNonUniqueCountMatrix(getReadClassDt(x), annotations, getSampleData(x)$id)
    })))

    colData <- do.call(rbind, unname(lapply(quantData, getSampleData)))

    se <- SummarizedExperiment(assays = SimpleList(counts = uniqueCounts))
    rowRanges(se) <- annotations
    colData(se) <- DataFrame(colData)
    metadata(se)$incompatibleCounts <- incompatibleCounts
    metadata(se)$nonuniqueCounts <- nonuniqueCounts
    metadata(se)$seType <- SE_TYPES[["uniqueCounts"]]
    return(se)
}
