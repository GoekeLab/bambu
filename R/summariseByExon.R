#' Summarise transcript expression to exon-level expression
#' @title summarise by exon
#' @param se a \code{SummarizedExperiment} object from \code{\link{bambu}}
#' @return A \code{RangedSummarizedExperiment} with exon-level counts
#' @details Counts are summed across all transcripts that share the same exon
#'   (defined by identical seqnames, start, end, and strand). The returned
#'   counts therefore represent the total evidence attributed to each unique
#'   exonic locus across all overlapping transcripts.
#' @import data.table
#' @importFrom Matrix sparseMatrix
#' @importFrom SummarizedExperiment assays rowRanges rowData colData SummarizedExperiment
#' @importFrom GenomicRanges GRanges seqnames start end strand
#' @importFrom IRanges IRanges
#' @export
summariseByExon <- function(se) {
    # Unlist GRangesList: one row per exon-transcript combination
    exonRanges <- unlist(rowRanges(se), use.names = TRUE)

    # Build a data.table of exon-transcript pairs
    txNames <- rownames(se)
    exonDt <- data.table(
        TXNAME   = names(exonRanges),
        seqnames = as.character(seqnames(exonRanges)),
        start    = start(exonRanges),
        end      = end(exonRanges),
        strand   = as.character(strand(exonRanges))
    )

    # Unique exon key: seqnames:start:end:strand
    exonDt[, exon_id := paste(seqnames, start, end, strand, sep = ":")]

    # Attach GENEID from rowData
    geneDt <- data.table(
        TXNAME = rownames(se),
        GENEID = rowData(se)$GENEID
    )
    exonDt <- geneDt[exonDt, on = "TXNAME"]

    # Collapse metadata per unique exon
    exonMeta <- exonDt[, .(
        seqnames = seqnames[1],
        start    = start[1],
        end      = end[1],
        strand   = strand[1],
        GENEID   = paste(sort(unique(GENEID)), collapse = ",")
    ), by = exon_id]

    # Build sparse binary matrix: unique_exons x transcripts
    # entry [i, j] = 1 if transcript j contains unique exon i
    uniqueExons <- exonMeta$exon_id
    exonIdx <- match(exonDt$exon_id, uniqueExons)
    txIdx   <- match(exonDt$TXNAME,  txNames)

    exonTxMat <- sparseMatrix(
        i    = exonIdx,
        j    = txIdx,
        x    = 1L,
        dims = c(length(uniqueExons), length(txNames)),
        dimnames = list(uniqueExons, txNames)
    )

    # Aggregate counts: unique_exons x samples
    txCounts   <- assays(se)$counts
    exonCounts <- exonTxMat %*% txCounts

    # Build GRanges for unique exons
    exonGRanges <- GRanges(
        seqnames = exonMeta$seqnames,
        ranges   = IRanges(start = exonMeta$start, end = exonMeta$end),
        strand   = exonMeta$strand
    )
    names(exonGRanges) <- exonMeta$exon_id
    mcols(exonGRanges)$GENEID <- exonMeta$GENEID

    # Return as SummarizedExperiment
    return(SummarizedExperiment(
        assays   = list(counts = exonCounts),
        rowRanges = exonGRanges,
        colData  = colData(se)
    ))
}
