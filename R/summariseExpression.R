#' Summarise transcript expression to exon-level expression
#' @title summarise by exon
#' @param se a \code{SummarizedExperiment} object from \code{\link{bambu}}
#' @return A \code{RangedSummarizedExperiment} with exon-level expression.
#'   Rows are named \code{EX1, EX2, ...} in order of first appearance.
#'   \code{rowRanges} is a \code{GRanges} with one entry per unique exon locus.
#'   \code{mcols} columns:
#'   \describe{
#'     \item{GENEID}{Comma-separated gene IDs for transcripts sharing the exon}
#'     \item{txNames}{Comma-separated transcript IDs that contain this exon}
#'     \item{exonClass}{\code{"annotated"} if the exon coordinates exist in the
#'       reference annotation, \code{"novel"} otherwise}
#'   }
#' @details Counts (and other assays) are summed across all transcripts sharing
#'   the same exon (identical seqnames, start, end, strand). CPM is recomputed
#'   from the aggregated counts. \code{fullLengthCounts} and
#'   \code{uniqueCounts} are summed where present.
#' @importFrom Matrix sparseMatrix
#' @importFrom SummarizedExperiment rowRanges rowData
#' @importFrom GenomicRanges GRanges seqnames start end strand
#' @importFrom IRanges IRanges
#' @importFrom dplyr left_join mutate group_by summarise filter distinct
#' @export
summariseByExon <- function(se) {
    exonRanges <- unlist(rowRanges(se), use.names = TRUE)
    txNames    <- rownames(se)

    # Build data.frame of exon-transcript pairs
    exonDf <- data.frame(
        TXNAME   = names(exonRanges),
        seqnames = as.character(seqnames(exonRanges)),
        start    = start(exonRanges),
        end      = end(exonRanges),
        strand   = as.character(strand(exonRanges)),
        stringsAsFactors = FALSE
    ) |>
        mutate(exon_id = paste(seqnames, start, end, strand, sep = ":"))

    # Attach GENEID and txClassDescription from rowData
    geneDf <- as.data.frame(rowData(se))[, c("GENEID", "txClassDescription"), drop = FALSE]
    geneDf$TXNAME <- rownames(se)
    exonDf <- left_join(exonDf, geneDf, by = "TXNAME")

    # Collapse coordinates and GENEID/txNames per unique exon
    exonMeta <- exonDf |>
        group_by(exon_id) |>
        summarise(
            seqnames = seqnames[1],
            start    = start[1],
            end      = end[1],
            strand   = strand[1],
            GENEID   = paste(sort(unique(GENEID)),  collapse = ","),
            txNames  = paste(sort(unique(TXNAME)),  collapse = ","),
            .groups  = "drop"
        )

    # exonClass: "annotated" if the exon coordinates (seqnames, start, end, strand)
    # exist in the reference annotation, "novel" otherwise.
    annotatedCoords <- exonDf |>
        filter(txClassDescription == "annotation") |>
        distinct(seqnames, start, end, strand) |>
        mutate(exonClass = "annotated")

    exonMeta <- exonMeta |>
        left_join(annotatedCoords, by = c("seqnames", "start", "end", "strand")) |>
        mutate(exonClass = ifelse(is.na(exonClass), "novel", "annotated"))

    # Sparse binary matrix: unique_exons x transcripts
    uniqueExons <- exonMeta$exon_id
    exonNames   <- paste0("EX", seq_along(uniqueExons))
    aggregationMat <- sparseMatrix(
        i    = match(exonDf$exon_id, uniqueExons),
        j    = match(exonDf$TXNAME,  txNames),
        x    = 1L,
        dims = c(length(uniqueExons), length(txNames)),
        dimnames = list(exonNames, txNames)
    )

    # Build output GRanges
    outRanges <- GRanges(
        seqnames = exonMeta$seqnames,
        ranges   = IRanges(start = exonMeta$start, end = exonMeta$end),
        strand   = exonMeta$strand
    )
    names(outRanges) <- exonNames
    mcols(outRanges)$GENEID  <- exonMeta$GENEID
    mcols(outRanges)$txNames <- exonMeta$txNames
    mcols(outRanges)$exonClass <- exonMeta$exonClass

    .aggregateAssays(se, aggregationMat, outRanges)
}
