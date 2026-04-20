#' Summarise transcript expression to a specified unit
#' @title summarise by expression
#' @param se a \code{SummarizedExperiment} object from \code{\link{bambu}}
#' @param type character, the unit to summarise by. Currently supports
#'   \code{"exon"} (default).
#' @return A \code{RangedSummarizedExperiment} with expression summarised to
#'   the requested unit. Assays \code{counts} and \code{CPM} are always
#'   present; \code{fullLengthCounts} and \code{uniqueCounts} are included
#'   when present in the input. For \code{type = "exon"}, rows are named
#'   \code{EX1, EX2, ...} and \code{rowRanges} mcols include \code{GENEID},
#'   \code{txNames}, and \code{exonClass} (\code{"annotated"} or
#'   \code{"novel"}).
#' @details Counts are summed across all transcripts belonging to the same
#'   group (e.g. identical exon locus). CPM is recomputed from the aggregated
#'   counts. An exon is \code{"novel"} only if all contributing transcripts
#'   have \code{novelTranscript = TRUE}.
#' @examples
#' se <- readRDS(system.file("extdata",
#'     "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",
#'     package = "bambu"
#' ))
#' summariseByExpression(se, type = "exon")
#' @importFrom Matrix sparseMatrix
#' @importFrom S4Vectors SimpleList
#' @importFrom SummarizedExperiment rowRanges rowData assays colData SummarizedExperiment
#' @importFrom GenomicRanges GRanges seqnames start end strand
#' @importFrom IRanges IRanges
#' @importFrom dplyr left_join mutate group_by summarise cur_group_id ungroup
#' @export
summariseByExpression <- function(se, type = "exon") {
    index <- switch(type,
        exon = .buildExonIndex(se),
        stop("Unsupported type: '", type, "'")
    )

    txNames <- rownames(se)
    nGroups <- length(index$outRanges)
    groupNames <- names(index$outRanges)

    aggregationMat <- sparseMatrix(
        i    = index$indexDf$group_idx,
        j    = match(index$indexDf$TXNAME, txNames),
        x    = 1L,
        dims = c(nGroups, length(txNames)),
        dimnames = list(groupNames, txNames)
    )

    counts <- aggregationMat %*% assays(se)$counts

    counts.total <- colSums(counts)
    counts.total[counts.total == 0] <- 1
    cpm <- counts / counts.total * 10^6

    outAssays <- SimpleList(counts = counts, CPM = cpm)

    for (nm in c("fullLengthCounts", "uniqueCounts")) {
        if (nm %in% names(assays(se))) {
            outAssays[[nm]] <- aggregationMat %*% assays(se)[[nm]]
        }
    }

    ColNames <- colnames(counts)
    ColData  <- colData(se)
    ColData@rownames      <- ColNames
    ColData@listData$name <- ColNames

    SummarizedExperiment(
        assays    = outAssays,
        rowRanges = index$outRanges,
        colData   = ColData
    )
}

#' @param se a \code{SummarizedExperiment} object from \code{\link{bambu}}
#' @return a data.frame with columns TXNAME, seqnames, start, end, strand,
#'   GENEID, and novelTranscript — one row per exon-transcript pair.
#' @noRd
.buildFlatDf <- function(se) {
    exonRanges <- unlist(rowRanges(se), use.names = TRUE)
    geneDf <- as.data.frame(rowData(se))[, c("GENEID", "novelTranscript"), drop = FALSE]
    geneDf$TXNAME <- rownames(se)

    data.frame(
        TXNAME   = names(exonRanges),
        seqnames = as.character(seqnames(exonRanges)),
        start    = start(exonRanges),
        end      = end(exonRanges),
        strand   = as.character(strand(exonRanges)),
        stringsAsFactors = FALSE
    ) %>%
        left_join(geneDf, by = "TXNAME")
}

#' @param se a \code{SummarizedExperiment} object from \code{\link{bambu}}
#' @return a list with:
#'   \code{indexDf} — data.frame with columns TXNAME and group_idx
#'     (one row per exon-transcript pair);
#'   \code{outRanges} — GRanges with one entry per unique exon locus,
#'     named EX1, EX2, ..., with mcols GENEID, txNames, and exonClass.
#' @noRd
.buildExonIndex <- function(se) {
    flatDf <- .buildFlatDf(se) %>%
        group_by(seqnames, start, end, strand) %>%
        mutate(group_idx = cur_group_id()) %>%
        ungroup()

    groupMeta <- flatDf %>%
        group_by(group_idx) %>%
        summarise(
            seqnames  = seqnames[1],
            start     = start[1],
            end       = end[1],
            strand    = strand[1],
            GENEID    = paste(sort(unique(GENEID)), collapse = ","),
            txNames   = paste(sort(unique(TXNAME)), collapse = ","),
            exonClass = ifelse(all(novelTranscript), "novel", "annotated"),
            .groups   = "drop"
        )

    outRanges <- GRanges(
        seqnames = groupMeta$seqnames,
        ranges   = IRanges(start = groupMeta$start, end = groupMeta$end),
        strand   = groupMeta$strand
    )
    names(outRanges) <- paste0("EX", seq_len(nrow(groupMeta)))
    mcols(outRanges)$GENEID    <- groupMeta$GENEID
    mcols(outRanges)$txNames   <- groupMeta$txNames
    mcols(outRanges)$exonClass <- groupMeta$exonClass

    list(
        indexDf   = flatDf[, c("TXNAME", "group_idx")],
        outRanges = outRanges
    )
}