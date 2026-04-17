# Core aggregation engine shared by summariseByExon (and any future summariseBy* functions).
# aggregationMat : sparse matrix (output_features x input_transcripts)
# outRanges      : GRanges / GRangesList for the output rows
#' @importFrom SummarizedExperiment assays colData SummarizedExperiment
#' @noRd
.aggregateAssays <- function(se, aggregationMat, outRanges) {
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
    ColData@rownames       <- ColNames
    ColData@listData$name  <- ColNames

    SummarizedExperiment(
        assays    = outAssays,
        rowRanges = outRanges,
        colData   = ColData
    )
}
