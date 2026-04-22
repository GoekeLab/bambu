#' Summarise transcript expression to a specified unit
#' @title summarise by expression
#' @param se a \code{SummarizedExperiment} object from \code{\link{bambu}}
#' @param type character, the unit to summarise by. Currently supports
#'   \code{"exon"} (default) and \code{"gene"}.
#' @return A \code{RangedSummarizedExperiment} with expression summarised to
#'   the requested unit. Assays \code{counts} and \code{CPM} are always
#'   present; \code{fullLengthCounts} and \code{uniqueCounts} are included
#'   when present in the input. For \code{type = "exon"}, rows are named
#'   \code{EX1, EX2, ...} and \code{rowRanges} mcols include \code{GENEID},
#'   \code{txNames}, and \code{exonClass} (\code{"annotated"} or
#'   \code{"novel"}). For \code{type = "gene"}, rows are named by GENEID and
#'   \code{rowRanges} is a \code{GRangesList} of reduced exon ranges per gene;
#'   mcols include \code{GENEID} and \code{newGeneClass}.
#' @details Counts are summed across all transcripts belonging to the same
#'   group (e.g. identical exon locus or gene). CPM is recomputed from the
#'   aggregated counts. For \code{type = "exon"}, an exon is \code{"novel"}
#'   only if all contributing transcripts have \code{novelTranscript = TRUE}.
#'   For \code{type = "gene"}, \code{newGeneClass} is derived from
#'   \code{txClassDescription}: Ensembl genes (GENEID matching \code{"ENSG"})
#'   are labelled \code{"annotation"}, others inherit the transcript class.
#'   Incompatible and non-unique counts stored in \code{metadata(se)} are
#'   added to the aggregated counts before CPM computation.
#' @examples
#' se <- readRDS(system.file("extdata",
#'     "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",
#'     package = "bambu"
#' ))
#' summariseByExpression(se, type = "exon")
#' summariseByExpression(se, type = "gene")
#' @importFrom Matrix sparseMatrix Matrix
#' @importFrom S4Vectors SimpleList metadata
#' @importFrom SummarizedExperiment rowRanges rowData assays colData SummarizedExperiment
#' @importFrom GenomicRanges GRanges seqnames start end strand mcols
#' @importFrom IRanges IRanges
#' @importFrom dplyr left_join mutate group_by summarise cur_group_id ungroup
#' @export
summariseByExpression <- function(se, type = "exon") {
    # build the grouping index and output ranges for the requested type
    index <- switch(type,
        exon = buildExonIndex(se),
        gene = buildGeneIndex(se),
        stop("Unsupported type: '", type, "'")
    )

    # sparse aggregation matrix: rows = groups, cols = transcripts
    # multiplying by transcript counts gives group-level counts
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

    # aggregate transcript counts to group level
    counts <- aggregationMat %*% assays(se)$counts

    # for gene level, add reads that could not be assigned to any transcript
    # but were localised to a gene, stored in metadata by bambu
    if (type == "gene" && !is.null(metadata(se)$incompatibleCounts)) {
        incompat <- metadata(se)$incompatibleCounts
        if ("nonuniqueCounts" %in% names(metadata(se)))
            incompat <- incompat + metadata(se)$nonuniqueCounts
        incompat <- Matrix(
            incompat[match(rownames(counts), rownames(incompat)), ],
            sparse = TRUE
        )
        counts <- counts + incompat
    }

    # compute CPM from aggregated counts
    counts.total <- colSums(counts)
    counts.total[counts.total == 0] <- 1
    cpm <- counts / counts.total * 10^6

    # aggregate optional assays if present
    outAssays <- SimpleList(counts = counts, CPM = cpm)
    for (i in c("fullLengthCounts", "uniqueCounts")) {
        if (i %in% names(assays(se))) {
            outAssays[[i]] <- aggregationMat %*% assays(se)[[i]]
        }
    }

    SummarizedExperiment(
        assays    = outAssays,
        rowRanges = index$outRanges,
        colData   = colData(se)
    )
}

#' @param se a \code{SummarizedExperiment} object from \code{\link{bambu}}
#' @return a list with:
#'   \code{indexDf} — data.frame with columns TXNAME and group_idx
#'     (one row per exon-transcript pair);
#'   \code{outRanges} — GRanges with one entry per unique exon locus,
#'     named EX1, EX2, ..., with mcols GENEID, txNames, and exonClass.
#' @noRd
buildExonIndex <- function(se) {
    # one row per exon-transcript pair, joined with transcript metadata
    txDf <- as.data.frame(rowData(se))[, c("TXNAME", "GENEID", "novelTranscript"), drop = FALSE]
    exonRanges <- unlist(rowRanges(se), use.names = TRUE)
    flatDf <- data.frame(
        TXNAME   = names(exonRanges),
        seqnames = as.character(seqnames(exonRanges)),
        start    = start(exonRanges),
        end      = end(exonRanges),
        strand   = as.character(strand(exonRanges)),
        stringsAsFactors = FALSE
    ) %>%
        left_join(txDf, by = "TXNAME") %>%
        # assign each unique exon locus a group index
        group_by(seqnames, start, end, strand) %>%
        mutate(group_idx = cur_group_id()) %>%
        ungroup()

    # one row per unique exon locus with summarised metadata
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

    # build output GRanges with exon metadata
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

#' @param se a \code{SummarizedExperiment} object from \code{\link{bambu}}
#' @return a list with:
#'   \code{indexDf} — data.frame with columns TXNAME and group_idx
#'     (one row per transcript);
#'   \code{outRanges} — GRangesList with one element per gene, named by
#'     GENEID, with mcols GENEID and newGeneClass.
#' @noRd
buildGeneIndex <- function(se) {
    # one row per transcript; group_idx preserves first-appearance order of GENEIDs
    rd <- as.data.frame(rowData(se))[, c("TXNAME", "GENEID", "txClassDescription"), drop = FALSE] %>%
        mutate(group_idx = match(GENEID, unique(GENEID)))

    # one row per gene with class label derived from txClassDescription
    groupMeta <- rd %>%
        group_by(group_idx) %>%
        summarise(
            GENEID       = GENEID[1],
            newGeneClass = ifelse(grepl("ENSG", GENEID[1]), "annotation", unique(txClassDescription)[1]),
            .groups      = "drop"
        )

    # build output GRangesList of reduced exon ranges per gene
    outRanges <- reducedRangesByGenes(rowRanges(se))[groupMeta$GENEID]
    mcols(outRanges)$GENEID       <- groupMeta$GENEID
    mcols(outRanges)$newGeneClass <- groupMeta$newGeneClass

    list(
        indexDf   = rd[, c("TXNAME", "group_idx")],
        outRanges = outRanges
    )
}
