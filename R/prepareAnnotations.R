#' This function takes a transcript annotation (path to gtf file/txdb
#' object) and extends the metadata important for \code{bambu} for each transcript.
#' @title prepare annotations from txdb object or gtf file
#' @param x A path to gtf file or a \code{TxDb} object
#' @details For each transcript, the exons are ranked based on the 
#' strands and their positions. The metadata contains 
#' information about which gene each transcript (tx) 
#' belongs to and its transcript equivalence class. A 
#' transcript is said to be in the transcript equivalence 
#' class of tx if its set of gaps between exon junctions 
#' contains the set of gaps between exon junctions of tx.
#' @return A \code{GRangesList} object
#' @importFrom methods is
#' @importFrom GenomicFeatures exonsBy
#' @export
#' @examples
#' gtf.file <- system.file("extdata",
#'     "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf",
#'     package = "bambu"
#' )
#' gtf_annotated <- prepareAnnotations(x = gtf.file)
#' gtf_annotated
#' 
#' metadata <- mcols(gtf_annotated) 
#' metadata
prepareAnnotations <- function(x) {
    if (is(x, "TxDb")) {
        exonsByTx <- exonsBy(x, by = "tx", use.names = FALSE)
        txNames <- values(transcripts(x, columns="tx_name"))$tx_name
        if (any(duplicated(txNames))) {
            message("transcript names are not unique,
            only one transcript per ID will be kept")
            exonsByTx <- exonsByTx[!duplicated(txNames)]
            txNames <- txNames[!duplicated(txNames)]
        }
        names(exonsByTx) <- txNames
        unlistedExons <- unlist(exonsByTx, use.names = FALSE)
        partitioning <- PartitioningByEnd(cumsum(elementNROWS(exonsByTx)),
            names = NULL)
        txIdForReorder <- togroup(PartitioningByWidth(exonsByTx))
        unlistedExons <- unlistedExons[order(txIdForReorder,
            unlistedExons$exon_rank)]
        # exonsByTx is always sorted by exon rank, not by strand, make sure that
        # this is the case here
        unlistedExons$exon_endRank <- unlist(lapply(elementNROWS(exonsByTx),
            seq, to = 1), use.names = FALSE)
        unlistedExons <- unlistedExons[order(txIdForReorder,
            start(unlistedExons))]
        mcols(unlistedExons) <- mcols(unlistedExons)[, c("exon_rank",
            "exon_endRank")]
        exonsByTx <- relist(unlistedExons, partitioning)

        mcols(exonsByTx) <- suppressMessages(AnnotationDbi::select(x,
            names(exonsByTx),
            columns = c("TXNAME", "GENEID"),
            keytype = "TXNAME"))
        mcols(exonsByTx)$txid <- seq_along(exonsByTx)
        minEqClasses <- getMinimumEqClassByTx(exonsByTx)
        if(!identical(names(exonsByTx),minEqClasses$queryTxId)) warning('eq classes might be incorrect')
        mcols(exonsByTx)$eqClass <- minEqClasses$eqClass
        mcols(exonsByTx)$eqClassById <- minEqClasses$eqClassById
    } else {
        tryCatch({
            exonsByTx = prepareAnnotationsFromGTF(x)
            },
        error = function(cond){
            stop("Input annotation file not readable. ",
                "Requires .gtf/.gff format or TxDb object")
        })
    }
    return(exonsByTx)
}
