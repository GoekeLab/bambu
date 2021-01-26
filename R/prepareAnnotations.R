#' Function to prepare tables and genomic ranges for
#' transript reconstruction using a txdb object
#' @title prepare annotations from txdb object or gtf file
#' @param x A \code{TxDb} object or a gtf file
#' @return A \code{GRangesList} object
#' @importFrom methods is
#' @importFrom GenomicFeatures exonsBy
#' @export
#' @examples
#' gtf.file <- system.file("extdata",
#'     "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf",
#'     package = "bambu"
#' )
#' prepareAnnotations(x = gtf.file)
prepareAnnotations <- function(x) {
    if (is(x, "TxDb")) {
        exonsByTx <- exonsBy(x, by = "tx", use.names = TRUE)
        if (any(duplicated(names(exonsByTx)))) {
            warning("transcript names are not unique,
            only one transcript per ID will be kept")
            exonsByTx <- exonsByTx[!duplicated(exonsByTx)]
        }
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
        minEqClasses <- getMinimumEqClassByTx(exonsByTx)
        mcols(exonsByTx)$eqClass <- minEqClasses$eqClass[match(
            names(exonsByTx),
            minEqClasses$queryTxId
        )]
        return(exonsByTx)
    } else {
        if (grepl(".gtf", x)) {
            return(prepareAnnotationsFromGTF(x))
        }
    }
}
