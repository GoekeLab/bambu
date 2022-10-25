#' @title prepare reference annotations for long read RNA-Seq analysis with Bambu
#' @param x A path to gtf file or a \code{TxDb} object.
#' @details This function creates a reference annotation object which is used for transcript 
#' discovery and quantification in Bambu. prepareAnnotations can use a path to a gtf file
#' or a TxDB object as input, and returns a annotation object that stores additional
#' information about transcripts which is used in Bambu. For each transcript, exons
#' are ranked from first to last exon in direction of transcription. 
#' @return A \code{GRangesList} object with additional details for each exon and transcript 
#' that are required by Bambu. Exons are ranked by the \code{exon_rank} column, corresponding 
#' to the rank in direction of transcription (from first to last exon). In addition to exon rank, 
#' gene id, transcript id, and the minimum transcript equivalent class is stored as well
#' (a transcript equivalence class of a transcript x is the collection of transcripts 
#' where their exon junctions contain, in a continuous way, the exon junctions of the transcript x). 
#' The object is designed to be used by Bambu, and the direct access of the metadata is not recommended.
#' @importFrom methods is
#' @importFrom GenomicFeatures exonsBy
#' @export
#' @seealso [bambu()] for transcript discovery and quantification from long read RNA-Seq.
#' @examples
#' gtf.file <- system.file("extdata",
#'     "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf",
#'     package = "bambu"
#' )
#' annotations <- prepareAnnotations(x = gtf.file)
#' 
#' # run bambu
#' test.bam <- system.file("extdata", 
#'     "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam", 
#'     package = "bambu")
#' fa.file <- system.file("extdata", 
#'     "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", 
#'     package = "bambu")
#' se <- bambu(reads = test.bam, annotations = annotations, 
#'     genome = fa.file, discovery = TRUE, quant = TRUE)
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
