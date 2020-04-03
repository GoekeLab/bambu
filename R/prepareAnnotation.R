#' Function to prepare tables and genomic ranges for transript reconstruction using a txdb object
#' @title PREPAREANNOTATIONS
#' @param txdb
#' @export
#' @examples
#' \dontrun{
#'  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#'  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#'  prepareAnnotations(txdb)
#'  }
prepareAnnotations <- function(txdb) {
  exonsByTx = exonsBy(txdb,by='tx', use.names=TRUE)
  unlistedExons <- unlist(exonsByTx, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(elementNROWS(exonsByTx)), names=NULL)
  txIdForReorder <- togroup(PartitioningByWidth(exonsByTx))
  unlistedExons <- unlistedExons[order(txIdForReorder, unlistedExons$exon_rank)]  #'exonsByTx' is always sorted by exon rank, not by strand, make sure that this is the case here
  unlistedExons$exon_endRank <- unlist(sapply(elementNROWS(exonsByTx),seq,to=1))
  unlistedExons <- unlistedExons[order(txIdForReorder, start(unlistedExons))]
  mcols(unlistedExons) <- mcols(unlistedExons)[,c('exon_rank','exon_endRank')]
  exonsByTx <- relist(unlistedExons, partitioning)

  mcols(exonsByTx) <-  AnnotationDbi::select(txdb, names(exonsByTx),
                                             columns=c("TXNAME", "GENEID"),
                                             keytype="TXNAME")
  minEqClasses <- getMinimumEqClassByTx(exonsByTx)
  mcols(exonsByTx)$eqClass <- minEqClasses$eqClass[match(names(exonsByTx),minEqClasses$queryTxId)]
  return(exonsByTx)
}


getMinimumEqClassByTx <- function(exonsByTranscripts) {

  exByTxAnnotated_singleBpStartEnd <- cutStartEndFromGrangesList(exonsByTranscripts)  # estimate overlap only based on junctions
  spliceOverlaps=findSpliceOverlapsQuick(exByTxAnnotated_singleBpStartEnd,exByTxAnnotated_singleBpStartEnd)  ## identify transcripts which are compatbile with other transcripts (subsets by splice sites)
  spliceOverlapsSelected =spliceOverlaps[mcols(spliceOverlaps)$compatible==TRUE,] ## select splicing compatible transcript matches

  minReadClassTable <- as_tibble(spliceOverlapsSelected) %>%
    dplyr::select(queryHits, subjectHits)
  minReadClassTable$queryTxId <- names(exByTxAnnotated_singleBpStartEnd)[minReadClassTable$queryHits]
  minReadClassTable$subjectTxId <- names(exByTxAnnotated_singleBpStartEnd)[minReadClassTable$subjectHits]
  minReadClassTable <- minReadClassTable %>%
    group_by(queryTxId) %>%
    arrange(queryTxId, subjectTxId) %>%
    mutate(eqClass = paste(subjectTxId, collapse='.'), minEqClassSize = n()) %>%
    dplyr::select(queryTxId, eqClass, minEqClassSize) %>%
    distinct()
  return(minReadClassTable)
}













