#' Function to prepare tables and genomic ranges for transript reconstruction using a txb object
#' @title PREPAREANNOTATIONS
#' @param txdb
#' @export
#' @examples
#' \dontrun{
#'  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#'  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#'  prepareAnnotations(txdb)
#'  }
# prepareAnnotationsAsList <- function(txdb) {
#   txdbTablesList <- list()
#   txdbTablesList[['intronsByTxEns']] <- intronsByTranscript(txdb,use.names=TRUE)
#
#
#   #txdbTablesList[['exonsByGene']]=exonsBy(txdb,by='gene')
#   txdbTablesList[['exonsByTx']]=exonsBy(txdb,by='tx', use.names=TRUE)
#
#   # add exon end rank, and reorder exonsByTx
#   unlistedExons <- unlist(txdbTablesList[['exonsByTx']], use.names = FALSE)
#   partitioning <- PartitioningByEnd(cumsum(elementNROWS(txdbTablesList[['exonsByTx']])), names=NULL)
#   txIdForReorder <- togroup(PartitioningByWidth(txdbTablesList[['exonsByTx']]))
#   unlistedExons <- unlistedExons[order(txIdForReorder, unlistedExons$exon_rank)]  # txdbTablesList[['exonsByTx']] is always sorted by exon rank, not by strand, make sure that this is the case here
#   unlistedExons$exon_endRank <- unlist(sapply(elementNROWS(txdbTablesList[['exonsByTx']]),seq,to=1))
#   unlistedExons <- unlistedExons[order(txIdForReorder, start(unlistedExons))]
#   txdbTablesList[['exonsByTx']] <- relist(unlistedExons, partitioning)
#
#
#
#   txdbTablesList[['txIdToGeneIdTable']] <- AnnotationDbi::select(txdb, names(txdbTablesList[['intronsByTxEns']]),
#                                                                  columns=c("GENEID","TXNAME", "TXCHROM", "TXSTART","TXEND","TXSTRAND","TXTYPE"),
#                                                                  keytype="TXNAME")
#   rownames(txdbTablesList[['txIdToGeneIdTable']]) <- txdbTablesList[['txIdToGeneIdTable']]$TXNAME
#
#   # add reference transcript Ids for transcripts which have identical splice patterns by choosing the longer one
#   annotationsExonsByCut <- cutStartEndFromGrangesList(txdbTablesList[['exonsByTx']])
#   spliceOverlapsAnnotations <- findSpliceOverlapsQuick(annotationsExonsByCut,annotationsExonsByCut)
#
#   identicalSplicedTx <- countQueryHits(spliceOverlapsAnnotations[mcols(spliceOverlapsAnnotations)$equal==TRUE])
#
#   multiMatchTxSubject <- names(annotationsExonsByCut)[subjectHits(spliceOverlapsAnnotations[mcols(spliceOverlapsAnnotations)$equal==TRUE])[queryHits(spliceOverlapsAnnotations[mcols(spliceOverlapsAnnotations)$equal==TRUE]) %in% which(identicalSplicedTx>1)]]
#   multiMatchTxQuery <- names(annotationsExonsByCut)[queryHits(spliceOverlapsAnnotations[mcols(spliceOverlapsAnnotations)$equal==TRUE])[queryHits(spliceOverlapsAnnotations[mcols(spliceOverlapsAnnotations)$equal==TRUE]) %in% which(identicalSplicedTx>1)]]
#   multiMatchTxSubjectWidth <- sum(width(txdbTablesList[['exonsByTx']]))[multiMatchTxSubject]
#   refTxIdTable <- tibble(sub=multiMatchTxQuery,que=multiMatchTxSubject,sum=multiMatchTxSubjectWidth) %>% group_by(sub) %>% filter(sum==max(sum))
#
#   txdbTablesList[['txIdToGeneIdTable']]$referenceTXNAME <- txdbTablesList[['txIdToGeneIdTable']]$TXNAME
#   txdbTablesList[['txIdToGeneIdTable']][refTxIdTable$sub,'referenceTXNAME'] <- refTxIdTable$que
#
#
#   ## add minimum equivalent class definitions per transcript
#   ## any transcript that is not a perfect subset of another transcript has a unique read class
#   ## for transcripts which are subsets, the minimal read class will be all transcripts that it is a subset of
#
#   spliceOverlapsCompatible =spliceOverlapsAnnotations[mcols(spliceOverlapsAnnotations)$compatible==TRUE,] ## select splicing compatible transcript matches
#
#   minEqClassTable <- as_tibble(spliceOverlapsCompatible) %>%
#     dplyr::select(queryHits, subjectHits)
#   minEqClassTable$queryTxId <- names(annotationsExonsByCut)[minEqClassTable$queryHits]
#   minEqClassTable$subjectTxId <- names(annotationsExonsByCut)[minEqClassTable$subjectHits]
#   minEqClassTable <- minEqClassTable %>%
#     group_by(queryTxId) %>%
#     arrange(queryTxId, subjectTxId) %>%
#     mutate(eqClass = paste(subjectTxId, collapse='.'), minEqClassSize = n()) %>%
#     dplyr::select(queryTxId, eqClass, minEqClassSize) %>%
#     distinct()
#
#   txdbTablesList[['txIdToGeneIdTable']]<- left_join(txdbTablesList[['txIdToGeneIdTable']],minEqClassTable, by=c('TXNAME'='queryTxId'))
#   rownames(txdbTablesList[['txIdToGeneIdTable']]) <- txdbTablesList[['txIdToGeneIdTable']]$TXNAME
#
#   show('test that rows are ordered by tx name')
#   show(identical(rownames((txdbTablesList[['txIdToGeneIdTable']])), names(txdbTablesList[['exonsByTx']])))
#
#   txdbTablesList[['unlisted_introns']] <- unlist(txdbTablesList[['intronsByTxEns']])
#   txdbTablesList[['unlisted_introns']]$txId <- names(txdbTablesList[['unlisted_introns']])
#   txdbTablesList[['unlisted_introns']]$geneId <- txdbTablesList[['txIdToGeneIdTable']][txdbTablesList[['unlisted_introns']]$txId,'GENEID']
#
#
#   return(txdbTablesList)
# }

prepareAnnotations <- function(txdb) {
  exonsByTx  <- ensembldb::exonsBy(txdb,by='tx', use.names=TRUE)
  unlistedExons <- unlist(exonsByTx, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(elementNROWS(exonsByTx)), names=NULL)
  txIdForReorder <- togroup(PartitioningByWidth(exonsByTx))
  unlistedExons <- unlistedExons[order(txIdForReorder, unlistedExons$exon_rank)]  # txdbTablesList[['exonsByTx']] is always sorted by exon rank, not by strand, make sure that this is the case here
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













