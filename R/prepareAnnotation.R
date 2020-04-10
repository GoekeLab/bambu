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




assignNewGeneIds <- function(exByTx, prefix='', minoverlap=5, ignore.strand=F){
  if(is.null(names(exByTx))){
    names(exByTx) <- 1:length(exByTx)
  }

  exonSelfOverlaps <- findOverlaps(exByTx,
                                   exByTx,
                                   select='all',
                                   minoverlap=minoverlap,
                                   ignore.strand=ignore.strand)
  hitObject <- tbl_df(exonSelfOverlaps) %>% arrange(queryHits, subjectHits)
  candidateList <- hitObject %>%
    group_by(queryHits) %>%
    filter(queryHits <= min(subjectHits), queryHits != subjectHits) %>%
    ungroup()

  filteredOverlapList <- hitObject %>% filter(queryHits < subjectHits)

  rm(list=c('exonSelfOverlaps','hitObject'))
  gc()
  length_tmp = 1
  while(nrow(candidateList) > length_tmp) {  # loop to include overlapping read classes which are not in order
    length_tmp <- nrow(candidateList)
    temp <- left_join(candidateList, filteredOverlapList, by=c("subjectHits"="queryHits")) %>%
      group_by(queryHits) %>%
      filter(! subjectHits.y %in% subjectHits, !is.na(subjectHits.y)) %>%
      ungroup %>%
      dplyr::select(queryHits, subjectHits.y) %>%
      distinct() %>%
      dplyr::rename(subjectHits=subjectHits.y)

    candidateList <- rbind(temp, candidateList)
    while(nrow(temp)>0) {
      ## annotated transcripts from unknown genes by new gene id
      temp= left_join(candidateList,filteredOverlapList,by=c("subjectHits"="queryHits")) %>%
        group_by(queryHits) %>%
        filter(! subjectHits.y %in% subjectHits, !is.na(subjectHits.y)) %>%
        ungroup %>%
        dplyr::select(queryHits,subjectHits.y) %>%
        distinct() %>%
        dplyr::rename(subjectHits=subjectHits.y)

      candidateList <- rbind(temp, candidateList)
    }
    ## second loop
    tst <- candidateList %>%
      group_by(subjectHits) %>%
      mutate(subjectCount = n()) %>%
      group_by(queryHits) %>%
      filter(max(subjectCount)>1) %>%
      ungroup()

    temp2 <- inner_join(tst, tst, by=c("subjectHits"="subjectHits")) %>%
      filter(queryHits.x!=queryHits.y)  %>%
      mutate(queryHits = if_else(queryHits.x > queryHits.y, queryHits.y, queryHits.x),
             subjectHits = if_else(queryHits.x > queryHits.y, queryHits.x, queryHits.y)) %>%
      dplyr::select(queryHits,subjectHits) %>%
      distinct()
    candidateList <-  distinct(rbind(temp2, candidateList))
  }

  candidateList <- candidateList %>%
    filter(! queryHits %in% subjectHits) %>%
    arrange(queryHits, subjectHits)
  idToAdd <- (which(!(1:length(exByTx) %in% unique(candidateList$subjectHits))))

  candidateList <- rbind(candidateList, tibble(queryHits=idToAdd, subjectHits=idToAdd)) %>%
    arrange(queryHits, subjectHits) %>%
    mutate(geneId = paste('gene', prefix, '.', queryHits, sep='')) %>%
    dplyr::select(subjectHits, geneId)
  candidateList$readClassId <- names(exByTx)[candidateList$subjectHits]

  candidateList <- dplyr::select(candidateList, readClassId, geneId)
  return(candidateList)
}











