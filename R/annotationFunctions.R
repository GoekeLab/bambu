#' Function to prepare tables and genomic ranges for transript reconstruction using a txdb object
#' @title prepare annotations from txdb object
#' @param txdb a \code{\link{TxDb}} object
#' @return A \code{\link{GrangesList}} object
#' @export
#' @examples
#' \dontrun{
#'  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#'  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#'  prepareAnnotations(txdb)
#'  }
prepareAnnotations <- function(txdb) {
  exonsByTx = exonsBy(txdb,by='tx', use.names=TRUE)
  if(any(duplicated(names(exonsByTx)))) {
    warning('transcript names are not unique, only one transcript per ID will be kept')
    exonsByTx <- exonsByTx[!duplicated(exonsByTx)]
  }

  unlistedExons <- unlist(exonsByTx, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(elementNROWS(exonsByTx)), names=NULL)
  txIdForReorder <- togroup(PartitioningByWidth(exonsByTx))
  unlistedExons <- unlistedExons[order(txIdForReorder, unlistedExons$exon_rank)]  #'exonsByTx' is always sorted by exon rank, not by strand, make sure that this is the case here
  unlistedExons$exon_endRank <- unlist(sapply(elementNROWS(exonsByTx),seq,to=1), use.names=FALSE)
  unlistedExons <- unlistedExons[order(txIdForReorder, start(unlistedExons))]
  mcols(unlistedExons) <- mcols(unlistedExons)[,c('exon_rank','exon_endRank')]
  exonsByTx <- relist(unlistedExons, partitioning)

  mcols(exonsByTx) <-  suppressMessages(AnnotationDbi::select(txdb, names(exonsByTx),
                                             columns=c("TXNAME", "GENEID"),
                                             keytype="TXNAME"))
  minEqClasses <- getMinimumEqClassByTx(exonsByTx)
  mcols(exonsByTx)$eqClass <- minEqClasses$eqClass[match(names(exonsByTx),minEqClasses$queryTxId)]
  return(exonsByTx)
}


#' Prepare annotations from gtf
#' @title prepare annotations from gtf file
#' @param gtf.file A string variable indicates the path to a gtf file.
#' @param organism as described in \code{\link{makeTxDbFromGFF}}.
#' @param dataSource as described in \code{\link{makeTxDbFromGFF}}.
#' @param taxonomyId as described in \code{\link{makeTxDbFromGFF}}.
#' @param chrominfo	as described in \code{\link{makeTxDbFromGFF}}.
#' @param miRBaseBuild as described in \code{\link{makeTxDbFromGFF}}.
#' @param metadata as described in \code{\link{makeTxDbFromGFF}}.
#' @param dbxrefTag as described in \code{\link{makeTxDbFromGFF}}.
#' @param ... see \code{\link{makeTxDbFromGFF}}.
#' @return A \code{\link{GrangesList}} object
#' @export
prepareAnnotationsFromGTF <- function(gtf.file, dataSource=NA,
                                     organism="Homo sapiens",
                                     taxonomyId=NA,
                                     chrominfo=NULL,
                                     miRBaseBuild=NA,
                                     metadata=NULL,
                                     dbxrefTag,...){
  return(prepareAnnotations(GenomicFeatures::makeTxDbFromGFF(gtf.file, format = "gtf",
                                                             organism = organism,
                                                             dataSource = dataSource,
                                                             taxonomyId = taxonomyId,
                                                             chrominfo = chrominfo,
                                                             miRBaseBuild = miRBaseBuild,
                                                             metadata = metadata,
                                                             dbxrefTag = dbxrefTag
  )))
}



#' Get minimum equivalent class by Transcript
#' @param exonsByTranscripts exonsByTranscripts
#' @noRd
getMinimumEqClassByTx <- function(exonsByTranscripts) {

  exByTxAnnotated_singleBpStartEnd <- cutStartEndFromGrangesList(exonsByTranscripts)  # estimate overlap only based on junctions
  spliceOverlaps <- findSpliceOverlapsQuick(exByTxAnnotated_singleBpStartEnd,exByTxAnnotated_singleBpStartEnd)  ## identify transcripts which are compatible with other transcripts (subsets by splice sites)
  spliceOverlapsSelected <- spliceOverlaps[mcols(spliceOverlaps)$compatible==TRUE,] ## select splicing compatible transcript matches

  minReadClassTable <- as_tibble(spliceOverlapsSelected) %>%
    dplyr::select(queryHits, subjectHits)
  minReadClassTable$queryTxId <- names(exByTxAnnotated_singleBpStartEnd)[minReadClassTable$queryHits]
  minReadClassTable$subjectTxId <- names(exByTxAnnotated_singleBpStartEnd)[minReadClassTable$subjectHits]
  minReadClassTable <- minReadClassTable %>%
    arrange(queryTxId, subjectTxId) %>%
    group_by(queryTxId) %>%
    mutate(eqClass = paste(subjectTxId, collapse='.'), minEqClassSize = n()) %>%
    dplyr::select(queryTxId, eqClass, minEqClassSize) %>%
    distinct()
  return(minReadClassTable)
}

#' Assign New Gene with Gene Ids
#' @param exByTx exByTx
#' @param prefix prefix, defaults to empty
#' @param minoverlap defaults to 5
#' @param ignore.strand defaults to FALSE
#' @noRd
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
  gc(verbose = FALSE)
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


#' Calculate distance from read class to annotation
#' @param exByTx exByTx
#' @param exByTxRef exByTxRef
#' @param maxDist defaults to 35
#' @param primarySecondaryDist defaults to 5
#' @param ignore.strand defaults to FALSE
#' @noRd
calculateDistToAnnotation <- function(exByTx, exByTxRef, maxDist = 35, primarySecondaryDist = 5, ignore.strand=FALSE) {

  ########## TODO: go through filter rules: (are these correct/up to date?)
  ## (1) select minimum distance match (note: allow for a few base pairs error?)
  ## (2) select hits within (minimum unique query sequence, no start match)
  ## (3) select hits with minimum unqiue start sequence/end sequence


  #(1)  find overlaps of read classes with annotated transcripts, allow for maxDist [b] distance for each exon; exons with size less than 35bp are dropped to find overlaps, but counted towards distance and compatibility
  spliceOverlaps <- findSpliceOverlapsByDist(exByTx,
                                             exByTxRef,
                                             maxDist=maxDist,
                                             firstLastSeparate=T,
                                             dropRangesByMinLength=T,
                                             cutStartEnd=TRUE,
                                             ignore.strand=ignore.strand)

  txToAnTable <- tbl_df(spliceOverlaps) %>%
    group_by(queryHits)  %>%
    mutate(dist = uniqueLengthQuery + uniqueLengthSubject) %>%
    mutate(txNumber = n())

  # first round of filtering should only exclude obvious mismatches
  txToAnTableFiltered <- txToAnTable %>%
    group_by(queryHits)  %>%
    arrange(queryHits, dist) %>%
    filter(dist <= (min(dist) + primarySecondaryDist)) %>%
    filter(queryElementsOutsideMaxDist + subjectElementsOutsideMaxDist == min(queryElementsOutsideMaxDist + subjectElementsOutsideMaxDist)) %>%
    filter((uniqueStartLengthQuery <= primarySecondaryDist & uniqueEndLengthQuery <= primarySecondaryDist) == max(uniqueStartLengthQuery <= primarySecondaryDist & uniqueEndLengthQuery <= primarySecondaryDist)) %>%
    mutate(txNumberFiltered = n())

  # (2) calculate splice overlap for any not in the list (all hits have a unique new exon of at least 35bp length, might be new candidates)
  setTMP <- unique(txToAnTableFiltered$queryHits)
  spliceOverlaps_rest <- findSpliceOverlapsByDist(exByTx[-setTMP],
                                                  exByTxRef,
                                                  maxDist=0,
                                                  type='any',
                                                  firstLastSeparate=T,
                                                  dropRangesByMinLength=F,
                                                  cutStartEnd=TRUE,
                                                  ignore.strand=ignore.strand)

  txToAnTableRest <- tbl_df(spliceOverlaps_rest) %>%
    group_by(queryHits) %>%
    mutate(dist=uniqueLengthQuery + uniqueLengthSubject) %>%
    mutate(txNumber=n())

  txToAnTableRest$queryHits <- (1:length(exByTx))[-setTMP][txToAnTableRest$queryHits]  # reassign IDs based on unfiltered list length

  # todo: check filters, what happens to reads with only start and end match?
  txToAnTableRest <- txToAnTableRest %>%
    group_by(queryHits)  %>%
    arrange(queryHits, dist) %>%
    filter(dist <= (min(dist) + primarySecondaryDist)) %>%
    filter(queryElementsOutsideMaxDist + subjectElementsOutsideMaxDist == min(queryElementsOutsideMaxDist + subjectElementsOutsideMaxDist)) %>%
    filter((uniqueStartLengthQuery <= primarySecondaryDist & uniqueEndLengthQuery <= primarySecondaryDist) == max(uniqueStartLengthQuery <= primarySecondaryDist & uniqueEndLengthQuery <= primarySecondaryDist)) %>%
    mutate(txNumberFiltered = n())

  # (3) find overlaps for remaining reads (reads which have start/end match, this time not cut and used to calculate distance)
  setTMPRest <- unique(c(txToAnTableRest$queryHits, setTMP))
  txToAnTableRestStartEnd <- NULL
  if(length(exByTx[-setTMPRest]) > 0) {
    spliceOverlaps_restStartEnd <- findSpliceOverlapsByDist(exByTx[-setTMPRest],
                                                            exByTxRef,
                                                            maxDist=0,
                                                            type='any',
                                                            firstLastSeparate=T,
                                                            dropRangesByMinLength=F,
                                                            cutStartEnd=F,
                                                            ignore.strand=ignore.strand)

    txToAnTableRestStartEnd <- tbl_df(spliceOverlaps_restStartEnd) %>%
      group_by(queryHits) %>%
      mutate(dist = uniqueLengthQuery + uniqueLengthSubject + uniqueStartLengthQuery + uniqueEndLengthQuery) %>%
      mutate(txNumber = n())

    txToAnTableRestStartEnd$queryHits <- (1:length(exByTx))[-setTMPRest][txToAnTableRestStartEnd$queryHits]  # reassign IDs based on unfiltered list length

    # todo: check filters, what happens to reads with only start and end match?
    txToAnTableRestStartEnd <- txToAnTableRestStartEnd %>%
      group_by(queryHits) %>%
      arrange(queryHits, dist) %>%
      filter(dist <= (min(dist) + primarySecondaryDist)) %>%
      mutate(txNumberFiltered = n())
  }

  txToAnTableFiltered <- rbind(txToAnTableFiltered, txToAnTableRest, txToAnTableRestStartEnd) %>% ungroup()

  txToAnTableFiltered$readClassId <- names(exByTx)[txToAnTableFiltered$queryHits]
  txToAnTableFiltered$annotationTxId <- names(exByTxRef)[txToAnTableFiltered$subjectHits]

  return(txToAnTableFiltered)
}

#' Get Empty Read Class From SE
#' @param se summarizedExperiment
#' @param annotationGrangesList defaults to NULL
#' @noRd
getEmptyClassFromSE <- function(se = se, annotationGrangesList = NULL){
  distTable <- data.table(metadata(se)$distTable)[,.(readClassId, annotationTxId, readCount, GENEID)]

  # filter out multiple geneIDs mapped to the same readClass based on rowData(se)
  compatibleData <- as.data.table(rowData(se), keep.rownames = TRUE)
  setnames(compatibleData, old = c("rn","geneId"),new = c("readClassId","GENEID"))
  distTable <- distTable[compatibleData[readClassId %in% unique(distTable$readClassId),.(readClassId,GENEID)], on = c("readClassId","GENEID")]

  distTable[, eqClass:=paste(sort(unique(annotationTxId)),collapse='.'), by = list(readClassId,GENEID)]

  rcTable <- unique(distTable[,.(readClassId, GENEID, eqClass, readCount)])
  rcTable[,eqClassReadCount:=sum(readCount), by = list(eqClass, GENEID)]
  rcTable <- unique(rcTable[,.(eqClass,eqClassReadCount, GENEID)])

  eqClassCountTable <- unique(distTable[,.(annotationTxId, GENEID, eqClass)][rcTable, on = c("GENEID","eqClass")])


  setnames(eqClassCountTable, c("annotationTxId"),c("TXNAME"))
  eqClassTable <- as.data.table(mcols(annotationGrangesList)[,c('GENEID','eqClass','TXNAME')])


  eqClassCountTable <- unique(merge(eqClassCountTable,eqClassTable,all = TRUE, on = c('GENEID','eqClass','TXNAME'))) # merge should be performed on both sides
  # if there exists new isoforms from eqClassCountTable, it would not be found in eqClassTable, keep them
  eqClassCountTable[is.na(eqClassReadCount), eqClassReadCount:=0]

  ## remove empty read class where there is no shared class found
  eqClassCountTable[,sum_nobs:=sum(eqClassReadCount), by = list(GENEID, TXNAME)]

  eqClassCountTable <- unique(eqClassCountTable[sum_nobs>0,.(GENEID, eqClass, eqClassReadCount,TXNAME)])

  setnames(eqClassCountTable,old = c("TXNAME","GENEID","eqClass","eqClassReadCount") , new = c("tx_id","gene_id","read_class_id","nobs"))


  return(eqClassCountTable)

}

#' From tx ranges to gene ranges
#' @noRd
txRangesToGeneRanges <- function(exByTx, TXNAMEGENEID_Map){
  # rename names to geneIDs
  names(exByTx) <- as.data.table(TXNAMEGENEID_Map)[match(names(exByTx),TXNAME)]$GENEID

  # combine gene exon ranges and reduce overlapping ones
  unlistData <- unlist(exByTx, use.names = TRUE)
  orderUnlistData <- unlistData[order(names(unlistData))]

  orderUnlistData$exon_rank <- NULL
  orderUnlistData$exon_endRank <- NULL

  exByGene <- splitAsList(orderUnlistData, names(orderUnlistData))

  exByGene <- GenomicRanges::reduce(exByGene)

  # add exon_rank and endRank
  unlistData <- unlist(exByGene, use.names = FALSE)
  partitionDesign <- cumsum(elementNROWS(exByGene))
  partitioning <- PartitioningByEnd(partitionDesign, names=NULL)
  geneStrand <- as.character(strand(unlistData))[partitionDesign]
  exon_rank <- sapply(width((partitioning)), seq, from=1)
  exon_rank[which(geneStrand == '-')] <- lapply(exon_rank[which(geneStrand == '-')], rev)  # * assumes positive for exon ranking
  exon_endRank <- lapply(exon_rank, rev)
  unlistData$exon_rank <- unlist(exon_rank)
  unlistData$exon_endRank <- unlist(exon_endRank)
  exByGene <- relist(unlistData, partitioning)

  return(exByGene)

}




