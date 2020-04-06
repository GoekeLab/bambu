#' reconstruct spliced transripts
#' @title CONSTRUCTSPLICEDREADCLASSTABLES
#' @param uniqueJunctions
#' @param unlisted_junctions
#' @param readGrglist
#' @param readNames
#' @importFrom unstrsplit getFromNamespace
constructSplicedReadClassTables <- function(uniqueJunctions, unlisted_junctions, readGrglist, readNames, quickMode = FALSE, verbose = FALSE){
  options(scipen = 999)

  allJunctionToUniqueJunctionOverlap <- findOverlaps(unlisted_junctions,
                                                     uniqueJunctions,type = 'equal',
                                                     ignore.strand = TRUE)

  junctionsByReadListCorrected <- splitAsList(uniqueJunctions$mergedHighConfJunctionId[subjectHits(allJunctionToUniqueJunctionOverlap)],names(unlisted_junctions))

  intronStartTMP <- start(uniqueJunctions[uniqueJunctions$mergedHighConfJunctionIdAll_noNA[subjectHits(allJunctionToUniqueJunctionOverlap)]])
  intronEndTMP <- end(uniqueJunctions[uniqueJunctions$mergedHighConfJunctionIdAll_noNA[subjectHits(allJunctionToUniqueJunctionOverlap)]])
  exon_0size <- which(intronStartTMP[-1] <= intronEndTMP[-length(intronEndTMP)] & names(unlisted_junctions[-1]) == names(unlisted_junctions[-length(unlisted_junctions)]))
  if(length(exon_0size) > 0) {
    intronStartTMP[-1][exon_0size] <- intronEndTMP[-length(intronEndTMP)][exon_0size]+1
  }
  intronStartCoordinates <- splitAsList(intronStartTMP, names(unlisted_junctions))
  intronEndCoordinates <- splitAsList(intronEndTMP, names(unlisted_junctions))

  #annotated strand of junctions for each read based on the infered read strand
  ## read strand assignment, is this necessary/helpful? Check again, add evaluation step?
  unlisted_junctions_strand <- uniqueJunctions$strand.mergedHighConfJunction[subjectHits(allJunctionToUniqueJunctionOverlap)]
  unlisted_junctions_strandList <- splitAsList(unlisted_junctions_strand, names(unlisted_junctions))
  strandJunctionSum <- sum(unlisted_junctions_strandList == '-') - sum(unlisted_junctions_strandList == '+')

  uniqueReadNames <- unique(names(unlisted_junctions))
  readStrand <- rep('*', length(uniqueReadNames))
  names(readStrand) <- uniqueReadNames
  readStrand[names(strandJunctionSum)][strandJunctionSum<0] <- '+'
  readStrand[names(strandJunctionSum)][strandJunctionSum>0] <- '-'
  strand(unlisted_junctions) <- readStrand[names(unlisted_junctions)]

  readTable <- tbl_df(data.frame(matrix(ncol = 12, nrow = length(uniqueReadNames))))
  colnames(readTable) <- c('chr', 'start', 'end', 'strand', 'blockCount', 'exonStarts', 'exonEnds', 'intronEnds', 'intronStarts', 'readClassId', 'readId', 'confidenceType')

  #uniqueReadNames <- names(which(highConfReadSet))
  readTable[, 'chr']    <-  as.character(unique(seqnames(readGrglist[uniqueReadNames])))  # as.character(unique(seqnames(readGrglist)))
  readTable[, 'start'] <- pmin(min(start(readGrglist[uniqueReadNames])),min(intronStartCoordinates[uniqueReadNames] -2))  # min(start(readGrglist))
  readTable[, 'end']   <- pmax(max(end(readGrglist[uniqueReadNames])),max(intronEndCoordinates[uniqueReadNames]+2))  # max(end(readGrglist))
  readTable[, 'intronEnds'] <- sapply(intronEndCoordinates[uniqueReadNames],paste0, collapse=',') #### replace with unstrsplit function (should be faster)
  readTable[, 'intronStarts'] <- sapply(intronStartCoordinates[uniqueReadNames],paste0, collapse=',')
  readTable[, 'strand'] <- readStrand[uniqueReadNames]
  readTable[, 'confidenceType'] <- 'highConfidenceJunctionReads'
  readTable[, 'readClassId'] <- NA


  #include reads with low confidence junctions, need to be assigned to real transcripts
  readTable[sum(is.na(junctionsByReadListCorrected[uniqueReadNames]))>0, 'confidenceType'] <- 'lowConfidenceJunctionReads'

  ## currently the 90% quantile is used to define the start and end ## this is the slowest part in the function.
  if(quickMode==FALSE){
    readClassTable <- readTable %>%
      group_by(chr, strand, intronEnds, intronStarts) %>%
      mutate(start = quantile(start, 0.1), end = quantile(end, 0.9), readCount = n()) %>%
      distinct(chr, strand, intronEnds, intronStarts, .keep_all = TRUE) %>%
      ungroup()
  } else {
    readClassTable <- readTable %>%
      group_by(chr, strand, intronEnds, intronStarts) %>%
      mutate(start = min(start), end = max(end), readCount = n()) %>%
      distinct(chr, strand, intronEnds, intronStarts, .keep_all = TRUE) %>%
      ungroup()
  }

  readClassTable[,'readClassId'] <- paste('rc', 1:nrow(readClassTable), sep = '.')
  gc()

  exonEndsShifted <- paste(readClassTable$intronStarts, as.integer(readClassTable$end) + 1, sep = ',')
  exonStartsShifted <- paste(as.integer(readClassTable$start) - 1, readClassTable$intronEnds, sep = ',')

  exonsByReadClass <- makeGRangesListFromFeatureFragments(seqnames = readClassTable$chr,
                                                          fragmentStarts = exonStartsShifted,
                                                          fragmentEnds = exonEndsShifted,
                                                          strand = readClassTable$strand)
  exonsByReadClass <- narrow(exonsByReadClass, start = 2, end = -2)  # correct junction to exon differences in coordinates
  names(exonsByReadClass) <- readClassTable$readClassId

  #add exon rank and exon_endRank

  ## combine new transcripts with annotated transcripts based on identical intron pattern
  unlistData <- unlist(exonsByReadClass, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(elementNROWS(exonsByReadClass)), names = NULL)

  exon_rank <- sapply(width((partitioning)), seq, from = 1)
  exon_rank[which(readClassTable$strand == '-')] <- lapply(exon_rank[which(readClassTable$strand == '-')], rev)  # * assumes positive for exon ranking
  exon_endRank <- lapply(exon_rank, rev)
  unlistData$exon_rank <- unlist(exon_rank)
  unlistData$exon_endRank <- unlist(exon_endRank)

  exonsByReadClass <- relist(unlistData, partitioning)

  readClassTable <- readClassTable %>%
    dplyr::select(chr.rc = chr, strand.rc = strand, intronStarts, intronEnds, confidenceType, readCount)

  mcols(exonsByReadClass) <- readClassTable
  return(exonsByReadClass)
}

######## FROM HERE STILL HAVE TO GO THROUGH ##########
#'@title CONSTRUCTUNSPLICEDREADCLASSES
#'@description reconstruct read classes using unspliced reads that fall within exons from annotations
#'@param granges
#'@param grangesReference
#'@param readNames
#'@param confidenceType
#'@param prefix
#'@param stranded
constructUnsplicedReadClasses <- function(granges, grangesReference,
                                          readNames, confidenceType='unspliced',
                                          prefix='unspliced', stranded=TRUE) {
  #unlistedExons <- unlist(exonsByTx)
  #uniqueExons <- unique(granges(unlistedExons))
  #singleExonReads <- unlist(granges)

  hitsWithin <- findOverlaps(granges,grangesReference, ignore.strand=!stranded, type='within', select='all')  # find reds
  hitsDF <- tbl_df(hitsWithin)
  hitsDF$chr <- as.character(seqnames(grangesReference)[subjectHits(hitsWithin)])
  hitsDF$start <- start(grangesReference)[subjectHits(hitsWithin)]
  hitsDF$end <- end(grangesReference)[subjectHits(hitsWithin)]
  if(stranded==FALSE) {
    hitsDF$strand = '*'
  } else {
    hitsDF$strand <- as.character(strand(grangesReference)[subjectHits(hitsWithin)])
  }
  ## create single exon read class by using the minimum end and maximum start of all overlapping exons (identical to minimum equivalent class)
  hitsDFGrouped <- hitsDF %>% group_by(queryHits) %>% mutate(maxStart=max(start), minEnd = min(end)) %>% dplyr::select(queryHits, chr, maxStart, minEnd, strand) %>% distinct() %>% group_by(chr, maxStart, minEnd, strand) %>% mutate(readClassId = paste0('rc',prefix,'.', group_indices())) %>% ungroup()

  readClassTableUnspliced <- hitsDFGrouped %>%
    dplyr::select(chr, start=maxStart, end=minEnd, strand, readClassId) %>%
    group_by(readClassId) %>%
    mutate(readCount=n()) %>%
    distinct() %>%
    ungroup() %>%
    mutate(confidenceType=confidenceType,
           intronStarts=NA,
           intronEnds=NA) %>%
    dplyr::select(chr, start, end, strand,intronStarts,intronEnds, confidenceType, readClassId, readCount)

  readTableUnspliced <-  dplyr::select(hitsDFGrouped, readClassId) %>%
    mutate(confidenceType=confidenceType, strand=as.character(strand(granges[hitsDFGrouped$queryHits])), readId=readNames[as.integer(names(granges[hitsDFGrouped$queryHits]))]) %>%
    dplyr::select(readId,  readClassId, confidenceType, strand)

  exByReadClassUnspliced <- GRanges(seqnames=readClassTableUnspliced$chr, ranges=IRanges(start=readClassTableUnspliced$start, end=readClassTableUnspliced$end), strand=readClassTableUnspliced$strand)


  #exByReadClassUnspliced$exon_id <- paste0('exId',prefix,'.',1:length(exByReadClassUnspliced))
  #exByReadClassUnspliced$exon_name <- paste0('ex',prefix,'.', 1:length(exByReadClassUnspliced))
  exByReadClassUnspliced$exon_rank <- 1
  exByReadClassUnspliced$exon_endRank <- 1
  partitioning <- PartitioningByEnd(1:length(exByReadClassUnspliced))
  exByReadClassUnspliced <- relist(exByReadClassUnspliced,partitioning)
  names(exByReadClassUnspliced) <- readClassTableUnspliced$readClassId

  readClassTableUnspliced <- readClassTableUnspliced %>% dplyr::select(chr.rc = chr, strand.rc = strand, intronStarts, intronEnds, confidenceType, readCount)#, readClassId, readCount)

  mcols(exByReadClassUnspliced) <- readClassTableUnspliced
  return(list(exonsByReadClass = exByReadClassUnspliced, readClassTable = readClassTableUnspliced, readTable = readTableUnspliced))
}

##### UNTIL HERE #####

#'@title COMBINEWITHANNOTATIONS
#'@description function to add annotations to reconstructed transcripts, and replace identical transcripts with reference
#'@param txList
#'@param txdbTablesList
#'@param matchOnly
combineWithAnnotations <- function(txList, txdbTablesList, matchOnly = T) {
  txExonsByCut <- cutStartEndFromGrangesList(txList$exonsByTx)
  annotationsExonsByCut <- cutStartEndFromGrangesList(txdbTablesList$exonsByTx)
  spliceOverlaps <- findSpliceOverlapsQuick(txExonsByCut,annotationsExonsByCut)

  annotatedTxId <- rep(NA, nrow(txList$txTable))

  annotatedTxId[unique(queryHits(spliceOverlaps[mcols(spliceOverlaps)$equal==TRUE]))] <- txdbTablesList$txIdToGeneIdTable$referenceTXNAME[subjectHits(spliceOverlaps[mcols(spliceOverlaps)$equal==TRUE])[!duplicated(queryHits(spliceOverlaps[mcols(spliceOverlaps)$equal==TRUE]))]]

  txList$txTable$isAnnotated <- !is.na(annotatedTxId)

  combinedTxId <- txList$txTable$txId
  combinedTxId[txList$txTable$isAnnotated] <- annotatedTxId[txList$txTable$isAnnotated]
  # readTxTable include reference ids
  txList$readTxTable$txId <- combinedTxId[match(txList$readTxTable$txId, txList$txTable$txId)]
  #tableTx include reference Ids
  txList$txTable$txId <- combinedTxId
  #replace reconstructed transcripts with annotated transcripts when equal splice sites are used
  txList$exonsByTx[txList$txTable$isAnnotated] <- txdbTablesList$exonsByTx[annotatedTxId[txList$txTable$isAnnotated]]
  names(txList$exonsByTx) <- combinedTxId

  #inlcude gene id
  #(1) based on intron match
  unlistedIntrons <- unlist(myGaps(txList$exonsByTx))
  overlapsNewIntronsAnnotatedIntrons <- findOverlaps(unlistedIntrons,txdbTablesList[['unlisted_introns']],type='equal',select='all', ignore.strand=FALSE)

  maxGeneCountPerNewTx <- tbl_df(data.frame(txId=names(unlistedIntrons)[queryHits(overlapsNewIntronsAnnotatedIntrons)],geneId=txdbTablesList[['unlisted_introns']]$geneId[subjectHits(overlapsNewIntronsAnnotatedIntrons)], stringsAsFactors=FALSE)) %>% group_by(txId, geneId) %>% mutate(geneCount = n()) %>% distinct() %>% group_by(txId) %>% filter(geneCount==max(geneCount)) %>% filter(!duplicated(txId)) %>% ungroup()

  geneIdByIntron <- rep(NA,nrow(txList$txTable))
  geneIdByIntron <- maxGeneCountPerNewTx$geneId[match(txList$txTable$txId, maxGeneCountPerNewTx$txId)]

  #(2) based on exon match

  exonMatchGene <- findOverlaps(txList$exonsByTx,txdbTablesList[['exonsByTx']],select = 'arbitrary',minoverlap = 20)
  geneIdByExon <- rep(NA,nrow(txList$txTable))
  geneIdByExon[!is.na(exonMatchGene)] <- txdbTablesList[['txIdToGeneIdTable']][exonMatchGene[!is.na(exonMatchGene)],'GENEID']
  geneIdByExon[!is.na(geneIdByIntron)] <-  geneIdByIntron[!is.na(geneIdByIntron)]

  exonMatchGene <- findOverlaps(txList$exonsByTx[is.na(geneIdByExon)],txList$exonsByTx[!is.na(geneIdByExon)],select = 'arbitrary',minoverlap = 20)
  while(any(!is.na(exonMatchGene))) {
    show('annoted new tx with existing gene id based on overlap with intermediate new tx')
    geneIdByExon[is.na(geneIdByExon)][!is.na(exonMatchGene)] <- geneIdByExon[!is.na(geneIdByExon)][exonMatchGene[!is.na(exonMatchGene)]]
    exonMatchGene <- findOverlaps(txList$exonsByTx[is.na(geneIdByExon)],txList$exonsByTx[!is.na(geneIdByExon)],select = 'arbitrary',minoverlap = 20)
  }

  geneLoci <- geneIdByExon ## will be used to annotate overlaping genes which do not share any exon or which are antisense

  exonSelfOverlaps <- findOverlaps(txList$exonsByTx[is.na(geneIdByExon)],txList$exonsByTx[is.na(geneIdByExon)],select = 'all',minoverlap = 20)
  hitObject = tbl_df(exonSelfOverlaps) %>% arrange(queryHits, subjectHits)
  length_tmp = 1
  while(nrow(hitObject)>length_tmp) {
    show('annotated transcripts from unknown genes by new gene id')

    length_tmp = nrow(hitObject)
    show(length_tmp)
    hitObject= inner_join(hitObject,hitObject,by=c("subjectHits"="queryHits")) %>% dplyr::select(queryHits,subjectHits.y) %>% distinct() %>%rename(subjectHits=subjectHits.y)
  }

  geneTxNames <- hitObject %>% group_by(queryHits) %>% mutate(geneId = paste('gene',first(subjectHits),sep='.')) %>% dplyr::select(queryHits, geneId) %>% distinct()
  geneIdByExon[is.na(geneIdByExon)] <- geneTxNames$geneId
  txList$txTable$geneId <- geneIdByExon

  #gene loci calculation
  rangeOverlap <- findOverlaps(range(txList$exonsByTx[is.na(geneLoci)]), txdbTablesList[['exonsByTx']], ignore.strand=TRUE, minoverlap = 20, select = 'arbitrary')
  geneLoci[is.na(geneLoci)][!is.na(rangeOverlap)] <- txdbTablesList[['txIdToGeneIdTable']][rangeOverlap[!is.na(rangeOverlap)],'GENEID']
  geneLoci[is.na(geneLoci)] <- geneIdByExon[is.na(geneLoci)]
  txList$txTable$geneLoci <- geneLoci

  # combine with annotated transcripts which are a superset transcript of reconstructed transcripts
  txList$txTable$compatibleAnnotatedTranscripts <- countQueryHits(spliceOverlaps[mcols(spliceOverlaps)$compatible==TRUE,])
  return(txList)
}
