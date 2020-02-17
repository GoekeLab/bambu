#'@title CONSTRUCTUNSPLICEDTRANSCRIPTS
#'@description function to construct unspliced transcripts which do not overlap with TSS/TES or fall in internal exons
#'@param txList
#'@param txdbTablesList
#'@param readGrglist
#'@param stranded
#'@param minReadCount
#'@param maxOverhang
constructUnsplicedTranscripts <- function(txList, txdbTablesList,readGrglist, stranded=FALSE, minReadCount = 2, maxOverhang = 5) {

  unlistedExons <- unlist(txdbTablesList[['exonsByTx']])
  combinedExonList <- c(unlist(txList$exonsByTx),unlistedExons)

  uniqueStartExons <-unique(combinedExonList[(combinedExonList$exon_rank==1 & combinedExonList$exon_endRank>1 & strand(combinedExonList)=='+') | (combinedExonList$exon_endRank==1 & combinedExonList$exon_rank>1 & strand(combinedExonList)=='-')])  ## allow for single exon genes to be included here

  uniqueEndExons <-unique(combinedExonList[(combinedExonList$exon_rank==1 & combinedExonList$exon_endRank>1 & strand(combinedExonList)=='-') | (combinedExonList$exon_endRank==1 & combinedExonList$exon_rank>1 & strand(combinedExonList)=='+')])## allow for single exon genes to be included here

  singleExonReads <- unlist(readGrglist[elementNROWS(readGrglist)==1])

  #allow for 5bp overhang at both ends
  start(singleExonReads) <-pmin(start(singleExonReads)+maxOverhang, end(singleExonReads))
  end(singleExonReads) <-pmax(end(singleExonReads)-maxOverhang, start(singleExonReads))

  #exclude reads which overlap with first exon and extend the TSS only, or which overlap with the last exon and extend the TES only, or which are within internal exons
  hitsStart <- findOverlapsOverhang(singleExonReads,uniqueStartExons, fix='end', ignore.strand=!stranded, maxgap=-1, select='arbitrary')  # find reds where the end lies within start exons
  hitsEnd <- findOverlapsOverhang(singleExonReads,uniqueEndExons, fix='start', ignore.strand=!stranded, maxgap=-1, select='arbitrary')  # find reds where the end lies within start exons
  hitsWithin <- findOverlaps(singleExonReads,unique(combinedExonList[(combinedExonList$exon_rank>1|combinedExonList$exon_endRank>1)]), ignore.strand=!stranded, type='within', select='arbitrary')  # find reds where the end lies within start exons

  singleExonReadsSelected <- singleExonReads[is.na(hitsWithin) & is.na(hitsEnd) & is.na(hitsStart)]
  if(stranded==FALSE) {
    strand(singleExonReadsSelected) <- '*'
  }
  reducedSingleExonReads <- reduce(singleExonReadsSelected, with.revmap=T)
  reducedSingleExonReads$readCount <- elementNROWS(reducedSingleExonReads$revmap)

  # filter by minimum read count
  reducedSingleExonReads <- reducedSingleExonReads[reducedSingleExonReads$readCount >=minReadCount]

  # add annotations to create exon granges list
  reducedSingleExonReads$exon_id <- 1:length(reducedSingleExonReads)+max(combinedExonList$exon_id)
  reducedSingleExonReads$exon_name <- paste0('exUnspliced.', reducedSingleExonReads$exon_id)
  reducedSingleExonReads$exon_rank <- 1
  reducedSingleExonReads$exon_endRank <- 1

  ##(1) create/Add to txTable

  txTableSingleExon <- data.frame(matrix(NA,ncol=ncol(txList$txTable), nrow=length(reducedSingleExonReads)))
  colnames(txTableSingleExon) <- colnames(txList$txTable)

  txTableSingleExon[,c('chr','start','end','strand','readCount')]<-as.data.frame(reducedSingleExonReads)[,c('seqnames','start','end','strand','readCount')]
  singleExonAnnotations <- unlistedExons[unlistedExons$exon_rank==1 & unlistedExons$exon_endRank==1]
  hitsAnnotations <- findOverlaps(reducedSingleExonReads,singleExonAnnotations, ignore.strand=!stranded, type='any', select='arbitrary')
  txTableSingleExon$annotatedTxId <- names(singleExonAnnotations)[hitsAnnotations]
  txTableSingleExon$annotatedTxStartDist <- txTableSingleExon$start - start(singleExonAnnotations)[hitsAnnotations]
  txTableSingleExon$annotatedTxEndDist<- txTableSingleExon$end - end(singleExonAnnotations)[hitsAnnotations]
  txTableSingleExon$firstJunctionCountPerRead <- NA
  txTableSingleExon$lastJunctionCountPerRead <- NA
  txTableSingleExon$firstJunctionInternalCountPerRead <- NA
  txTableSingleExon$lastJunctionInternalCountPerRead <- NA
  txTableSingleExon$confidenceType    <- 'singleExonReads'
  txTableSingleExon$combinedTxId<- paste0('txUnspliced.',1:nrow(txTableSingleExon))
  txTableSingleExon$combinedTxId[!is.na(txTableSingleExon$annotatedTxId )] <- txTableSingleExon$annotatedTxId[!is.na(txTableSingleExon$annotatedTxId )]
  txTableSingleExon$exonOverlapGeneId <-  paste0('geneUnspliced.',1:nrow(txTableSingleExon))
  txTableSingleExon$exonOverlapGeneId[!is.na(txTableSingleExon$annotatedTxId )] <- txdbTablesList$txIdToGeneIdTable[txTableSingleExon$annotatedTxId[!is.na(txTableSingleExon$annotatedTxId)],'GENEID']

  txList$txTable <- rbind(txList$txTable, txTableSingleExon)

  ##(2) Add to granges list object (exonsByTx)

  partitioning <- PartitioningByEnd(1:length(reducedSingleExonReads), names=NULL)
  singleExonsGeneRangesList <- relist(reducedSingleExonReads[,c('exon_id','exon_name','exon_rank','exon_endRank')], partitioning)
  names(singleExonsGeneRangesList) <- txTableSingleExon$combinedTxId  # avoid including names twice

  txList$exonsByTx <- c(txList$exonsByTx,singleExonsGeneRangesList)

  ##(3) add to readTxTable

  hitsReadToTx <- findOverlaps(singleExonReadsSelected,reducedSingleExonReads, ignore.strand= ! stranded, type='within', select='all')
  singleExonReadTxTable <- data.frame(readId=names(readIdToName[match(names(singleExonReadsSelected[queryHits(hitsReadToTx)]), readIdToName)]),
                                      combinedTxId=txTableSingleExon$combinedTxId[subjectHits(hitsReadToTx)],
                                      confidenceType=txTableSingleExon$confidenceType[subjectHits(hitsReadToTx)],
                                      strand=txTableSingleExon$strand[subjectHits(hitsReadToTx)])

  txList$readTxTable <- rbind(txList$readTxTable,singleExonReadTxTable)

  return(txList)
}

#'@title findOverlapOverhang
#'@description function to find ranges where the start or end is within another range (eg start/end overlaps)
#'@param query
#'@param subject
#'@param fix
#'@param select
#'@param maxgap
#'@param ignore.strand
findOverlapsOverhang <- function(query, subject, fix, select='all', maxgap = -1, ignore.strand = TRUE) {  # type=c('start','end'); maxgap=distance to start/end
  if(class(query)!='GRanges') {
    stop("query has to be GRanges")
  }
  hits <- findOverlaps(resize(query,width = 1,fix = fix, ignore.strand=T), subject, type='any', select=select, maxgap = -1, ignore.strand=ignore.strand)
  if(maxgap> (-1)) {
    hits2 <- findOverlaps(resize(query,width = 1,fix = fix, ignore.strand=T), resize(subject,width = 1,fix = fix, ignore.strand=T), type='any', select=select, maxgap=maxgap, ignore.strand=ignore.strand)
    if(select=='all') {
      hits <- sort(unique(Hits(from=c(queryHits(hits), queryHits(hits2)), to=c(subjectHits(hits), subjectHits(hits2)), nLnode = queryLength(hits),subjectLength(hits))))
    }
  }
  return(hits)
}

#'@title CUTSTARTENDFROMGRANGESLIST
#'@description function to reduce the start end end of the first and last elements in a granges list objects to a single basepair, helper to identify overlaps based on splicing only (allow for flexible TSS/TES)
#'@param grangesList
cutStartEndFromGrangesList <- function(grangesList) {
  unlistedExons <- unlist(grangesList, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(elementNROWS(grangesList)), names=NULL)
  startExonsSet <- (which((unlistedExons$exon_rank==1 & as.character(strand(unlistedExons))!='-')|(unlistedExons$exon_endRank==1 & as.character(strand(unlistedExons))=='-')))
  endExonsSet <- (which((unlistedExons$exon_rank==1 & as.character(strand(unlistedExons))=='-')|(unlistedExons$exon_endRank==1 & as.character(strand(unlistedExons))!='-')))

  start(unlistedExons[startExonsSet]) <-  end(unlistedExons[startExonsSet])-1
  end(unlistedExons[endExonsSet]) <-  start(unlistedExons[endExonsSet])+1

  return(relist(unlistedExons, partitioning))
}

#'@title CUTGRANGESLISTELEMENTS
#'@description function to reduce the start end end of all elements in a granges list objects to a single basepair
#'@param grangesList
#'@param by
cutGrangesListElements <- function(grangesList, by=5) {
  unlistedExons <- unlist(grangesList, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(elementNROWS(grangesList)), names=NULL)

  start(unlistedExons) <-  pmin(start(unlistedExons)+by,end(unlistedExons))
  end(unlistedExons) <-  pmax(start(unlistedExons),end(unlistedExons)-by)

  return(relist(unlistedExons, partitioning))
}

#'@title EXTENDGRANGESLISTELEMENTS
#'@param grangesList
#'@param by
extendGrangesListElements <- function(grangesList, by=5) {
  unlistedExons <- unlist(grangesList, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(elementNROWS(grangesList)), names=NULL)
  # unlistedExons
  start(unlistedExons) <-  pmax(0,start(unlistedExons)-by)
  end(unlistedExons) <-  end(unlistedExons)+by

  return(relist(unlistedExons, partitioning))
}

#'@title DROPGRANGESLISTELEMENTSBYWIDTH
#'@param grangesList
#'@param minWidth
#'@param cutStartEnd
dropGrangesListElementsByWidth <- function(grangesList, minWidth=5, cutStartEnd=FALSE) {
  unlistedExons <- unlist(grangesList, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(sum(width(grangesList)>=minWidth)), names=NULL)
  exonWidth <- width(unlistedExons)
  if(cutStartEnd) {
    startExonsSet <- (which((unlistedExons$exon_rank==1 & as.character(strand(unlistedExons))!='-')|(unlistedExons$exon_endRank==1 & as.character(strand(unlistedExons))=='-')))
    endExonsSet <- (which((unlistedExons$exon_rank==1 & as.character(strand(unlistedExons))=='-')|(unlistedExons$exon_endRank==1 & as.character(strand(unlistedExons))!='-')))

    start(unlistedExons[startExonsSet]) <-  end(unlistedExons[startExonsSet])-1
    end(unlistedExons[endExonsSet]) <-  start(unlistedExons[endExonsSet])+1
  }
  unlistedExons <- unlistedExons[exonWidth>= minWidth]

  return(relist(unlistedExons, partitioning))
}

#'@title SELECTSTARTEXONSFROMGRANGESLIST
#'@description function that selects the first N exons from a grangeslist object (exon_rank is required)
#'@param grangesList
#'@param exonNumber
selectStartExonsFromGrangesList <- function(grangesList, exonNumber=2) {
  unlisted_granges <- unlist(grangesList, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(pmin(elementNROWS(grangesList), exonNumber)), names=NULL)
  startExonsSet <- which(unlisted_granges$exon_rank<=exonNumber)
  return(relist(unlisted_granges[startExonsSet], partitioning))
}

#'@title SELECTENDEXONSFROMGRANGESLIST
#'@description function that selects the last N exons from a grangeslist object (exon_endRank is required)
#'@param grangesList
#'@param exonNumber
selectEndExonsFromGrangesList <- function(grangesList, exonNumber=2) {
  unlisted_granges <- unlist(grangesList, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(pmin(elementNROWS(grangesList), exonNumber)), names=NULL)
  endExonsSet <- which(unlisted_granges$exon_endRank<=exonNumber)
  return(relist(unlisted_granges[endExonsSet], partitioning))
}


### Todo: integrate into code above
#'@title FINDSPLICEOVERLAPSBYDIST
#'@description This function calcualtes compatible splice overlaps allowing for a distance threshold, and returns distance in bp between query and subject. Can be used to assign more transcripts to annotations and reads to transcripts.
#'@param query
#'@param subject
#'@param ignore.strand
#'@param maxDist
#'@param type
#'@param firstLastSeparate
#'@param dropRangesByMinLength
#'@param cutStartEnd
findSpliceOverlapsByDist <-function(query, subject, ignore.strand=FALSE, maxDist = 5, type='within', firstLastSeparate = T, dropRangesByMinLength=F, cutStartEnd = T) {

  #  with this option the first and last exons are stored and the distance for each between query and subject hits is returned
  if(firstLastSeparate) {
    queryStart <- selectStartExonsFromGrangesList(query, exonNumber = 1)
    queryEnd <- selectEndExonsFromGrangesList(query, exonNumber = 1)
    subjectStart <- selectStartExonsFromGrangesList(subject, exonNumber = 1)
    subjectEnd <- selectEndExonsFromGrangesList(subject, exonNumber = 1)
    subjectFull <- subject
  }


  #this list is used to find candidate matches. allows ranges to be dropped from query (to find matches not strictly within). more matches>slower
  if(dropRangesByMinLength) {
    queryForOverlap <- dropGrangesListElementsByWidth(query, minWidth=maxDist, cutStartEnd=cutStartEnd)
  } else if(cutStartEnd) {
    queryForOverlap <- cutStartEndFromGrangesList(query)
  } else {
    queryForOverlap <- query
  }


  # first and last exons are reduced to 2bp, internal exon and splice matches are emphasised
  query <- cutStartEndFromGrangesList(query)
  #  subject <- cutStartEndFromGrangesList(subject)

  #  queryExonsRemoved <- dropGrangesListElementsByWidth(query, minWidth=maxDist)
  subjectExtend <- extendGrangesListElements(subject, by=maxDist)

  olap <- findOverlaps(queryForOverlap, subjectExtend, ignore.strand=ignore.strand, type=type)
  olapEqual <- findOverlaps(query, cutStartEndFromGrangesList(subject), ignore.strand=ignore.strand, type='equal')

  query <- query[queryHits(olap)]
  subject <- subject[subjectHits(olap)]
  splice <- myGaps(query)


  compatible <- rangesDist(query, subject, splice, maxDist)
  equal <- (!is.na(S4Vectors::match(olap, olapEqual)))
  unique <- myOneMatch(compatible$compatible, queryHits(olap))
  strandSpecific <- all(strand(query) != "*")
  strandedMatch <- ((all(strand(query)=='-') & all(strand(subject)=='-'))| (all(strand(query)=='+') & all(strand(subject)=='+')))
  mcols(olap) <- DataFrame(compatible, equal, unique, strandSpecific, strandedMatch)

  if(firstLastSeparate) { ##################### HERE IS AN ERROR WITH THE START SEQUENCE !!! ####
    if(length(olap)>0) {
      queryStart <- ranges(queryStart[queryHits(olap)])
      subjectStart <- ranges(subjectStart[subjectHits(olap)])
      queryEnd <- ranges(queryEnd[queryHits(olap)])
      subjectEnd <- ranges(subjectEnd[subjectHits(olap)])
      subjectFull <- ranges(subjectFull[subjectHits(olap)])
      #queryFirstJunction <- ranges(selectStartExonsFromGrangesList(splice, exonNumber = 1))
      #queryLastJunction <- ranges(selectEndExonsFromGrangesList(splice, exonNumber = 1))

      startMatch <- poverlaps(unlist(queryStart), unlist(subjectStart))
      endMatch <- poverlaps(unlist(queryEnd), unlist(subjectEnd))

      #calculate distance between first and last exon matches:
      # distances corresponds to length of unique sequences in all matching exons
      # distance = width(element in query) + width(all matching elements in subject) - 2x width(intersection(query, subject)
      subjectList <- unlist(subjectFull)
      queryStartList <- rep(unlist(queryStart),elementNROWS(subjectFull))
      myId <- rep(1:length(queryStart),elementNROWS(subjectFull))
      byExonIntersect=pintersect(queryStartList,subjectList, resolve.empty='start.x')
      startDist <- as.integer(tapply(width(subjectList)*(width(byExonIntersect)>0)-2*width(byExonIntersect),myId,sum) + width(unlist(queryStart)))
      uniqueStartLengthQuery <-   sum(width(GenomicRanges::setdiff(queryStart, subjectFull)))
      uniqueStartLengthSubject <-  startDist - uniqueStartLengthQuery



      queryEndList <- rep(unlist(queryEnd),elementNROWS(subjectFull))
      myId <- rep(1:length(queryEnd),elementNROWS(subjectFull))
      byExonIntersect=pintersect(queryEndList,subjectList, resolve.empty='start.x')
      endDist <- as.integer(tapply(width(subjectList)*(width(byExonIntersect)>0)-2*width(byExonIntersect),myId,sum) + width(unlist(queryEnd)))
      uniqueEndLengthQuery <-   sum(width(GenomicRanges::setdiff(queryEnd, subjectFull)))
      uniqueEndLengthSubject <-  endDist - uniqueEndLengthQuery
    } else {
      startMatch <- NULL
      uniqueStartLengthQuery <- NULL
      uniqueStartLengthSubject <- NULL
      endMatch <- NULL
      uniqueEndLengthQuery <- NULL
      uniqueEndLengthSubject <- NULL
    }

    mcols(olap) <- DataFrame( mcols(olap), startMatch, uniqueStartLengthQuery,uniqueStartLengthSubject,endMatch, uniqueEndLengthQuery,uniqueEndLengthSubject)
  }

  olap
}
