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
  equal <- olap %in% olapEqual
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
      uniqueStartLengthQuery <-   sum(width(setdiff(queryStart, subjectFull)))
      uniqueStartLengthSubject <-  startDist - uniqueStartLengthQuery



      queryEndList <- rep(unlist(queryEnd),elementNROWS(subjectFull))
      myId <- rep(1:length(queryEnd),elementNROWS(subjectFull))
      byExonIntersect=pintersect(queryEndList,subjectList, resolve.empty='start.x')
      endDist <- as.integer(tapply(width(subjectList)*(width(byExonIntersect)>0)-2*width(byExonIntersect),myId,sum) + width(unlist(queryEnd)))
      uniqueEndLengthQuery <-   sum(width(setdiff(queryEnd, subjectFull)))
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

#'@title RANGEDIST
#'@param query
#'@param subject
#'@param splice
#'@param maxDist
rangesDist <- function(query, subject, splice, maxDist)
{
  qrng <- ranges(query)
  srng <- ranges(subject)
  sprng <- ranges(splice)

  setDiffQ <-   width(setdiff(qrng, srng))
  interesectS <- width(intersect(srng, sprng))
  uniqueLengthQuery <- sum(setDiffQ)
  uniqueLengthSubject <- sum(interesectS)

  queryElementsOutsideMaxDist <- sum(setDiffQ>=maxDist)
  subjectElementsOutsideMaxDist <-  sum(interesectS>=maxDist)
  compatible <- (queryElementsOutsideMaxDist == 0) & (subjectElementsOutsideMaxDist == 0)
  DataFrame(uniqueLengthQuery,uniqueLengthSubject, compatible, queryElementsOutsideMaxDist, subjectElementsOutsideMaxDist) ##### HERE ####
}

#'@title FINDSPLICEOVERLAPSQUICK
#'@description the following functions are implemented in R, I just included the within option to make them significantly faster+memorey friendly for this purpose (original code copied from https://rdrr.io/bioc/GenomicAlignments/src/R/findSpliceOverlaps-methods.R)
#'@param query
#'@param subject
#'@param ignore.strand
findSpliceOverlapsQuick <- function(query, subject, ignore.strand=FALSE) {
  olap <- findOverlaps(query, subject, ignore.strand=ignore.strand, type='within')
  olapEqual <- findOverlaps(query, subject, ignore.strand=ignore.strand, type='equal')
  if (length(olap) == 0L)
    return(.result(olap))


  query <- query[queryHits(olap)]
  subject <- subject[subjectHits(olap)]
  splice <- myGaps(query)


  compatible <- myCompatibleTranscription(query, subject, splice)
  equal <- olap %in% olapEqual
  unique <- myOneMatch(compatible, queryHits(olap))
  strandSpecific <- all(strand(query) != "*")
  mcols(olap) <- DataFrame(compatible, equal, unique, strandSpecific)
  olap
}

#'@title MYCOMPATIBLETRANSCRIPTION
#'@param query
#'@param subject
#'@param splice
myCompatibleTranscription <- function(query, subject, splice)
{
  qrng <- ranges(query)
  srng <- ranges(subject)
  sprng <- ranges(splice)

  bnds <- elementNROWS(setdiff(qrng, srng)) == 0L
  splc <- elementNROWS(intersect(srng, sprng)) == 0L
  bnds & splc
}

#'@title MYCOMPATIBLETRANSCRIPTION
#'@param query
#'@param subject
#'@param splice
myOneMatch <- function(x, idx)
{
  xcnt <- rowsum(as.integer(x), idx)[,1]
  oneMatch <- rep((xcnt == 1L), table(idx))
  unname(x & oneMatch)
}
.isNumericOrNAs <- S4Vectors:::isNumericOrNAs
myGaps <- function(x, start=NA, end=NA)
{
  if (!.isNumericOrNAs(start))
    stop("'start' must be an integer vector or NA")
  if (!is.integer(start))
    start <- as.integer(start)
  if (!.isNumericOrNAs(end))
    stop("'end' must be an integer vector or NA")
  if (!is.integer(end))
    end <- as.integer(end)

  ## seqname and strand consistent in list elements
  if (all(elementNROWS(runValue(seqnames(x))) == 1L) &&
      all(elementNROWS(runValue(strand(x))) == 1L)) {
    flat <- unlist(x, use.names=FALSE)
    gaps <- gaps(ranges(x), start, end)
    ### FIXME: this makes this function more of an 'introns' than a .gaps.
    ### FIXME: this breaks when the GRangesList is not ordered by position
    if (!is.null(mcols(x, use.names=FALSE)$query.break)) {
      insert_gaps <- as(ranges(.insertGaps(x)), "CompressedIRangesList")
      gaps <- setdiff(gaps, insert_gaps)
    }

    idx <- elementNROWS(gaps) != 0
    ## FIXME : can't handle lists with empty elements
    ##         'start' and 'end' not quite right here
    firstseg <- start(PartitioningByWidth(x))
    seqnms <- rep(seqnames(flat)[firstseg], elementNROWS(gaps))
    strand <- rep(strand(flat)[firstseg], elementNROWS(gaps))
    gr <- relist(GRanges(seqnms, unlist(gaps, use.names=FALSE), strand), gaps)
    gr
  } else {
    ### FIXME: does not handle query.break column yet
    setdiff(range(x), x)
  }

}

.extract_unoriented_intron_motif <- function(genome, junctions)
{
  mcols(junctions) <- NULL
  junctions_len <- length(junctions)
  Ldinucl_gr <- Rdinucl_gr <- junctions
  end(Ldinucl_gr) <- start(Ldinucl_gr) + 1L
  start(Rdinucl_gr) <- end(Rdinucl_gr) - 1L
  all_dinucl <- getSeq(genome, c(Ldinucl_gr, Rdinucl_gr))
  Ldinucl <- head(all_dinucl, n=junctions_len)
  Rdinucl <- tail(all_dinucl, n=junctions_len)
  xscat(Ldinucl, "-", Rdinucl)
}


spliceStrand <- function(motif){
  NATURAL_INTRON_MOTIFS_RC <- as.character(reverseComplement(DNAStringSet(NATURAL_INTRON_MOTIFS)))

  motifStrand <- ifelse(motif %in% NATURAL_INTRON_MOTIFS,'+','*')
  motifStrand[motif %in% NATURAL_INTRON_MOTIFS_RC] <- '-'
  return(motifStrand)
}



evalAnnotationOverlap <- function(intronRanges, intronsByTx, ignore.strand=FALSE)
{
  return(table(!is.na(match(intronRanges, unique(unlist(intronsByTx)),ignore.strand=ignore.strand))))
}

# function returns exon granges list where the first and last exons are trimmed to width 2, all introns are preserved, requires exon_rank and exon_endRank
trimFirstLastExons <- function(grangesListWithExonRanks) {
  unlistedGrangesList <- unlist(grangesListWithExonRanks, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(elementNROWS(grangesListWithExonRanks)), names=NULL)
  startExonsSet <- (which((unlistedGrangesList$exon_rank==1 & as.character(strand(unlistedGrangesList))!='-')|(unlistedGrangesList$exon_endRank==1 & as.character(strand(unlistedGrangesList))=='-')))
  endExonsSet <- (which((unlistedGrangesList$exon_rank==1 & as.character(strand(unlistedGrangesList))=='-')|(unlistedGrangesList$exon_endRank==1 & as.character(strand(unlistedGrangesList))!='-')))
  start(unlistedGrangesList[startExonsSet]) <-  end(unlistedGrangesList[startExonsSet])-1
  end(unlistedGrangesList[endExonsSet]) <-  start(unlistedGrangesList[endExonsSet])+1
  return(relist(unlistedGrangesList, partitioning))
}

fitBinomialModel <- function(labels.train, data.train, data.test, show.cv=TRUE, maxSize.cv=10000, ...)
{
  if(show.cv)
  {
    mySample=sample(1:length(labels.train),min(floor(length(labels.train)/2),maxSize.cv))
    data.train.cv=data.train[mySample,]#[1:floor(length(mySample)/2)],]
    labels.train.cv=labels.train[mySample]#[1:floor(length(mySample)/2)]]
    data.train.cv.test=data.train[-mySample,]#[-(1:floor(length(mySample)/2))],]
    labels.train.cv.test=labels.train[-mySample]#[-(1:floor(length(mySample)/2))]]

    cv.fit=cv.glmnet(x=data.train.cv,y=labels.train.cv,family='binomial', ...)
    predictions=predict(cv.fit,newx=data.train.cv.test,s='lambda.min')
    show('prediction accuracy (CV) (higher for splice donor than splice acceptor)')

    show( fisher.test(table(predictions>0,labels.train.cv.test)))
    show(myPerformance(labels.train.cv.test==1,predictions)$AUC	)
  }

  cv.fit=cv.glmnet(x=data.train,y=labels.train,family='binomial', ...)
  predictions= predict(cv.fit,newx=data.test,s='lambda.min')
  return(list(predictions,cv.fit))
}

## Function to predict splice site as true or false positive based on annotations, requires annotated junctions object, optional list of models learned on the data
predictSpliceJunctions <- function(annotatedJunctions, junctionModel=NULL)
{

  annotatedJunctionsStart=unique(GRanges(seqnames=seqnames(annotatedJunctions),ranges=IRanges(start=start(annotatedJunctions),end=start(annotatedJunctions)), strand='*', mcols(annotatedJunctions)[,c('startScore','junctionStartName','annotatedStart','spliceStrand','spliceMotif')]))

  annotatedJunctionsStart$distStart.start=c(0,(start(annotatedJunctionsStart[-1])-start(annotatedJunctionsStart[-length(annotatedJunctionsStart)]))*as.integer((seqnames(annotatedJunctionsStart[-1])==seqnames(annotatedJunctionsStart[-length(annotatedJunctionsStart)]))))
  annotatedJunctionsStart$distStart.end=c((end(annotatedJunctionsStart[-length(annotatedJunctionsStart)])-end(annotatedJunctionsStart[-1]))*as.integer((seqnames(annotatedJunctionsStart[-length(annotatedJunctionsStart)])==seqnames(annotatedJunctionsStart[-1]))),0)
  annotatedJunctionsStart$annotatedStart.start = c(FALSE,annotatedJunctionsStart$annotatedStart[-length(annotatedJunctionsStart)])
  annotatedJunctionsStart$annotatedStart.end = c(annotatedJunctionsStart$annotatedStart[-1],FALSE)
  annotatedJunctionsStart$startScore.start = c(FALSE,annotatedJunctionsStart$startScore[-length(annotatedJunctionsStart)])
  annotatedJunctionsStart$startScore.end = c(annotatedJunctionsStart$startScore[-1],FALSE)
  annotatedJunctionsStart$spliceStrand.start = c(FALSE,annotatedJunctionsStart$spliceStrand[-length(annotatedJunctionsStart)])
  annotatedJunctionsStart$spliceStrand.end = c(annotatedJunctionsStart$spliceStrand[-1],FALSE)
  annotatedJunctionsStart$spliceMotif.start = c(FALSE,annotatedJunctionsStart$spliceMotif[-length(annotatedJunctionsStart)])
  annotatedJunctionsStart$spliceMotif.end = c(annotatedJunctionsStart$spliceMotif[-1],FALSE)

  annotatedJunctionsEnd=sort(unique(GRanges(seqnames=seqnames(annotatedJunctions),ranges=IRanges(start=end(annotatedJunctions),end=end(annotatedJunctions)), strand='*', mcols(annotatedJunctions)[,c('endScore','junctionEndName','annotatedEnd','spliceStrand','spliceMotif')])))

  annotatedJunctionsEnd$distEnd.start=c(0,(start(annotatedJunctionsEnd[-1])-start(annotatedJunctionsEnd[-length(annotatedJunctionsEnd)]))*as.integer((seqnames(annotatedJunctionsEnd[-1])==seqnames(annotatedJunctionsEnd[-length(annotatedJunctionsEnd)]))))
  annotatedJunctionsEnd$distEnd.end=c((end(annotatedJunctionsEnd[-length(annotatedJunctionsEnd)])-end(annotatedJunctionsEnd[-1]))*as.integer((seqnames(annotatedJunctionsEnd[-length(annotatedJunctionsEnd)])==seqnames(annotatedJunctionsEnd[-1]))),0)
  annotatedJunctionsEnd$annotatedEnd.start = c(FALSE,annotatedJunctionsEnd$annotatedEnd[-length(annotatedJunctionsEnd)])
  annotatedJunctionsEnd$annotatedEnd.end = c(annotatedJunctionsEnd$annotatedEnd[-1],FALSE)
  annotatedJunctionsEnd$endScore.start = c(FALSE,annotatedJunctionsEnd$endScore[-length(annotatedJunctionsEnd)])
  annotatedJunctionsEnd$endScore.end = c(annotatedJunctionsEnd$endScore[-1],FALSE)
  annotatedJunctionsEnd$spliceStrand.start = c(FALSE,annotatedJunctionsEnd$spliceStrand[-length(annotatedJunctionsEnd)])
  annotatedJunctionsEnd$spliceStrand.end = c(annotatedJunctionsEnd$spliceStrand[-1],FALSE)
  annotatedJunctionsEnd$spliceMotif.start = c(FALSE,annotatedJunctionsEnd$spliceMotif[-length(annotatedJunctionsEnd)])
  annotatedJunctionsEnd$spliceMotif.end = c(annotatedJunctionsEnd$spliceMotif[-1],FALSE)




  ## test start splice site given close by splice site (left/5')
  if(is.null(junctionModel))
  {
    junctionModelList <- list()
  }

  mySet.all=((annotatedJunctionsStart$distStart.start!=0)&annotatedJunctionsStart$spliceStrand!='*'&annotatedJunctionsStart$startScore>0&(annotatedJunctionsStart$distStart.start<15))
  mySet.training=(annotatedJunctionsStart$annotatedStart.start|annotatedJunctionsStart$annotatedStart)[mySet.all]

  myData=data.frame(annotatedJunctionsStart$startScore/(annotatedJunctionsStart$startScore.start+annotatedJunctionsStart$startScore),annotatedJunctionsStart$startScore,annotatedJunctionsStart$distStart.start,(annotatedJunctionsStart$spliceStrand.start=='+'),annotatedJunctionsStart$spliceStrand.start=='-',(annotatedJunctionsStart$spliceStrand=='+'))[mySet.all,]#,(annotatedJunctionsStart$spliceMotif.start),(annotatedJunctionsStart$spliceMotif))[mySet.training,]
  colnames(myData) <- paste('A',1:ncol(myData),sep='.')

  modelmatrix=model.matrix(~A.1+A.2+A.3+A.4+A.5, data=data.frame((myData))) #+A.6+A.7+A.8+A.9+A.10

  if(is.null(junctionModel))
  {


    myResults= fitBinomialModel(labels.train=as.integer(annotatedJunctionsStart$annotatedStart)[mySet.all][mySet.training], data.train=modelmatrix[mySet.training,], data.test=modelmatrix, show.cv=TRUE, maxSize.cv=10000)

    junctionModelList[['spliceSitePredictionStart.start']] <- myResults[[2]]
    predictions = myResults[[1]]
  }
  else{

    predictions= predict(junctionModel[['spliceSitePredictionStart.start']],newx=modelmatrix,s='lambda.min')
  }
  spliceSitePredictionStart.start <- rep(NA, length(annotatedJunctionsStart))
  names(spliceSitePredictionStart.start) <- annotatedJunctionsStart$junctionStartName
  spliceSitePredictionStart.start[mySet.all]=predictions
  annotatedJunctions$spliceSitePredictionStart.start <- rep(NA, length(annotatedJunctions))
  annotatedJunctions$spliceSitePredictionStart.start <- spliceSitePredictionStart.start[annotatedJunctions$junctionStartName]


  ## test start splice site given close by splice site (right/3')

  mySet.all=((annotatedJunctionsStart$distStart.end!=0)&annotatedJunctionsStart$spliceStrand!='*'&annotatedJunctionsStart$startScore>0&abs(annotatedJunctionsStart$distStart.end)<15)
  mySet.training=(annotatedJunctionsStart$annotatedStart.end|annotatedJunctionsStart$annotatedStart)[mySet.all]


  myData=data.frame(annotatedJunctionsStart$startScore/(annotatedJunctionsStart$startScore.end+annotatedJunctionsStart$startScore),annotatedJunctionsStart$startScore,annotatedJunctionsStart$distStart.end,(annotatedJunctionsStart$spliceStrand.end=='+'),(annotatedJunctionsStart$spliceStrand.end=='-'),(annotatedJunctionsStart$spliceStrand=='+'))[mySet.all,]#,(annotatedJunctionsStart$spliceMotif.start),(annotatedJunctionsStart$spliceMotif))[mySet.training,]
  colnames(myData) <- paste('A',1:ncol(myData),sep='.')

  modelmatrix=model.matrix(~A.1+A.2+A.3+A.4+A.5, data=data.frame((myData))) #+A.6+A.7+A.8+A.9+A.10
  if(is.null(junctionModel))
  {


    myResults= fitBinomialModel(labels.train=as.integer(annotatedJunctionsStart$annotatedStart)[mySet.all][mySet.training], data.train=modelmatrix[mySet.training,], data.test=modelmatrix, show.cv=TRUE, maxSize.cv=10000)

    junctionModelList[['spliceSitePredictionStart.end']] <- myResults[[2]]
    predictions =myResults[[1]]
  }
  else{

    predictions= predict(junctionModel[['spliceSitePredictionStart.end']],newx=modelmatrix,s='lambda.min')
  }
  spliceSitePredictionStart.end <- rep(NA, length(annotatedJunctionsStart))
  names(spliceSitePredictionStart.end) <- annotatedJunctionsStart$junctionStartName
  spliceSitePredictionStart.end[mySet.all]=predictions
  annotatedJunctions$spliceSitePredictionStart.end <- rep(NA, length(annotatedJunctions))
  annotatedJunctions$spliceSitePredictionStart.end <- spliceSitePredictionStart.end[annotatedJunctions$junctionStartName]


  ## test end splice site given close by splice site (start/5')

  mySet.all=(annotatedJunctionsEnd$distEnd.start!=0&annotatedJunctionsEnd$spliceStrand!='*'&annotatedJunctionsEnd$endScore>0&(annotatedJunctionsEnd$distEnd.start<15))
  mySet.training=(annotatedJunctionsEnd$annotatedEnd.start|annotatedJunctionsEnd$annotatedEnd)[mySet.all]


  myData=data.frame(annotatedJunctionsEnd$endScore/(annotatedJunctionsEnd$endScore.start+annotatedJunctionsEnd$endScore),annotatedJunctionsEnd$endScore,annotatedJunctionsEnd$distEnd.start,(annotatedJunctionsEnd$spliceStrand.start=='+'),annotatedJunctionsEnd$spliceStrand.start=='-',(annotatedJunctionsEnd$spliceStrand=='+'))[mySet.all,]

  colnames(myData) <- paste('A',1:ncol(myData),sep='.')
  modelmatrix=model.matrix(~A.1+A.2+A.3+A.4+A.5, data=data.frame((myData))) #+A.6+A.7+A.8+A.9+A.10
  if(is.null(junctionModel))
  {


    myResults= fitBinomialModel(labels.train=as.integer(annotatedJunctionsEnd$annotatedEnd)[mySet.all][mySet.training], data.train=modelmatrix[mySet.training,], data.test=modelmatrix, show.cv=TRUE, maxSize.cv=10000)
    junctionModelList[['spliceSitePredictionEnd.start']] <- myResults[[2]]
    predictions =myResults[[1]]
  }
  else{

    predictions= predict(junctionModel[['spliceSitePredictionEnd.start']],newx=modelmatrix,s='lambda.min')
  }
  spliceSitePredictionEnd.start <- rep(NA, length(annotatedJunctionsEnd))
  names(spliceSitePredictionEnd.start) <- annotatedJunctionsEnd$junctionEndName
  spliceSitePredictionEnd.start[mySet.all]=predictions
  annotatedJunctions$spliceSitePredictionEnd.start <- rep(NA, length(annotatedJunctions))
  annotatedJunctions$spliceSitePredictionEnd.start <- spliceSitePredictionEnd.start[annotatedJunctions$junctionEndName]

  ## test end splice site given close by splice site (right/3')



  mySet.all=(annotatedJunctionsEnd$distEnd.end!=0&annotatedJunctionsEnd$spliceStrand!='*'&annotatedJunctionsEnd$endScore>0&abs(annotatedJunctionsEnd$distEnd.end)<15)
  mySet.training=(annotatedJunctionsEnd$annotatedEnd.end|annotatedJunctionsEnd$annotatedEnd)[mySet.all]


  myData=data.frame(annotatedJunctionsEnd$endScore/(annotatedJunctionsEnd$endScore.end+annotatedJunctionsEnd$endScore),annotatedJunctionsEnd$endScore,annotatedJunctionsEnd$distEnd.end ,annotatedJunctionsEnd$distEnd.end,(annotatedJunctionsEnd$spliceStrand.end=='+'),(annotatedJunctionsEnd$spliceStrand.end=='-'),(annotatedJunctionsEnd$spliceStrand=='+'))[mySet.all,]

  colnames(myData) <- paste('A',1:ncol(myData),sep='.')
  modelmatrix=model.matrix(~A.1+A.2+A.3+A.4+A.5, data=data.frame((myData))) #+A.6+A.7+A.8+A.9+A.10
  if(is.null(junctionModel))
  {


    myResults= fitBinomialModel(labels.train=as.integer(annotatedJunctionsEnd$annotatedEnd)[mySet.all][mySet.training], data.train=modelmatrix[mySet.training,], data.test=modelmatrix, show.cv=TRUE, maxSize.cv=10000)
    junctionModelList[['spliceSitePredictionEnd.end']] <- myResults[[2]]
    predictions =myResults[[1]]
  }
  else{
    predictions= predict(junctionModel[['spliceSitePredictionEnd.end']],newx=modelmatrix,s='lambda.min')
  }
  spliceSitePredictionEnd.end <- rep(NA, length(annotatedJunctionsEnd))
  names(spliceSitePredictionEnd.end) <- annotatedJunctionsEnd$junctionEndName
  spliceSitePredictionEnd.end[mySet.all]=predictions
  annotatedJunctions$spliceSitePredictionEnd.end <- rep(NA, length(annotatedJunctions))
  annotatedJunctions$spliceSitePredictionEnd.end <- spliceSitePredictionEnd.end[annotatedJunctions$junctionEndName]
  if(is.null(junctionModel))
  { junctionModel=junctionModelList}
  return(list(annotatedJunctions, junctionModel))
}


