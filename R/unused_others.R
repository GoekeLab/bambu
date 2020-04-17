###########ALL FUNCTIONS HERE ARE NON ESSENTIAL #############

#'@title PLOTGRANGESLISTBYNAMES
#'@param grangesList
plotGRangesListByNames<-function(grangesList)
{

  plotTrack = unlist(grangesList)
  plotTrack$grouping <- names(plotTrack)
  plotTrack=AnnotationTrack(plotTrack,group=plotTrack$grouping, id=plotTrack$grouping,fill=2,col.line=1,fontcolor.group=1,fontcolor.item=1)
  plotTracks(list(plotTrack),showId=TRUE)
}

#' @describeIn plotGRangesListByNames
makeTrackFromGrangesList <- function(grangesList)
{

  plotTrack = unlist(grangesList)
  plotTrack$grouping <- names(plotTrack)
  plotTrack=AnnotationTrack(plotTrack,group=plotTrack$grouping, id=plotTrack$grouping,fill=2,col.line=1,fontcolor.group=1,fontcolor.item=1)
  return(plotTrack)
}


#' @title GRANGESLISTTOBED
#' @param x
#' @noRd
grangesListToBed<-function(x) # note: 0 based coordinates
{
  # useful: generate UCSC custom track with the following steps


  xUnlist <- unlist(range(x))

  bedData = data.frame(as.character(seqnames(xUnlist)),start(xUnlist)-1,end(xUnlist),names(xUnlist),rep(1000,length(xUnlist)),as.character(strand(xUnlist)),start(xUnlist)-1,end(xUnlist),rep(0,length(xUnlist)),elementNROWS(x),unlist(lapply(width(x),paste,collapse=',')),unlist(lapply((start(x)-min(start(x))),paste,collapse=','))) #, stringsAsFactors=FALSE
  colnames(bedData) <-  c('chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts')
  return(bedData)
}

#'@title findOverlapsOverhang
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



######### ABOVE FUNCTIONS MIGHT BE USEFUL BUT ARE NOT USED IN THIS PACKAGE ######

##### TODO: Check if this might be useful to filter reads from new transcripts?
classifyReadClasses <- function(readClassList) {

  exByTx_singleBpStartEnd <- cutStartEndFromGrangesList(readClassListFull$exonsByReadClass)
  spliceOverlaps=findSpliceOverlapsQuick(exByTx_singleBpStartEnd,exByTx_singleBpStartEnd)
  spliceOverlapsSelected =spliceOverlaps[mcols(spliceOverlaps)$compatible==TRUE,]
  txIsSubsetOf <- countQueryHits(spliceOverlapsSelected)-1
  txHasSubsetIn <- countSubjectHits(spliceOverlapsSelected)-1

  readClassListFull$readClassTable$txIsSubsetOf <- txIsSubsetOf


  setIncompatible <- (readClassListFull$readClassTable$readClassId %in% distTable$readClassId[distTable$compatible==FALSE])
  spliceOverlapsUnexplained=findSpliceOverlapsQuick(exByTx_singleBpStartEnd[setIncompatible],exByTx_singleBpStartEnd[setIncompatible])
  spliceOverlapsSelected =spliceOverlapsUnexplained[mcols(spliceOverlapsUnexplained)$compatible==TRUE,]

  #     txIsSubsetOfUnexplained <- rep(0, length(exByTx_singleBpStartEnd))
  #     txIsSubsetOfUnexplained[setIncompatible] <- countQueryHits(spliceOverlapsSelected)-1
  txHasSubsetInUnexplained <- rep(0, length(exByTx_singleBpStartEnd))
  txHasSubsetInUnexplained[setIncompatible] <- countSubjectHits(spliceOverlapsSelected)-1

  # readClassListFull$readClassTable$txIsSubsetOfUnexplained <- txIsSubsetOfUnexplained
  readClassListFull$readClassTable$txHasSubsetInUnexplained <- txHasSubsetInUnexplained


  distTableByGene <- left_join(distTable, dplyr::select(readClassListFull$readClassTable, readClassId,txIsSubsetOf, txHasSubsetInUnexplained)) %>% group_by(GENEID) %>% mutate(geneCount = sum(relativeReadCount), geneCountFractionCompatible = sum(relativeReadCount * compatible)/ geneCount)

  ## here ##
  ## here: add longest transcripts first, then validate with number of explained reads, predict tss/tes is separate task (?)
  ## return only new transcripts/genes
  ## add visualisation for new transcripts
  ##as.data.frame(filter(distTableByGene, GENEID == 'ENSG00000107104', txIsSubsetOf==0, confidenceType=='highConfidenceJunctionReads', compatible==FALSE))


}

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
  #reducedSingleExonReads$exon_id <- 1:length(reducedSingleExonReads)+max(combinedExonList$exon_id)
  #reducedSingleExonReads$exon_name <- paste0('exUnspliced.', reducedSingleExonReads$exon_id)
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
  singleExonsGeneRangesList <- relist(reducedSingleExonReads[,c('exon_rank','exon_endRank')], partitioning)
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

######### ABOVE FUNCTIONS HAVE USEFUL CODE CHUCNKS, LOOK THROUGH ######
