#' Function to prepare tables and genomic ranges for transript reconstruction using a txb object
#'@title PREPAREANNOTATIONS
#'@param txdb
#'@export
#'@examples
#' \dontrun{
#'  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#'  txdb <- as.list(TxDb.Hsapiens.UCSC.hg38.knownGene)
#'  prepareAnnotations(txdb)
#'  }
prepareAnnotations <- function(txdb) {
  txdbTablesList <- list()
  txdbTablesList[['intronsByTxEns']] <- intronsByTranscript(txdb,use.names=T)


  #txdbTablesList[['exonsByGene']]=exonsBy(txdb,by='gene')
  txdbTablesList[['exonsByTx']]=exonsBy(txdb,by='tx', use.names=TRUE)

  # add exon end rank, and reorder exonsByTx
  unlistedExons <- unlist(txdbTablesList[['exonsByTx']], use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(elementNROWS(txdbTablesList[['exonsByTx']])), names=NULL)
  txIdForReorder <- togroup(PartitioningByWidth(txdbTablesList[['exonsByTx']]))
  unlistedExons <- unlistedExons[order(txIdForReorder, unlistedExons$exon_rank)]  # txdbTablesList[['exonsByTx']] is always sorted by exon rank, not by strand, make sure that this is the case here
  unlistedExons$exon_endRank <- unlist(sapply(elementNROWS(txdbTablesList[['exonsByTx']]),seq,to=1))
  unlistedExons <- unlistedExons[order(txIdForReorder, start(unlistedExons))]
  txdbTablesList[['exonsByTx']] <- relist(unlistedExons, partitioning)



  txdbTablesList[['txIdToGeneIdTable']] <- AnnotationDbi::select(txdb, names(txdbTablesList[['intronsByTxEns']]),
                                                                 columns=c("GENEID","TXNAME", "TXCHROM", "TXSTART","TXEND","TXSTRAND","TXTYPE"),
                                                                 keytype="TXNAME")
  rownames(txdbTablesList[['txIdToGeneIdTable']]) <- txdbTablesList[['txIdToGeneIdTable']]$TXNAME

  # add reference transcript Ids for transcripts which have identical splice patterns by choosing the longer one
  annotationsExonsByCut <- cutStartEndFromGrangesList(txdbTablesList[['exonsByTx']])
  spliceOverlapsAnnotations <- findSpliceOverlapsQuick(annotationsExonsByCut,annotationsExonsByCut)

  identicalSplicedTx <- countQueryHits(spliceOverlapsAnnotations[mcols(spliceOverlapsAnnotations)$equal==TRUE])

  multiMatchTxSubject <- names(annotationsExonsByCut)[subjectHits(spliceOverlapsAnnotations[mcols(spliceOverlapsAnnotations)$equal==TRUE])[queryHits(spliceOverlapsAnnotations[mcols(spliceOverlapsAnnotations)$equal==TRUE]) %in% which(identicalSplicedTx>1)]]
  multiMatchTxQuery <- names(annotationsExonsByCut)[queryHits(spliceOverlapsAnnotations[mcols(spliceOverlapsAnnotations)$equal==TRUE])[queryHits(spliceOverlapsAnnotations[mcols(spliceOverlapsAnnotations)$equal==TRUE]) %in% which(identicalSplicedTx>1)]]
  multiMatchTxSubjectWidth <- sum(width(txdbTablesList[['exonsByTx']]))[multiMatchTxSubject]
  refTxIdTable <- tibble(sub=multiMatchTxQuery,que=multiMatchTxSubject,sum=multiMatchTxSubjectWidth) %>% group_by(sub) %>% filter(sum==max(sum))

  txdbTablesList[['txIdToGeneIdTable']]$referenceTXNAME <- txdbTablesList[['txIdToGeneIdTable']]$TXNAME
  txdbTablesList[['txIdToGeneIdTable']][refTxIdTable$sub,'referenceTXNAME'] <- refTxIdTable$que


  ## add minimum equivalent class definitions per transcript
  ## any transcript that is not a perfect subset of another transcript has a unique read class
  ## for transcripts which are subsets, the minimal read class will be all transcripts that it is a subset of

  spliceOverlapsCompatible =spliceOverlapsAnnotations[mcols(spliceOverlapsAnnotations)$compatible==TRUE,] ## select splicing compatible transcript matches

  minEqClassTable <- as_tibble(spliceOverlapsCompatible) %>%
    select(queryHits, subjectHits)
  minEqClassTable$queryTxId <- names(annotationsExonsByCut)[minEqClassTable$queryHits]
  minEqClassTable$subjectTxId <- names(annotationsExonsByCut)[minEqClassTable$subjectHits]
  minEqClassTable <- minEqClassTable %>%
    group_by(queryTxId) %>%
    arrange(queryTxId, subjectTxId) %>%
    mutate(eqClass = paste(subjectTxId, collapse='.'), minEqClassSize = n()) %>%
    select(queryTxId, eqClass, minEqClassSize) %>%
    distinct()

  txdbTablesList[['txIdToGeneIdTable']]<- left_join(txdbTablesList[['txIdToGeneIdTable']],minEqClassTable, by=c('TXNAME'='queryTxId'))
  rownames(txdbTablesList[['txIdToGeneIdTable']]) <- txdbTablesList[['txIdToGeneIdTable']]$TXNAME

  show('test that rows are ordered by tx name')
  show(identical(rownames((txdbTablesList[['txIdToGeneIdTable']])), names(txdbTablesList[['exonsByTx']])))

  txdbTablesList[['unlisted_introns']] <- unlist(txdbTablesList[['intronsByTxEns']])
  txdbTablesList[['unlisted_introns']]$txId <- names(txdbTablesList[['unlisted_introns']])
  txdbTablesList[['unlisted_introns']]$geneId <- txdbTablesList[['txIdToGeneIdTable']][txdbTablesList[['unlisted_introns']]$txId,'GENEID']


  return(txdbTablesList)
}


#' Function to prepare reads for processing from bam file
#'@title PREPAREDATAFROMBAM
#'@param bamFile
#'@param yieldSize
#'@export
prepareDataFromBam <- function(bamFile, yieldSize=NULL) {
  ## Todo: optimise which data is essential
  ## Todo: don't use names for reads, use index instead? Names are require a lot of memory
  if(class(bamFile)=='BamFile') {
    if(!is.null(yieldSize)) {
      yieldSize(bamFile) <- yieldSize
    } else {
      yieldSize <- yieldSize(bamFile)
    }

  } else {
    if(is.null(yieldSize)) {
      yieldSize <- NA
    }
    bamFile <- BamFile(bamFile, yieldSize = yieldSize)
  }
  bf <- open(bamFile)
  #readJunctions <- GRangesList()
  readGrglist <- GRangesList()
  # startRanges <- GRanges()
  # endRanges <- GRanges()
  # readAlignmentStrand <- NULL
  # readIdToName <- NULL
  counter <- 1

  while(isIncomplete(bf)){
    # reads=readGAlignments(bamFile, param=ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE),what = scanBamWhat()), use.names=T)
    reads=readGAlignments(bf,
                          param=ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE)), ## whether supplementary alignments would be considered as well?? To confirm with Jonathan.
                          use.names=TRUE)
    #readJunctions <- c(readJunctions,junctions(reads))
    readGrglist<-c(readGrglist,grglist(reads))
    # startRanges <- c(startRanges,GRanges(seqnames=seqnames(reads),ranges=IRanges(start=start(reads),end=start(reads))))  ## can be extracted from readGrglist, remove here? range(readGrglist[1:10])
    # endRanges <- c(endRanges,GRanges(seqnames=seqnames(reads),ranges=IRanges(start=end(reads),end=end(reads))))## can be extracted from readGrglist, remove here? range(readGrglist[1:10])
    # readIdToName <- c(readIdToName,names(reads))## replace with integer Ids or use real names
    # readAlignmentStrand <- c(readAlignmentStrand,as.character(strand(reads)))
    show(counter* yieldSize)
    counter <- counter + 1
  }

  close(bf)
  readGrglist <- readGrglist[width(readGrglist)>1]  # remove microexons Of width 1bp from list
  readNames <- names(readGrglist)
  names(readGrglist) <- 1:length(readGrglist)  # names needed to be replaced as some reads are multiple times mapped (distinct parts of the read which are compatible)

  # readJunctions <- myGaps(readGrglist)
  # names(readJunctions) <- paste('read',1:length(readJunctions),sep='.')  ## later: remove names (they blow up memory usage!), maybe hashmap? hm <- hashmap(1:nrow(data), data$names)
  # names(startRanges) <- names(readJunctions)## later: remove names (they blow up memory usage!) remove this object and replace with range(readGrglist) whenever needed
  # names(endRanges) <- names(readJunctions)## later: remove names (they blow up memory usage!)remove this object and replace with range(readGrglist) whenever needed
  # names(readIdToName) <- names(readJunctions)  ## requires a lot of memory, replace later
  # names(readAlignmentStrand) <- names(readJunctions)

  return(list(readGrglist=readGrglist, readNames=readNames)) # startRanges=startRanges, endRanges=endRanges))
}



#'@title CREATEJUNCTIONTABLE
#'@description This function creates a table of all junctions from a list grangeslist of junctions, it add strand and splice motif information
#'@param unlisted_junction_granges
#'@param genomeDB
#'@param genomeFA
#'@export
createJunctionTable <- function(unlisted_junction_granges, genomeDB=NULL, genomeFA=NULL) {
  #code below adopted from ?summarizeJunctions (credit+license)
  ##Todo: don't create junction names, instead work with indices/intergers (names are memory intensive)
  ##Todo: add code acknowledgment for summarizeJunctions

  #unlisted_junctions <- unlist(grangesList, use.names=TRUE)
  unstranded_unlisted_junctions <- unstrand(unlisted_junction_granges)
  uniqueJunctions <- sort(unique(unstranded_unlisted_junctions))
  names(uniqueJunctions) <- paste('junc',1:length(uniqueJunctions),sep='.') ##replace with integer later (or only use indices if possible)

  #calculate stranded read counts
  junctionMatchList <- as(findMatches(uniqueJunctions, unstranded_unlisted_junctions), "List")
  uniqueJunctions_score <- elementNROWS(junctionMatchList)
  junctionStrandList <- extractList(strand(unlisted_junction_granges), junctionMatchList)
  uniqueJunctions_plus_score <- sum(junctionStrandList == "+")
  uniqueJunctions_minus_score <- sum(junctionStrandList == "-")
  uniqueJunctions_mcols <- DataFrame(score=uniqueJunctions_score,
                                     plus_score=uniqueJunctions_plus_score,
                                     minus_score=uniqueJunctions_minus_score)
  ## CY/JG if we want to avoid using files, then this has to be implemented again
  # if (!is.null(genomeDB)) {
  #   if (!suppressWarnings(require(BSgenome, quietly=TRUE)))
  #     stop("you need to install the BSgenome package in order ",
  #          "to use the 'genomeDB' argument")
  #   genomeDB <- BSgenome::getBSgenome(genomeDB)
  #   seqlevelsStyle(genomeDB) <- seqlevelsStyle(uniqueJunctions)[1]
  #
  #   unoriented_intron_motif <- as.character(.extract_unoriented_intron_motif(genomeDB,
  #                                                                            uniqueJunctions))
  #
  #   uniqueJunctions_intron_strand <- spliceStrand(unoriented_intron_motif)
  #   uniqueJunctions_mcols <- cbind(uniqueJunctions_mcols,
  #                                  DataFrame(spliceMotif = unoriented_intron_motif,
  #                                            spliceStrand = uniqueJunctions_intron_strand))
  # } else
  #
  if (!is.null(genomeFA)) {
    junctionSeqStart<-getSeq(genomeFA,IRanges::shift(flank(uniqueJunctions,width=2),2)) # shift: from IRanges
    junctionSeqEnd<-getSeq(genomeFA,IRanges::shift(flank(uniqueJunctions,width=2,start=FALSE),-2)) # shift: from IRanges

    junctionMotif <- paste(junctionSeqStart,junctionSeqEnd,sep='-')
    uniqueJunctions_mcols <- cbind(uniqueJunctions_mcols,
                                   DataFrame(spliceMotif = junctionMotif,
                                             spliceStrand = spliceStrand(junctionMotif)))
  }else{
    stop('Please provide genome annotation .fa file!') ## is this necessary??
  }


  mcols(uniqueJunctions) <- uniqueJunctions_mcols
  strand(uniqueJunctions)<-uniqueJunctions$spliceStrand

  uniqueJunctions$junctionStartName <- paste(seqnames(uniqueJunctions),start(uniqueJunctions),sep=':')
  uniqueJunctions$junctionEndName <- paste(seqnames(uniqueJunctions),end(uniqueJunctions),sep=':')

  startScore <- tapply(uniqueJunctions$score,  uniqueJunctions$junctionStartName ,sum)
  uniqueJunctions$startScore <- startScore[uniqueJunctions$junctionStartName]
  endScore <- tapply(uniqueJunctions$score,  uniqueJunctions$junctionEndName ,sum)
  uniqueJunctions$endScore <- endScore[uniqueJunctions$junctionEndName]

  return(uniqueJunctions)
}

#' JUNCTIONSTRANDCORRECTION
#' @title JUNCTIONSTRANDCORRECTION
#' @param uniqueJunctions
#' @param unlisted_junction_granges
#' @param intronsByTx
#' @param stranded
junctionStrandCorrection <- function(uniqueJunctions, unlisted_junction_granges, intronsByTx, stranded) {
  ##note: the strand is not always correctly infered! based on motif introduces errors from alignment, should be based on overall read + stranded read prediction scores, first only use conservative estimates, than improve less certain junction strand predictions


  #unlisted_junction_granges <- unlist(readJunctions, use.names=TRUE)
  allJunctionToUniqueJunctionOverlap <- findOverlaps(unlisted_junction_granges, uniqueJunctions,type='equal',ignore.strand=T)

  uniqueJunctionsUpdate <- uniqueJunctions # make a copy to revert to if strand correction does not improve results

  annotatedIntronNumber <- evalAnnotationOverlap(uniqueJunctions, intronsByTx, ignore.strand=F)['TRUE']
  show('before strand correction, annotated introns:')
  show(annotatedIntronNumber)
  show(annotatedIntronNumber/length(uniqueJunctions))

  #infer strand for each read based on strand of junctions
  if(stranded==FALSE) {## todo: check that this works for stranded and unstranded data!  implemented here already but maybe check, or include evaluation
    strandStep <- TRUE
    while(strandStep) { # iterate this 2 times to predict strand more accurately using the mean junction counts

      #annotated strand of jucntions for each read based on the infered read strand
      unlisted_junction_granges_strand <- as.character(strand(uniqueJunctionsUpdate)[subjectHits(allJunctionToUniqueJunctionOverlap)])
      unlisted_junction_granges_strandList = splitAsList(unlisted_junction_granges_strand,names(unlisted_junction_granges))
      strandJunctionSum = sum(unlisted_junction_granges_strandList=='-')-sum(unlisted_junction_granges_strandList=='+')

      uniqueReadNames <- unique(names(unlisted_junction_granges))
      readStrand <- rep('*',length(uniqueReadNames))
      names(readStrand) <- uniqueReadNames
      readStrand[names(strandJunctionSum)][strandJunctionSum<0] <- '+'
      readStrand[names(strandJunctionSum)][strandJunctionSum>0] <- '-'
      strand(unlisted_junction_granges) <- readStrand[names(unlisted_junction_granges)]

      #
      #       readStrand <- rep('*',length(readJunctions))
      #       names(readStrand) <- names(readJunctions)
      #       readStrand[names(strandJunctionSum)][strandJunctionSum<0] <- '+'
      #       readStrand[names(strandJunctionSum)][strandJunctionSum>0] <- '-'
      #       strand(unlisted_junction_granges) <- readStrand[names(unlisted_junction_granges)]

      unstranded_unlisted_junction_granges <- unstrand(unlisted_junction_granges)
      junctionMatchList <- as(findMatches(uniqueJunctions, unstranded_unlisted_junction_granges), "List")
      tmp <- extractList(strand(unlisted_junction_granges), junctionMatchList)
      uniqueJunctionsUpdate$plus_score_inferedByRead <- sum(tmp == "+")
      uniqueJunctionsUpdate$minus_score_inferedByRead <- sum(tmp == "-")
      strandScoreByRead <- uniqueJunctionsUpdate$minus_score_inferedByRead - uniqueJunctionsUpdate$plus_score_inferedByRead
      strand(uniqueJunctionsUpdate[strandScoreByRead< 0 ]) = '+' ## note: here I overwrite the information from the motif which increases overlap with known junctions. is this good?? has to be evaluated
      strand(uniqueJunctionsUpdate[strandScoreByRead>0 ]) = '-'
      annotatedIntronNumberNew <- evalAnnotationOverlap(uniqueJunctionsUpdate, intronsByTx, ignore.strand=F)['TRUE']
      if(annotatedIntronNumberNew > annotatedIntronNumber & !is.na(annotatedIntronNumber)) # update junctions object if strand prediction improves overlap with annotations
      {
        show('after strand correction, annotated introns:')
        show(annotatedIntronNumberNew)
        show(annotatedIntronNumberNew/length(uniqueJunctionsUpdate))

        annotatedIntronNumber <- annotatedIntronNumberNew
        uniqueJunctions <- uniqueJunctionsUpdate
      }
      else {
        strandStep <- FALSE
      }
    }
  } else { # correct strand for stranded data using read counts
    strandStep <- TRUE
    while(strandStep) {
      strandScoreByRead <- uniqueJunctionsUpdate$minus_score - uniqueJunctionsUpdate$plus_score
      strand(uniqueJunctionsUpdate[strandScoreByRead>0 ]) = '-' ## note: here I overwrite the information from the motif which increases overlap with known junctions. is this good?? has to be evaluated
      strand(uniqueJunctionsUpdate[strandScoreByRead<(0) ]) = '+'
      annotatedIntronNumberNew <- evalAnnotationOverlap(uniqueJunctionsUpdate, intronsByTx, ignore.strand=F)['TRUE']
      if(annotatedIntronNumberNew > annotatedIntronNumber & !is.na(annotatedIntronNumber)) # update junctions object if strand prediction improves overlap with annotations
      {
        show('after strand correction, annotated introns:')
        show(annotatedIntronNumberNew)
        show(annotatedIntronNumberNew/length(uniqueJunctionsUpdate))
        annotatedIntronNumber <- annotatedIntronNumberNew
        uniqueJunctions <- uniqueJunctionsUpdate
      }
      else {
        strandStep <- FALSE
      }
    }
  }
  return(list(uniqueJunctions = uniqueJunctions, unlisted_junctions = unlisted_junction_granges))
}

#'@title FINDHIGHCONFIDENTJUNCTIONS
#'@description this function adds "mergedHighConfJunctionId" to the junciton list which contains the ID of the most likely high confident junction that each junction originates from
#'@param junctions
#'@param junctionModel
findHighConfidenceJunctions <- function(junctions, junctionModel) {
  show('reads count for all annotated junctions')
  show(sum(junctions$score[junctions$annotatedJunction]))
  show(sum(junctions$score[junctions$annotatedJunction])/ sum(junctions$score))

  ## next part: find reference junction and correction junctions to assume they are similar to their closest reference junctions
  ##note this part should be visualised (bed/bigbed track? how does correction change the reads? visualise strand prediction, all in bigbed)
  #calculate high confiddence junctions which are used as reference
  #calculation is based on distance and properties of next junctions


  candidateJunctionsPlus <- junctions[strand(junctions)=='+'|strand(junctions)=='*']
  candidateJunctionsMinus <- junctions[strand(junctions)=='-'|strand(junctions)=='*']

  highConfidentJunctionSetPlus <- candidateJunctionsPlus$score>1&candidateJunctionsPlus$spliceStrand=='+'
  highConfidentJunctionSetMinus <- candidateJunctionsMinus$score>1&candidateJunctionsMinus$spliceStrand=='-'

  if( sum(highConfidentJunctionSetPlus)>0 &  sum(highConfidentJunctionSetMinus)>0) {
    highConfJunctionsPlus <- predictSpliceJunctions(candidateJunctionsPlus[highConfidentJunctionSetPlus] ,junctionModel = junctionModel)[[1]]
    highConfJunctionsMinus <- predictSpliceJunctions(candidateJunctionsMinus[highConfidentJunctionSetMinus] ,junctionModel = junctionModel)[[1]]

    candidateJunctionsPlus$highConfJunctionPrediction = rep(FALSE, length(candidateJunctionsPlus))
    candidateJunctionsMinus$highConfJunctionPrediction = rep(FALSE, length(candidateJunctionsMinus))

    ## note: for evaluation, annotations should not be used here ## they seem to only add very little extra information
    setReferenceJunctionsPlus=    ((highConfJunctionsPlus$spliceSitePredictionStart.start>0|is.na(highConfJunctionsPlus$spliceSitePredictionStart.start)) & (highConfJunctionsPlus$spliceSitePredictionStart.end>0 | is.na(highConfJunctionsPlus$spliceSitePredictionStart.end)) & (highConfJunctionsPlus$spliceSitePredictionEnd.start>0 | is.na(highConfJunctionsPlus$spliceSitePredictionEnd.start)) & (highConfJunctionsPlus$spliceSitePredictionEnd.end>0 | is.na(highConfJunctionsPlus$spliceSitePredictionEnd.end))) | highConfJunctionsPlus$annotatedJunction

    setReferenceJunctionsMinus=   ((highConfJunctionsMinus$spliceSitePredictionStart.start>0|is.na(highConfJunctionsMinus$spliceSitePredictionStart.start)) & (highConfJunctionsMinus$spliceSitePredictionStart.end>0 | is.na(highConfJunctionsMinus$spliceSitePredictionStart.end)) & (highConfJunctionsMinus$spliceSitePredictionEnd.start>0 | is.na(highConfJunctionsMinus$spliceSitePredictionEnd.start)) & (highConfJunctionsMinus$spliceSitePredictionEnd.end>0 | is.na(highConfJunctionsMinus$spliceSitePredictionEnd.end))) |  highConfJunctionsMinus$annotatedJunction

    #
    #        highConfJunctions = highConfJunctions[setReferenceJunctions]
    candidateJunctionsPlus$highConfJunctionPrediction[highConfidentJunctionSetPlus] <- setReferenceJunctionsPlus
    candidateJunctionsMinus$highConfJunctionPrediction[highConfidentJunctionSetMinus] <- setReferenceJunctionsMinus

    mergedHighConfJunctionIdPlus <- rep(NA,length(candidateJunctionsPlus))
    mergedHighConfJunctionIdMinus <- rep(NA,length(candidateJunctionsMinus))



    ##for comparison: it seem like only 1 annotated Tx is lost with junction correction, but 210000 non annotated ones are lost (reads are counted towards annotated Tx!) systemtic evalution requried.
    # mergedHighConfJunctionId <-1:length(candJunc) ### this is for comparison: what is the difference between corrected and non corrected tx assembly?

    ##max Distance should be a parameter that can be set
    #here: assign reference junction to all junctions based on prediciton score
    for(maxDist in 0:10)
    {
      show(maxDist)
      overlapByDistPlus=findOverlaps(candidateJunctionsPlus[is.na(mergedHighConfJunctionIdPlus)], candidateJunctionsPlus[candidateJunctionsPlus$highConfJunctionPrediction], maxgap=maxDist,type='equal')
      overlapByDistMinus=findOverlaps(candidateJunctionsMinus[is.na(mergedHighConfJunctionIdMinus)], candidateJunctionsMinus[candidateJunctionsMinus$highConfJunctionPrediction], maxgap=maxDist,type='equal')

      mergedHighConfJunctionIdPlus[is.na(mergedHighConfJunctionIdPlus)][queryHits(overlapByDistPlus)] <- names(candidateJunctionsPlus)[candidateJunctionsPlus$highConfJunctionPrediction][subjectHits(overlapByDistPlus)]

      mergedHighConfJunctionIdMinus[is.na(mergedHighConfJunctionIdMinus)][queryHits(overlapByDistMinus)] <- names(candidateJunctionsMinus)[candidateJunctionsMinus$highConfJunctionPrediction][subjectHits(overlapByDistMinus)]

    }


    candidateJunctionsPlus$mergedHighConfJunctionId <- mergedHighConfJunctionIdPlus
    candidateJunctionsMinus$mergedHighConfJunctionId <- mergedHighConfJunctionIdMinus
    #in case of conflict (very rare) use reference junction with higher read count/score
    conflictJunctions <- junctions[names(candidateJunctionsMinus[!is.na(candidateJunctionsMinus$mergedHighConfJunctionId)][which(names(candidateJunctionsMinus)[!is.na(candidateJunctionsMinus$mergedHighConfJunctionId)] %in% names(candidateJunctionsPlus)[!is.na(candidateJunctionsPlus$mergedHighConfJunctionId)])])]
    scoreDiff <- junctions[candidateJunctionsPlus[names(conflictJunctions)]$mergedHighConfJunctionId]$score - junctions[candidateJunctionsMinus[names(conflictJunctions)]$mergedHighConfJunctionId]$score

    if(any(scoreDiff<0)){
      candidateJunctionsPlus[names(conflictJunctions)]$mergedHighConfJunctionId[scoreDiff <0 ] <- candidateJunctionsMinus[names(conflictJunctions)]$mergedHighConfJunctionId[scoreDiff <0 ]
    }
    if(any(scoreDiff>0)){
      candidateJunctionsMinus[names(conflictJunctions)]$mergedHighConfJunctionId[scoreDiff >0 ] <- candidateJunctionsPlus[names(conflictJunctions)]$mergedHighConfJunctionId[scoreDiff >0 ]
    }
    if(0 %in% scoreDiff){
      candidateJunctionsPlus[names(conflictJunctions)]$mergedHighConfJunctionId[scoreDiff ==0 ] <- NA
      candidateJunctionsMinus[names(conflictJunctions)]$mergedHighConfJunctionId[scoreDiff ==0 ] <- NA
    }

    #   sum(candidateJunctionsMinus$score[candidateJunctionsMinus$annotatedJunction])
    #   sumByJuncId <- tapply(candidateJunctionsMinus$score, candidateJunctionsMinus$mergedHighConfJunctionId, sum)
    #   sum(sumByJuncId[junctions[names(sumByJuncId)]$annotatedJunction])
    #   sum(candidateJunctionsPlus$score[candidateJunctionsPlus$annotatedJunction])
    #   sumByJuncId <- tapply(candidateJunctionsPlus$score, candidateJunctionsPlus$mergedHighConfJunctionId, sum)
    #   sum(sumByJuncId[junctions[names(sumByJuncId)]$annotatedJunction])

    ## from here... part above needs to be evaluated in more detail
    mergedHighConfJunctionId <- rep(NA, length(junctions))
    names(mergedHighConfJunctionId) <- names(junctions)
    mergedHighConfJunctionId[names(candidateJunctionsPlus)] <- candidateJunctionsPlus$mergedHighConfJunctionId
    mergedHighConfJunctionId[names(candidateJunctionsMinus)] <- candidateJunctionsMinus$mergedHighConfJunctionId

    junctions$mergedHighConfJunctionId <-  as.character(mergedHighConfJunctionId)
  } else {
    show('no junction correction, no high confidence reference junctions found')
    junctions$mergedHighConfJunctionId <-  names(junctions)
  }
  rm(candidateJunctionsPlus, candidateJunctionsMinus, highConfJunctionsPlus, highConfJunctionsMinus) ## clean up, more sytematically later for memory improvements
  show('reads count for all annotated junctions after correction to reference junction')
  sumByJuncId <- tapply(junctions$score, junctions$mergedHighConfJunctionId, sum)
  show(sum(sumByJuncId[junctions[names(sumByJuncId)]$annotatedJunction]))
  show(sum(sumByJuncId[junctions[names(sumByJuncId)]$annotatedJunction])/sum(junctions$score))

  #########

  return(junctions)
}




#' reconstruct spliced transripts
#' @title CONSTRUCTSPLICEDREADCLASSTABLES
#' @param uniqueJunctions
#' @param unlisted_junctions
#' @param readGrglist
#' @param readNames
constructSplicedReadClassTables <- function(uniqueJunctions, unlisted_junctions, readGrglist, readNames){ #, startRanges, endRanges) {
  options(scipen = 999)

  show('### create transcript models ###')
  allJunctionToUniqueJunctionOverlap <- findOverlaps(unlisted_junctions, uniqueJunctions,type='equal',ignore.strand=T)

  junctionsByReadListCorrected = splitAsList(uniqueJunctions$mergedHighConfJunctionId[subjectHits(allJunctionToUniqueJunctionOverlap)],names(unlisted_junctions))

  #   intronStartCoordinates = splitAsList(start(uniqueJunctions[uniqueJunctions$mergedHighConfJunctionIdAll_noNA[subjectHits(allJunctionToUniqueJunctionOverlap)]]),names(unlisted_junctions))
  #   intronEndCoordinates = splitAsList(end(uniqueJunctions[uniqueJunctions$mergedHighConfJunctionIdAll_noNA[subjectHits(allJunctionToUniqueJunctionOverlap)]]),names(unlisted_junctions))
  #

  intronStartTMP=start(uniqueJunctions[uniqueJunctions$mergedHighConfJunctionIdAll_noNA[subjectHits(allJunctionToUniqueJunctionOverlap)]])
  intronEndTMP=end(uniqueJunctions[uniqueJunctions$mergedHighConfJunctionIdAll_noNA[subjectHits(allJunctionToUniqueJunctionOverlap)]])
  exon_0size =  which(intronStartTMP[-1]<=intronEndTMP[-length(intronEndTMP)] & names(unlisted_junctions[-1])==names(unlisted_junctions[-length(unlisted_junctions)]))
  if(length(exon_0size)>0) {
    intronStartTMP[-1][exon_0size] <- intronEndTMP[-length(intronEndTMP)][exon_0size]+1
  }
  intronStartCoordinates = splitAsList(intronStartTMP,names(unlisted_junctions))
  intronEndCoordinates = splitAsList(intronEndTMP,names(unlisted_junctions))


  #for comparison, use uncorrected reads
  #intronStartCoordinates = splitAsList(start(uniqueJunctions[subjectHits(allJunctionToUniqueJunctionOverlap)]),names(unlisted_junctions))
  #intronEndCoordinates = splitAsList(end(uniqueJunctions[subjectHits(allJunctionToUniqueJunctionOverlap)]),names(unlisted_junctions))
  #intronEndCoordinates = splitAsList(end(uniqueJunctions)[subjectHits(allJunctionToUniqueJunctionOverlap)],names(unlisted_junctions))


  #annotated strand of junctions for each read based on the infered read strand
  #unlisted_junctions_strand <- as.character(strand(uniqueJunctions)[subjectHits(allJunctionToUniqueJunctionOverlap)])
  ## read strand assignment, is this necessary/helpful? Check again, add evaluation step?
  unlisted_junctions_strand <- uniqueJunctions$strand.mergedHighConfJunction[subjectHits(allJunctionToUniqueJunctionOverlap)]
  unlisted_junctions_strandList = splitAsList(unlisted_junctions_strand,names(unlisted_junctions))
  strandJunctionSum = sum(unlisted_junctions_strandList=='-')-sum(unlisted_junctions_strandList=='+')

  uniqueReadNames <- unique(names(unlisted_junctions))
  readStrand <- rep('*',length(uniqueReadNames))
  names(readStrand) <- uniqueReadNames
  readStrand[names(strandJunctionSum)][strandJunctionSum<0] <- '+'
  readStrand[names(strandJunctionSum)][strandJunctionSum>0] <- '-'
  strand(unlisted_junctions) <- readStrand[names(unlisted_junctions)]

  #
  #
  #   readStrand <- rep('*',length(startRanges))
  #   names(readStrand) <- names(startRanges)
  #   readStrand[names(strandJunctionSum)][strandJunctionSum<0] <- '+'
  #   readStrand[names(strandJunctionSum)][strandJunctionSum>0] <- '-'
  #   strand(unlisted_junctions) <- readStrand[names(unlisted_junctions)]

  # highConfReadSet <- sum(is.na(junctionsByReadListCorrected))>=0  # include all reads for transcript reconstruction, remove unnecessary ones later, or re-assign to other transcripts
  ##### HERE: USE ALL READS, DEAL WITH READ ID, REMOVE START AND END RANGES ####
  # show('number of reads with only high confidence junctions')
  # show(table(highConfReadSet))

  readTable <- tbl_df(data.frame(matrix(ncol=12,nrow=length(uniqueReadNames))))
  colnames(readTable) <- c('chr','start','end','strand', 'blockCount','exonStarts','exonEnds','intronEnds','intronStarts','readClassId', 'readId','confidenceType')

  #uniqueReadNames <- names(which(highConfReadSet))
  readTable[, 'chr']    <-  as.character(unique(seqnames(readGrglist[uniqueReadNames])))  # as.character(unique(seqnames(readGrglist)))
  readTable[, 'start'] <- pmin(min(start(readGrglist[uniqueReadNames])),min(intronStartCoordinates[uniqueReadNames] -2))  # min(start(readGrglist))
  readTable[, 'end']   <- pmax(max(end(readGrglist[uniqueReadNames])),max(intronEndCoordinates[uniqueReadNames]+2))  # max(end(readGrglist))
  readTable[, 'blockCount'] <- elementNROWS(junctionsByReadListCorrected[uniqueReadNames])
  readTable[, 'exonStarts'] <- paste(intronEndCoordinates[uniqueReadNames]+1,collapse=',') #### replace with unstrsplit function (should be faster)
  readTable[, 'exonEnds'] <- paste(intronStartCoordinates[uniqueReadNames]-1,collapse=',')
  readTable[, 'intronEnds'] <- paste(intronEndCoordinates[uniqueReadNames],collapse=',') #### replace with unstrsplit function (should be faster)
  readTable[, 'intronStarts'] <- paste(intronStartCoordinates[uniqueReadNames],collapse=',')
  readTable[, 'strand'] <- readStrand[uniqueReadNames]
  readTable[, 'readId'] <- readNames[as.integer(uniqueReadNames)]
  readTable[, 'confidenceType'] <- 'highConfidenceJunctionReads'
  readTable[, 'readClassId'] <- NA

  #   candReads <- names(which(highConfReadSet))
  #   readTable[, 'chr']    <- as.character(seqnames(startRanges[candReads]) )  # as.character(unique(seqnames(readGrglist)))
  #   readTable[, 'start'] <- pmin(start(startRanges[candReads]),min(intronStartCoordinates[candReads] -2))  # min(start(readGrglist))
  #   readTable[, 'end']   <- pmax(end(endRanges[candReads]),max(intronEndCoordinates[candReads]+2))  # max(end(readGrglist))
  #   readTable[, 'blockCount'] <- elementNROWS(junctionsByReadListCorrected[candReads])
  #   readTable[, 'exonStarts'] <- paste(intronEndCoordinates[candReads]+1,collapse=',') #### replace with unstrsplit function (should be faster)
  #   readTable[, 'exonEnds'] <- paste(intronStartCoordinates[candReads]-1,collapse=',')
  #   readTable[, 'intronEnds'] <- paste(intronEndCoordinates[candReads],collapse=',') #### replace with unstrsplit function (should be faster)
  #   readTable[, 'intronStarts'] <- paste(intronStartCoordinates[candReads],collapse=',')
  #   readTable[, 'strand'] <- readStrand[candReads]
  #   readTable[, 'readId'] <- candReads
  #   readTable[, 'confidenceType'] <- 'highConfidenceJunctionReads'
  #   readTable[, 'readClassId'] <- NA
  #
  #################HERE @################
  #include reads with low confidence junctions, need to be assigned to real transcripts
  readTable[sum(is.na(junctionsByReadListCorrected[uniqueReadNames]))>0, 'confidenceType'] <- 'lowConfidenceJunctionReads'

  ## currently the 90% quantile is used to define the start and end ## this is the slowest part in the function.
  readClassTable <-readTable %>% group_by(chr, strand, intronEnds, intronStarts) %>% mutate(start = quantile(start, 0.1), end = quantile(end, 0.9), readCount=n())%>% distinct(chr, strand, intronEnds, intronStarts, .keep_all=TRUE) %>% ungroup()

  #readClassTable <-readTable %>% group_by(chr, strand, intronEnds, intronStarts) %>% mutate(start = quantile(start, 0.1), end = quantile(end, 0.9), readCount=n())%>% distinct(chr, strand, intronEnds, intronStarts, .keep_all=TRUE) %>% ungroup()

  # required to create exon ranges to build exon By Tx object
  readClassTable <- readClassTable %>%  mutate(exonStarts = paste(as.integer(start), exonStarts,sep=','), exonEnds = paste(exonEnds, as.integer(end),sep=','))


  ### todo here: create txTable after merging all idetnical intron-reads, then reconstruct transcrtipts, and compare with annoated transcripts (using intron ranges only). COmpare this approach to more naive approach, how many distinct RNAs are found, how many reads supprt them etc. It seem like correcting junctions greatly reduves the complexity of transcripts. This appraoch should be very fast in the end
  readClassTable[,'readClassId'] <- paste('rc',1:nrow(readClassTable),sep='.')

  readTable=left_join(select(readTable,readId,chr,strand,intronEnds,intronStarts,confidenceType),select(readClassTable,chr,strand,intronEnds,intronStarts,readClassId)) %>% select(readId,readClassId,confidenceType, strand)

  #remove unecessary information from readTable
  #readTable <- readTable %>% select(readId,readClassId,confidenceType, strand)
  gc()



  exonsByReadClass= with(readClassTable,makeGRangesListFromFeatureFragments(seqnames=chr,
                                                                            fragmentStarts=exonStarts,
                                                                            fragmentEnds=exonEnds,
                                                                            strand=strand))
  names(exonsByReadClass) <- readClassTable$readClassId

  #add exon rank and exon_endRank

  ## combine new transcripts with annotated transcripts based on identical intron pattern
  unlistData = unlist(exonsByReadClass, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(elementNROWS(exonsByReadClass)), names=NULL)

  unlistData$exon_id <- paste0('exId.',1:length(unlistData))
  unlistData$exon_name <- paste0('ex.',1:length(unlistData))
  exon_rank <- sapply(width((partitioning)),seq, from=1)
  exon_rank[which(readClassTable$strand=='-')] <- lapply(exon_rank[which(readClassTable$strand=='-')], rev)  # * assumes positive for exon ranking
  exon_endRank <- lapply(exon_rank, rev)
  unlistData$exon_rank <- unlist(exon_rank)
  unlistData$exon_endRank <- unlist(exon_endRank)

  exonsByReadClass <- relist(unlistData, partitioning)

  readClassTable <- readClassTable %>% select(chr, start, end, strand, readCount, confidenceType, readClassId)

  return(list(exonsByReadClass = exonsByReadClass, readClassTable = readClassTable, readTable = readTable))
}

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
  hitsDFGrouped <- hitsDF %>% group_by(queryHits) %>% mutate(maxStart=max(start), minEnd = min(end)) %>% select(queryHits, chr, maxStart, minEnd, strand) %>% distinct() %>% group_by(chr, maxStart, minEnd, strand) %>% mutate(readClassId = paste0('rc',prefix,'.', group_indices())) %>% ungroup()

  readClassTableUnspliced <- hitsDFGrouped %>%
    select(chr, start=maxStart, end=minEnd, strand, readClassId) %>%
    group_by(readClassId) %>%
    mutate(readCount=n()) %>%
    distinct() %>%
    ungroup() %>%
    mutate(confidenceType=confidenceType) %>%
    select(chr, start, end, strand, readCount, confidenceType, readClassId)

  readTableUnspliced <-  select(hitsDFGrouped, readClassId) %>%
    mutate(confidenceType=confidenceType, strand=as.character(strand(granges[hitsDFGrouped$queryHits])), readId=readNames[as.integer(names(granges[hitsDFGrouped$queryHits]))]) %>%
    select(readId,  readClassId, confidenceType, strand)

  exByReadClassUnspliced <- GRanges(seqnames=readClassTableUnspliced$chr, ranges=IRanges(start=readClassTableUnspliced$start, end=readClassTableUnspliced$end), strand=readClassTableUnspliced$strand)
  exByReadClassUnspliced$exon_id <- paste0('exId',prefix,'.',1:length(exByReadClassUnspliced))
  exByReadClassUnspliced$exon_name <- paste0('ex',prefix,'.', 1:length(exByReadClassUnspliced))
  exByReadClassUnspliced$exon_rank <- 1
  exByReadClassUnspliced$exon_endRank <- 1
  partitioning <- PartitioningByEnd(1:length(exByReadClassUnspliced))
  exByReadClassUnspliced <- relist(exByReadClassUnspliced,partitioning)
  names(exByReadClassUnspliced) <- readClassTableUnspliced$readClassId

  return(list(exonsByReadClass = exByReadClassUnspliced, readClassTable = readClassTableUnspliced, readTable = readTableUnspliced))
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
    hitObject= inner_join(hitObject,hitObject,by=c("subjectHits"="queryHits")) %>% select(queryHits,subjectHits.y) %>% distinct() %>%rename(subjectHits=subjectHits.y)
  }

  geneTxNames <- hitObject %>% group_by(queryHits) %>% mutate(geneId = paste('gene',first(subjectHits),sep='.')) %>% select(queryHits, geneId) %>% distinct()
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
