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

