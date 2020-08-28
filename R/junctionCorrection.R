#' Create Junction tables from unlisted junction granges
#' @importFrom BiocParallel bppram bpvec
#' @noRd
createJunctionTable <- function(unlisted_junction_granges, genomeSequence=NULL, ncore=1) {
  # License note: This function is adopted from the GenomicAlignments package (Author: Hervé Pagès, Valerie Obenchain, Martin Morgan)
  # https://doi.org/doi:10.18129/B9.bioc.GenomicAlignments

  if(is.null(genomeSequence)){
    stop("Reference genome sequence is missing, please provide fasta file or BSgenome name, see available.genomes()")
  }
  
  original_seqlevelstyle <- seqlevelsStyle(unlisted_junction_granges)[1]
  
  if(is(genomeSequence,'character')){
    if(grepl('.fa',genomeSequence)){
      
      if(.Platform$OS.type == "windows"){
        genomeSequence <- Biostrings::readDNAStringSet(genomeSequence)
        newlevels <- unlist(lapply(strsplit(names(genomeSequence)," "),"[[",1))  
        names(genomeSequence) <- newlevels
      }else{
        genomeSequence <- Rsamtools::FaFile(genomeSequence)
      }
      if(seqlevelsStyle(genomeSequence)[1]  != seqlevelsStyle(unlisted_junction_granges)[1]){
        seqlevelsStyle(unlisted_junction_granges) <- seqlevelsStyle(genomeSequence)[1] 
      }
    } else {
      genomeSequence <- BSgenome::getBSgenome(genomeSequence)
      seqlevelsStyle(genomeSequence) <- seqlevelsStyle(unlisted_junction_granges)[1]
    }
  } else if(is(genomeSequence,'BSgenome')){
    seqlevelsStyle(genomeSequence) <- seqlevelsStyle(unlisted_junction_granges)[1]
  } 
  
  if(is(genomeSequence,'FaFile')){
    if(seqlevelsStyle(genomeSequence)[1]  != seqlevelsStyle(unlisted_junction_granges)[1]){
      seqlevelsStyle(unlisted_junction_granges) <- seqlevelsStyle(genomeSequence)[1] 
    }
  }
  if(!all(seqlevels(unlisted_junction_granges) %in% seqlevels(genomeSequence))) {
    message("not all chromosomes present in reference genome sequence, ranges are dropped")
    unlisted_junction_granges <- keepSeqlevels(unlisted_junction_granges,
                                               value = seqlevels(unlisted_junction_granges)[seqlevels(unlisted_junction_granges) %in% seqlevels(genomeSequence)],
                                               pruning.mode = 'coarse')
  }
  
  ##Todo: don't create junction names, instead work with indices/intergers (names are memory intensive)

  unstranded_unlisted_junctions <- unstrand(unlisted_junction_granges)
  uniqueJunctions <- sort(unique(unstranded_unlisted_junctions))
  names(uniqueJunctions) <- paste('junc',seq_along(uniqueJunctions),sep='.') ##replace with integer later (or only use indices if possible)

  #calculate stranded read counts
  junctionMatchList <- as(findMatches(uniqueJunctions, unstranded_unlisted_junctions), "List")
  uniqueJunctions_score <- elementNROWS(junctionMatchList)
  junctionStrandList <- extractList(strand(unlisted_junction_granges), junctionMatchList)


  # the code now includes a parallel implementation, which is only helpful when the BSgenome package is used

  junctionSeqStart<-BSgenome::getSeq(genomeSequence,IRanges::shift(flank(uniqueJunctions,width=2),2)) # shift: from IRanges
  junctionSeqEnd<-BSgenome::getSeq(genomeSequence,IRanges::shift(flank(uniqueJunctions,width=2,start=FALSE),-2)) # shift: from IRanges

  # bpParameters <- BiocParallel::bpparam()
  # bpParameters$workers <- ncore
  # junctionSeqStart <- BiocParallel::bpvec(IRanges::shift(flank(uniqueJunctions,width=2),2),
  #                           BSgenome::getSeq,
  #                           x = genomeSequence,
  #                           BPPARAM=bpParameters)
  # 
  # junctionSeqEnd <- BiocParallel::bpvec(IRanges::shift(flank(uniqueJunctions,width=2,start=FALSE),-2),
  #                                       BSgenome::getSeq,
  #                         x = genomeSequence,
  #                         BPPARAM=bpParameters)


  junctionMotif <- paste(junctionSeqStart,junctionSeqEnd,sep='-')

  junctionStartName <- paste(seqnames(uniqueJunctions),start(uniqueJunctions),sep=':')
  junctionEndName <- paste(seqnames(uniqueJunctions),end(uniqueJunctions),sep=':')

  startScore <- as.integer(tapply(uniqueJunctions_score,  junctionStartName ,sum)[junctionStartName])
  endScore <- as.integer(tapply(uniqueJunctions_score,  junctionEndName ,sum)[junctionEndName])

  mcols(uniqueJunctions) <- DataFrame(score=uniqueJunctions_score,
                                      plus_score=sum(junctionStrandList == "+"),
                                      minus_score=sum(junctionStrandList == "-"),
                                      spliceMotif = junctionMotif,
                                      spliceStrand = spliceStrand(junctionMotif),
                                      junctionStartName = junctionStartName,
                                      junctionEndName = junctionEndName,
                                      startScore = startScore,
                                      endScore=endScore,
                                      id=seq_along(uniqueJunctions))

  strand(uniqueJunctions)<-uniqueJunctions$spliceStrand
  if(original_seqlevelstyle != seqlevelsStyle(unlisted_junction_granges)[1]){
    seqlevelsStyle(uniqueJunctions) <- original_seqlevelstyle
  }
  return(uniqueJunctions)
}



#' JUNCTIONSTRANDCORRECTION
#' @noRd
junctionStrandCorrection <- function(uniqueJunctions, unlisted_junction_granges, uniqueAnnotatedIntrons, stranded, verbose = FALSE) {
  ##note: the strand is not always correctly infered based on motifs, it might introduce systematic errors due to alignment (which is biased towards splice motifs)

  allJunctionToUniqueJunctionOverlap <- findOverlaps(unlisted_junction_granges, uniqueJunctions,type='equal',ignore.strand=TRUE)

  uniqueJunctionsUpdate <- uniqueJunctions # make a copy to revert to if strand correction does not improve results

  annotatedIntronNumber <- evalAnnotationOverlap(uniqueJunctions, uniqueAnnotatedIntrons, ignore.strand=FALSE)['TRUE']
  if(verbose) {
    message('before strand correction, annotated introns:')
    message(annotatedIntronNumber)
    message(annotatedIntronNumber/length(uniqueJunctions))
  }

  #infer strand for each read based on strand of junctions
  if(stranded==FALSE) {
    strandStep <- TRUE
    while(strandStep) { # iterate this 2 times to predict strand more accurately using the mean junction counts

      #annotated strand of jucntions for each read based on the infered read strand
      unlisted_junction_granges_strand <- as.character(strand(uniqueJunctionsUpdate)[subjectHits(allJunctionToUniqueJunctionOverlap)])
      unlisted_junction_granges_strandList = splitAsList(unlisted_junction_granges_strand,mcols(unlisted_junction_granges)$id)
      strandJunctionSum = sum(unlisted_junction_granges_strandList=='-')-sum(unlisted_junction_granges_strandList=='+')

      readStrand <- rep('*',length(unlisted_junction_granges_strandList))
      readStrand[strandJunctionSum<0] <- '+'
      readStrand[strandJunctionSum>0] <- '-'
      strand(unlisted_junction_granges) <- readStrand[match(mcols(unlisted_junction_granges)$id,
                                                            as.integer(names(unlisted_junction_granges_strandList)))]


      unstranded_unlisted_junction_granges <- unstrand(unlisted_junction_granges)
      junctionMatchList <- as(findMatches(uniqueJunctions, unstranded_unlisted_junction_granges), "List")
      tmp <- extractList(strand(unlisted_junction_granges), junctionMatchList)
      uniqueJunctionsUpdate$plus_score_inferedByRead <- sum(tmp == "+")
      uniqueJunctionsUpdate$minus_score_inferedByRead <- sum(tmp == "-")
      strandScoreByRead <- uniqueJunctionsUpdate$minus_score_inferedByRead - uniqueJunctionsUpdate$plus_score_inferedByRead
      strand(uniqueJunctionsUpdate[strandScoreByRead< 0 ]) = '+' ## here we overwrite the information from the motif which increases overlap with known junctions
      strand(uniqueJunctionsUpdate[strandScoreByRead>0 ]) = '-'
      annotatedIntronNumberNew <- evalAnnotationOverlap(uniqueJunctionsUpdate, uniqueAnnotatedIntrons, ignore.strand=FALSE)['TRUE']
      if(annotatedIntronNumberNew > annotatedIntronNumber & !is.na(annotatedIntronNumber)) # update junctions object if strand prediction improves overlap with annotations
      {
        if(verbose) {
          message('after strand correction, annotated introns:')
          message(annotatedIntronNumberNew)
          message(annotatedIntronNumberNew/length(uniqueJunctionsUpdate))
        }
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
      strand(uniqueJunctionsUpdate[strandScoreByRead>0 ]) = '-' ## here we verwrite the information from the motif which increases overlap with known junctions
      strand(uniqueJunctionsUpdate[strandScoreByRead<(0) ]) = '+'
      annotatedIntronNumberNew <- evalAnnotationOverlap(uniqueJunctionsUpdate, uniqueAnnotatedIntrons, ignore.strand=FALSE)['TRUE']
      if(annotatedIntronNumberNew > annotatedIntronNumber & !is.na(annotatedIntronNumber)) # update junctions object if strand prediction improves overlap with annotations
      {
        if(verbose) {
          message('after strand correction, annotated introns:')
          message(annotatedIntronNumberNew)
          message(annotatedIntronNumberNew/length(uniqueJunctionsUpdate))
        }
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


#' Evaluate annoation overlap
#' @noRd
evalAnnotationOverlap <- function(intronRanges, uniqueAnnotatedIntrons, ignore.strand=FALSE)
{
  return(table(!is.na(GenomicRanges::match(intronRanges, uniqueAnnotatedIntrons,ignore.strand=ignore.strand))))
}


#' Predict splicing junctions
#' @noRd
predictSpliceJunctions <- function(annotatedJunctions, junctionModel=NULL, verbose = FALSE)
{
  ##note: readibility can be improved using dplyr
  annotatedJunctionsStart <- unique(GRanges(seqnames=seqnames(annotatedJunctions),
                                            ranges=IRanges(start=start(annotatedJunctions),
                                                           end=start(annotatedJunctions)),
                                            strand='*',
                                            mcols(annotatedJunctions)[,c('startScore','junctionStartName','annotatedStart','spliceStrand','spliceMotif')]))

  annotatedJunctionsStart$distStart.start <- c(0,(start(annotatedJunctionsStart[-1])-start(annotatedJunctionsStart[-length(annotatedJunctionsStart)]))*as.integer((seqnames(annotatedJunctionsStart[-1])==seqnames(annotatedJunctionsStart[-length(annotatedJunctionsStart)]))))
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
  annotatedJunctions$spliceSitePredictionStart.start <- rep(NA, length(annotatedJunctions))
  if(any(mySet.all)){
    mySet.training=(annotatedJunctionsStart$annotatedStart.start|annotatedJunctionsStart$annotatedStart)[mySet.all]
    
    myData=data.frame(annotatedJunctionsStart$startScore/(annotatedJunctionsStart$startScore.start+annotatedJunctionsStart$startScore),annotatedJunctionsStart$startScore,annotatedJunctionsStart$distStart.start,(annotatedJunctionsStart$spliceStrand.start=='+'),annotatedJunctionsStart$spliceStrand.start=='-',(annotatedJunctionsStart$spliceStrand=='+'))[mySet.all,]
    
    colnames(myData) <- paste('A',seq_len(ncol(myData)),sep='.')
    
    modelmatrix=model.matrix(~A.1+A.2+A.3+A.4+A.5, data=data.frame((myData)))
    
    if(is.null(junctionModel))
    {
      
      
      myResults= fitBinomialModel(labels.train=as.integer(annotatedJunctionsStart$annotatedStart)[mySet.all][mySet.training], data.train=modelmatrix[mySet.training,], data.test=modelmatrix, show.cv=verbose, maxSize.cv=10000)
      
      junctionModelList[['spliceSitePredictionStart.start']] <- myResults[[2]]
      predictions = myResults[[1]]
    }
    else{
      predictions= glmnet:::predict.cv.glmnet(junctionModel[['spliceSitePredictionStart.start']],newx=modelmatrix,s='lambda.min')
    }
    spliceSitePredictionStart.start <- rep(NA, length(annotatedJunctionsStart))
    names(spliceSitePredictionStart.start) <- annotatedJunctionsStart$junctionStartName
    spliceSitePredictionStart.start[mySet.all]=predictions
    
    annotatedJunctions$spliceSitePredictionStart.start <- spliceSitePredictionStart.start[annotatedJunctions$junctionStartName]
  }
  
  ## test start splice site given close by splice site (right/3')
  
  mySet.all=((annotatedJunctionsStart$distStart.end!=0)&annotatedJunctionsStart$spliceStrand!='*'&annotatedJunctionsStart$startScore>0&abs(annotatedJunctionsStart$distStart.end)<15)
  annotatedJunctions$spliceSitePredictionStart.end <- rep(NA, length(annotatedJunctions))
  if(any(mySet.all)){
    mySet.training=(annotatedJunctionsStart$annotatedStart.end|annotatedJunctionsStart$annotatedStart)[mySet.all]
    
    
    myData=data.frame(annotatedJunctionsStart$startScore/(annotatedJunctionsStart$startScore.end+annotatedJunctionsStart$startScore),annotatedJunctionsStart$startScore,annotatedJunctionsStart$distStart.end,(annotatedJunctionsStart$spliceStrand.end=='+'),(annotatedJunctionsStart$spliceStrand.end=='-'),(annotatedJunctionsStart$spliceStrand=='+'))[mySet.all,]
    colnames(myData) <- paste('A',seq_len(ncol(myData)),sep='.')
    
    modelmatrix=model.matrix(~A.1+A.2+A.3+A.4+A.5, data=data.frame((myData)))
    if(is.null(junctionModel))
    {
      
      
      myResults= fitBinomialModel(labels.train=as.integer(annotatedJunctionsStart$annotatedStart)[mySet.all][mySet.training], data.train=modelmatrix[mySet.training,], data.test=modelmatrix, show.cv=verbose, maxSize.cv=10000)
      
      junctionModelList[['spliceSitePredictionStart.end']] <- myResults[[2]]
      predictions =myResults[[1]]
    }
    else{
      
      predictions= glmnet:::predict.cv.glmnet(junctionModel[['spliceSitePredictionStart.end']],newx=modelmatrix,s='lambda.min')
    }
    spliceSitePredictionStart.end <- rep(NA, length(annotatedJunctionsStart))
    names(spliceSitePredictionStart.end) <- annotatedJunctionsStart$junctionStartName
    spliceSitePredictionStart.end[mySet.all]=predictions
    
    annotatedJunctions$spliceSitePredictionStart.end <- spliceSitePredictionStart.end[annotatedJunctions$junctionStartName]
  }
  
  ## test end splice site given close by splice site (start/5')
  
  mySet.all=(annotatedJunctionsEnd$distEnd.start!=0&annotatedJunctionsEnd$spliceStrand!='*'&annotatedJunctionsEnd$endScore>0&(annotatedJunctionsEnd$distEnd.start<15))
  annotatedJunctions$spliceSitePredictionEnd.start <- rep(NA, length(annotatedJunctions))
  
  if(any(mySet.all)){
    mySet.training=(annotatedJunctionsEnd$annotatedEnd.start|annotatedJunctionsEnd$annotatedEnd)[mySet.all]
    
    
    myData=data.frame(annotatedJunctionsEnd$endScore/(annotatedJunctionsEnd$endScore.start+annotatedJunctionsEnd$endScore),annotatedJunctionsEnd$endScore,annotatedJunctionsEnd$distEnd.start,(annotatedJunctionsEnd$spliceStrand.start=='+'),annotatedJunctionsEnd$spliceStrand.start=='-',(annotatedJunctionsEnd$spliceStrand=='+'))[mySet.all,]
    
    colnames(myData) <- paste('A',seq_len(ncol(myData)),sep='.')
    modelmatrix=model.matrix(~A.1+A.2+A.3+A.4+A.5, data=data.frame((myData))) #+A.6+A.7+A.8+A.9+A.10
    if(is.null(junctionModel))
    {
      
      
      myResults= fitBinomialModel(labels.train=as.integer(annotatedJunctionsEnd$annotatedEnd)[mySet.all][mySet.training], data.train=modelmatrix[mySet.training,], data.test=modelmatrix, show.cv=verbose, maxSize.cv=10000)
      junctionModelList[['spliceSitePredictionEnd.start']] <- myResults[[2]]
      predictions =myResults[[1]]
    }
    else{
      predictions= glmnet:::predict.cv.glmnet(junctionModel[['spliceSitePredictionEnd.start']],newx=modelmatrix,s='lambda.min')
    }
    spliceSitePredictionEnd.start <- rep(NA, length(annotatedJunctionsEnd))
    names(spliceSitePredictionEnd.start) <- annotatedJunctionsEnd$junctionEndName
    spliceSitePredictionEnd.start[mySet.all]=predictions
    annotatedJunctions$spliceSitePredictionEnd.start <- spliceSitePredictionEnd.start[annotatedJunctions$junctionEndName]
  }
  ## test end splice site given close by splice site (right/3')
  
  mySet.all=(annotatedJunctionsEnd$distEnd.end!=0&annotatedJunctionsEnd$spliceStrand!='*'&annotatedJunctionsEnd$endScore>0&abs(annotatedJunctionsEnd$distEnd.end)<15)
  annotatedJunctions$spliceSitePredictionEnd.end <- rep(NA, length(annotatedJunctions))
  
  if(any(mySet.all)){
    mySet.training=(annotatedJunctionsEnd$annotatedEnd.end|annotatedJunctionsEnd$annotatedEnd)[mySet.all]
    
    
    myData=data.frame(annotatedJunctionsEnd$endScore/(annotatedJunctionsEnd$endScore.end+annotatedJunctionsEnd$endScore),annotatedJunctionsEnd$endScore,annotatedJunctionsEnd$distEnd.end ,annotatedJunctionsEnd$distEnd.end,(annotatedJunctionsEnd$spliceStrand.end=='+'),(annotatedJunctionsEnd$spliceStrand.end=='-'),(annotatedJunctionsEnd$spliceStrand=='+'))[mySet.all,]
    
    colnames(myData) <- paste('A',seq_len(ncol(myData)),sep='.')
    modelmatrix=model.matrix(~A.1+A.2+A.3+A.4+A.5, data=data.frame((myData)))
    if(is.null(junctionModel))
    {
      
      
      myResults= fitBinomialModel(labels.train=as.integer(annotatedJunctionsEnd$annotatedEnd)[mySet.all][mySet.training], data.train=modelmatrix[mySet.training,], data.test=modelmatrix, show.cv=verbose, maxSize.cv=10000)
      junctionModelList[['spliceSitePredictionEnd.end']] <- myResults[[2]]
      predictions =myResults[[1]]
    }
    else{
      predictions= glmnet:::predict.cv.glmnet(junctionModel[['spliceSitePredictionEnd.end']],newx=modelmatrix,s='lambda.min')
    }
    spliceSitePredictionEnd.end <- rep(NA, length(annotatedJunctionsEnd))
    names(spliceSitePredictionEnd.end) <- annotatedJunctionsEnd$junctionEndName
    spliceSitePredictionEnd.end[mySet.all]=predictions
    
    annotatedJunctions$spliceSitePredictionEnd.end <- spliceSitePredictionEnd.end[annotatedJunctions$junctionEndName]
  }
  
  if(is.null(junctionModel))
  { junctionModel=junctionModelList}
  return(list(annotatedJunctions, junctionModel))
}


#' Fit binomial model
#' @noRd
fitBinomialModel <- function(labels.train, data.train, data.test, show.cv=TRUE, maxSize.cv=10000)
{
  if(show.cv)
  {
    mySample=sample(seq_along(labels.train),min(floor(length(labels.train)/2),maxSize.cv))
    data.train.cv=data.train[mySample,]
    labels.train.cv=labels.train[mySample]
    data.train.cv.test=data.train[-mySample,]
    labels.train.cv.test=labels.train[-mySample]

    cv.fit=glmnet::cv.glmnet(x=data.train.cv,y=labels.train.cv,family='binomial')
    predictions=glmnet:::predict.cv.glmnet(cv.fit,newx=data.train.cv.test,s='lambda.min')
    message('prediction accuracy (CV) (higher for splice donor than splice acceptor)')

    testResults <- fisher.test(table(predictions>0,labels.train.cv.test))
    show(testResults$estimate)
    show(testResults$p.value)
    show(evalutePerformance(labels.train.cv.test==1,predictions)$AUC)
  }

  cv.fit=glmnet::cv.glmnet(x=data.train,y=labels.train,family='binomial')
  predictions= glmnet:::predict.cv.glmnet(cv.fit,newx=data.test,s='lambda.min')
  return(list(predictions,cv.fit))
}




#'  this function adds "mergedHighConfJunctionId" to the junciton list which contains the ID of the most likely high confident junction that each junction originates from
#' @noRd
findHighConfidenceJunctions <- function(junctions, junctionModel, verbose = FALSE) {
  if(verbose) {
    message('reads count for all annotated junctions')
    message(sum(junctions$score[junctions$annotatedJunction]))
    message(sum(junctions$score[junctions$annotatedJunction])/ sum(junctions$score))
  }
  ##note: the output can be visualised (bed/bigbed track)
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

    ## note: for evaluation/comparison, annotations should not be used here, they only add little extra information
    setReferenceJunctionsPlus=    ((highConfJunctionsPlus$spliceSitePredictionStart.start>0|is.na(highConfJunctionsPlus$spliceSitePredictionStart.start)) & (highConfJunctionsPlus$spliceSitePredictionStart.end>0 | is.na(highConfJunctionsPlus$spliceSitePredictionStart.end)) & (highConfJunctionsPlus$spliceSitePredictionEnd.start>0 | is.na(highConfJunctionsPlus$spliceSitePredictionEnd.start)) & (highConfJunctionsPlus$spliceSitePredictionEnd.end>0 | is.na(highConfJunctionsPlus$spliceSitePredictionEnd.end))) | highConfJunctionsPlus$annotatedJunction

    setReferenceJunctionsMinus=   ((highConfJunctionsMinus$spliceSitePredictionStart.start>0|is.na(highConfJunctionsMinus$spliceSitePredictionStart.start)) & (highConfJunctionsMinus$spliceSitePredictionStart.end>0 | is.na(highConfJunctionsMinus$spliceSitePredictionStart.end)) & (highConfJunctionsMinus$spliceSitePredictionEnd.start>0 | is.na(highConfJunctionsMinus$spliceSitePredictionEnd.start)) & (highConfJunctionsMinus$spliceSitePredictionEnd.end>0 | is.na(highConfJunctionsMinus$spliceSitePredictionEnd.end))) |  highConfJunctionsMinus$annotatedJunction


    candidateJunctionsPlus$highConfJunctionPrediction[highConfidentJunctionSetPlus] <- setReferenceJunctionsPlus
    candidateJunctionsMinus$highConfJunctionPrediction[highConfidentJunctionSetMinus] <- setReferenceJunctionsMinus

    mergedHighConfJunctionIdPlus <- rep(NA,length(candidateJunctionsPlus))
    mergedHighConfJunctionIdMinus <- rep(NA,length(candidateJunctionsMinus))



    ##max Distance can be a parameter that can be set by users
    #here: assign reference junction to all junctions based on prediciton score
    for(maxDist in 0:10)
    {
      if(verbose) message(maxDist)
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

    mergedHighConfJunctionId <- rep(NA, length(junctions))
    names(mergedHighConfJunctionId) <- names(junctions)
    mergedHighConfJunctionId[names(candidateJunctionsPlus)] <- candidateJunctionsPlus$mergedHighConfJunctionId
    mergedHighConfJunctionId[names(candidateJunctionsMinus)] <- candidateJunctionsMinus$mergedHighConfJunctionId

    junctions$mergedHighConfJunctionId <-  as.character(mergedHighConfJunctionId)
  } else {
    warning('no junction correction as no high confidence reference junctions found')
    junctions$mergedHighConfJunctionId <-  names(junctions)
  }
  rm(candidateJunctionsPlus, candidateJunctionsMinus, highConfJunctionsPlus, highConfJunctionsMinus)

  if(verbose) {
    message('reads count for all annotated junctions after correction to reference junction')
    sumByJuncId <- tapply(junctions$score, junctions$mergedHighConfJunctionId, sum)
    message(sum(sumByJuncId[junctions[names(sumByJuncId)]$annotatedJunction]))
    message(sum(sumByJuncId[junctions[names(sumByJuncId)]$annotatedJunction])/sum(junctions$score))
  }
  return(junctions[,'mergedHighConfJunctionId'])
}


#' Evaluate performance
#' @noRd
evalutePerformance <- function(labels, scores, decreasing = TRUE){
  labels <- labels[order(scores, decreasing = decreasing)]
  results <- list()
  results[['TPR']] <- cumsum(labels)/sum(labels)  # TP/(TP+FP); True Positive Rate;Sensitivity; recall
  results[['FPR']] <- cumsum(!labels)/sum(!labels)  # FP/(FP+TN); False Positive Rate;1-Specificity
  results[['precision']] <- cumsum(labels)/(seq_along(labels))  # TP/(TP+FP); positive predictive value;precision
  results[['AUC']] <- sum(results[['TPR']][!duplicated( results[['FPR']], fromLast=TRUE)]/sum(!duplicated( results[['FPR']], fromLast=TRUE)))
  return(results)
}
