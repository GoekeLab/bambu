#'@title CREATEJUNCTIONTABLE
#'@description This function creates a table of all junctions from a list grangeslist of junctions, it add strand and splice motif information
#'@param unlisted_junction_granges
#'@param genomeDB
#'@param genomeFA
#'@export
createJunctionTable <- function(unlisted_junction_granges, genomeSequence=NULL) {
  # License note: This function is adopted from the GenomicAlignments package (Author: Hervé Pagès, Valerie Obenchain, Martin Morgan)
  # https://doi.org/doi:10.18129/B9.bioc.GenomicAlignments

  if(is.null(genomeSequence)){
    stop("Reference genome sequence is missing, please provide fasta file or BSgenome name, see available.genomes()")
  }else if(class(genomeSequence) != 'FaFile'){
    if(grepl('.fa',genomeSequence)){
      genomeSequence <- Rsamtools::FaFile(genomeSequence)
    }else {
      if (!suppressWarnings(require(BSgenome, quietly=TRUE)))
        stop("Please install the BSgenome package")

      genomeSequence <- BSgenome::getBSgenome(genomeSequence)
      seqlevelsStyle(genomeSequence) <- seqlevelsStyle(unlisted_junction_granges)[1]
      if(!all(seqlevels(unlisted_junction_granges) %in% seqlevels(genomeSequence))) {
        warning("not all chromosomes present in reference, ranges are dropped")
        unlisted_junction_granges <- keepSeqlevels(unlisted_junction_granges,
                                                   value = seqlevels(unlisted_junction_granges)[seqlevels(unlisted_junction_granges) %in% seqlevels(genomeSequence)],
                                                   pruning.mode = 'coarse')
      }
    }
  }

  ##Todo: don't create junction names, instead work with indices/intergers (names are memory intensive)

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


  junctionSeqStart<-BSgenome::getSeq(genomeSequence,IRanges::shift(flank(uniqueJunctions,width=2),2)) # shift: from IRanges
  junctionSeqEnd<-BSgenome::getSeq(genomeSequence,IRanges::shift(flank(uniqueJunctions,width=2,start=FALSE),-2)) # shift: from IRanges

  junctionMotif <- paste(junctionSeqStart,junctionSeqEnd,sep='-')
  uniqueJunctions_mcols <- cbind(uniqueJunctions_mcols,
                                 DataFrame(spliceMotif = junctionMotif,
                                           spliceStrand = spliceStrand(junctionMotif)))


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
junctionStrandCorrection <- function(uniqueJunctions, unlisted_junction_granges, intronsByTx, stranded, verbose = FALSE) {
  ##note: the strand is not always correctly infered based on motifs, it might introduce systematic errors due to alignment (which is biased towards splice motifs)

  allJunctionToUniqueJunctionOverlap <- findOverlaps(unlisted_junction_granges, uniqueJunctions,type='equal',ignore.strand=T)

  uniqueJunctionsUpdate <- uniqueJunctions # make a copy to revert to if strand correction does not improve results

  annotatedIntronNumber <- evalAnnotationOverlap(uniqueJunctions, intronsByTx, ignore.strand=F)['TRUE']
  if(verbose) {
    show('before strand correction, annotated introns:')
    show(annotatedIntronNumber)
    show(annotatedIntronNumber/length(uniqueJunctions))
  }

  #infer strand for each read based on strand of junctions
  if(stranded==FALSE) {
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

      unstranded_unlisted_junction_granges <- unstrand(unlisted_junction_granges)
      junctionMatchList <- as(findMatches(uniqueJunctions, unstranded_unlisted_junction_granges), "List")
      tmp <- extractList(strand(unlisted_junction_granges), junctionMatchList)
      uniqueJunctionsUpdate$plus_score_inferedByRead <- sum(tmp == "+")
      uniqueJunctionsUpdate$minus_score_inferedByRead <- sum(tmp == "-")
      strandScoreByRead <- uniqueJunctionsUpdate$minus_score_inferedByRead - uniqueJunctionsUpdate$plus_score_inferedByRead
      strand(uniqueJunctionsUpdate[strandScoreByRead< 0 ]) = '+' ## ere we overwrite the information from the motif which increases overlap with known junctions
      strand(uniqueJunctionsUpdate[strandScoreByRead>0 ]) = '-'
      annotatedIntronNumberNew <- evalAnnotationOverlap(uniqueJunctionsUpdate, intronsByTx, ignore.strand=F)['TRUE']
      if(annotatedIntronNumberNew > annotatedIntronNumber & !is.na(annotatedIntronNumber)) # update junctions object if strand prediction improves overlap with annotations
      {
        if(verbose) {
          show('after strand correction, annotated introns:')
          show(annotatedIntronNumberNew)
          show(annotatedIntronNumberNew/length(uniqueJunctionsUpdate))
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
      annotatedIntronNumberNew <- evalAnnotationOverlap(uniqueJunctionsUpdate, intronsByTx, ignore.strand=F)['TRUE']
      if(annotatedIntronNumberNew > annotatedIntronNumber & !is.na(annotatedIntronNumber)) # update junctions object if strand prediction improves overlap with annotations
      {
        if(verbose) {
          show('after strand correction, annotated introns:')
          show(annotatedIntronNumberNew)
          show(annotatedIntronNumberNew/length(uniqueJunctionsUpdate))
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



#'@title EVALANNOTATIONOVERLAP
#'@param intronRanges
#'@param intronsByTx
#'@param ignore.strand
evalAnnotationOverlap <- function(intronRanges, intronsByTx, ignore.strand=FALSE)
{
  return(table(!is.na(GenomicRanges::match(intronRanges, unique(unlist(intronsByTx)),ignore.strand=ignore.strand))))
}



#'@title PREDICTSPLICEJUNCTIONS
#'@description Function to predict splice site as true or false positive based on annotations, requires annotated junctions object, optional list of models learned on the data
#'@param annotatedJunctions
#'@param junctionModel
predictSpliceJunctions <- function(annotatedJunctions, junctionModel=NULL, verbose = FALSE)
{

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
  mySet.training=(annotatedJunctionsStart$annotatedStart.start|annotatedJunctionsStart$annotatedStart)[mySet.all]

  myData=data.frame(annotatedJunctionsStart$startScore/(annotatedJunctionsStart$startScore.start+annotatedJunctionsStart$startScore),annotatedJunctionsStart$startScore,annotatedJunctionsStart$distStart.start,(annotatedJunctionsStart$spliceStrand.start=='+'),annotatedJunctionsStart$spliceStrand.start=='-',(annotatedJunctionsStart$spliceStrand=='+'))[mySet.all,]

  colnames(myData) <- paste('A',1:ncol(myData),sep='.')

  modelmatrix=model.matrix(~A.1+A.2+A.3+A.4+A.5, data=data.frame((myData)))

  if(is.null(junctionModel))
  {


    myResults= fitBinomialModel(labels.train=as.integer(annotatedJunctionsStart$annotatedStart)[mySet.all][mySet.training], data.train=modelmatrix[mySet.training,], data.test=modelmatrix, show.cv=verbose, maxSize.cv=10000)

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


  myData=data.frame(annotatedJunctionsStart$startScore/(annotatedJunctionsStart$startScore.end+annotatedJunctionsStart$startScore),annotatedJunctionsStart$startScore,annotatedJunctionsStart$distStart.end,(annotatedJunctionsStart$spliceStrand.end=='+'),(annotatedJunctionsStart$spliceStrand.end=='-'),(annotatedJunctionsStart$spliceStrand=='+'))[mySet.all,]
  colnames(myData) <- paste('A',1:ncol(myData),sep='.')

  modelmatrix=model.matrix(~A.1+A.2+A.3+A.4+A.5, data=data.frame((myData)))
  if(is.null(junctionModel))
  {


    myResults= fitBinomialModel(labels.train=as.integer(annotatedJunctionsStart$annotatedStart)[mySet.all][mySet.training], data.train=modelmatrix[mySet.training,], data.test=modelmatrix, show.cv=verbose, maxSize.cv=10000)

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


    myResults= fitBinomialModel(labels.train=as.integer(annotatedJunctionsEnd$annotatedEnd)[mySet.all][mySet.training], data.train=modelmatrix[mySet.training,], data.test=modelmatrix, show.cv=verbose, maxSize.cv=10000)
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
  modelmatrix=model.matrix(~A.1+A.2+A.3+A.4+A.5, data=data.frame((myData)))
  if(is.null(junctionModel))
  {


    myResults= fitBinomialModel(labels.train=as.integer(annotatedJunctionsEnd$annotatedEnd)[mySet.all][mySet.training], data.train=modelmatrix[mySet.training,], data.test=modelmatrix, show.cv=verbose, maxSize.cv=10000)
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



#'@title FITBINOMIALMODEL
#'@param labels.train
#'@param data.train
#'@param data.test
#'@param show.cv
#'@param maxSize.cv
fitBinomialModel <- function(labels.train, data.train, data.test, show.cv=TRUE, maxSize.cv=10000, ...)
{
  if(show.cv)
  {
    mySample=sample(1:length(labels.train),min(floor(length(labels.train)/2),maxSize.cv))
    data.train.cv=data.train[mySample,]
    labels.train.cv=labels.train[mySample]
    data.train.cv.test=data.train[-mySample,]
    labels.train.cv.test=labels.train[-mySample]

    cv.fit=glmnet::cv.glmnet(x=data.train.cv,y=labels.train.cv,family='binomial', ...)
    predictions=predict(cv.fit,newx=data.train.cv.test,s='lambda.min')
    show('prediction accuracy (CV) (higher for splice donor than splice acceptor)')

    show( fisher.test(table(predictions>0,labels.train.cv.test)))
    show(myPerformance(labels.train.cv.test==1,predictions)$AUC	)
  }

  cv.fit=glmnet::cv.glmnet(x=data.train,y=labels.train,family='binomial', ...)
  predictions= predict(cv.fit,newx=data.test,s='lambda.min')
  return(list(predictions,cv.fit))
}





#'@title FINDHIGHCONFIDENTJUNCTIONS
#'@description this function adds "mergedHighConfJunctionId" to the junciton list which contains the ID of the most likely high confident junction that each junction originates from
#'@param junctions
#'@param junctionModel
findHighConfidenceJunctions <- function(junctions, junctionModel, verbose = FALSE) {
  if(verbose) {
    show('reads count for all annotated junctions')
    show(sum(junctions$score[junctions$annotatedJunction]))
    show(sum(junctions$score[junctions$annotatedJunction])/ sum(junctions$score))
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
      if(verbose) show(maxDist)
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
    show('reads count for all annotated junctions after correction to reference junction')
    sumByJuncId <- tapply(junctions$score, junctions$mergedHighConfJunctionId, sum)
    show(sum(sumByJuncId[junctions[names(sumByJuncId)]$annotatedJunction]))
    show(sum(sumByJuncId[junctions[names(sumByJuncId)]$annotatedJunction])/sum(junctions$score))
  }
  return(junctions)
}
