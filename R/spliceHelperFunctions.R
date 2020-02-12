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
  equal <- (!is.na(S4Vectors::match(olap ,olapEqual)))
  unique <- myOneMatch(compatible, queryHits(olap))
  strandSpecific <- all(strand(query) != "*")
  mcols(olap) <- DataFrame(compatible, equal, unique, strandSpecific)
  olap
}

#'@title MYCOMPATIBLETRANSCRIPTION
#'@param splice
#'@param query
#'@param subject
myCompatibleTranscription <- function(query, subject, splice)
{
  qrng <- ranges(query)
  srng <- ranges(subject)
  sprng <- ranges(splice)

  bnds <- elementNROWS(GenomicRanges::setdiff(qrng, srng)) == 0L
  splc <- elementNROWS(GenomicRanges::intersect(srng, sprng)) == 0L
  return(bnds & splc)
}

#'@title MYONEMATCH
#'@param x
#'@param idx
myOneMatch <- function(x, idx)
{
  xcnt <- rowsum(as.integer(x), idx)[,1]
  oneMatch <- rep((xcnt == 1L), table(idx))
  unname(x & oneMatch)
}
.isNumericOrNAs <- S4Vectors:::isNumericOrNAs

#'@title MYGAPS
#'@param x
#'@param start
#'@param end
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


#'@title .EXTRACT_UNORIENTED_INTRON_MOTIF
#'@param genome
#'@param junctions
.extract_unoriented_intron_motif <- function(genome, junctions)
{
  mcols(junctions) <- NULL
  junctions_len <- length(junctions)
  Ldinucl_gr <- Rdinucl_gr <- junctions
  end(Ldinucl_gr) <- start(Ldinucl_gr) + 1L
  start(Rdinucl_gr) <- end(Rdinucl_gr) - 1L
  all_dinucl <- getSeq(genome, c(Ldinucl_gr, Rdinucl_gr))
  Rdinucl <- tail(all_dinucl, n=junctions_len)
  xscat(Ldinucl, "-", Rdinucl)
}

#'@title SPLICESTRAND
#'@param motif
spliceStrand <- function(motif){
  NATURAL_INTRON_MOTIFS_RC <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(GenomicAlignments::NATURAL_INTRON_MOTIFS)))

  motifStrand <- ifelse(motif %in% GenomicAlignments::NATURAL_INTRON_MOTIFS,'+','*')
  motifStrand[motif %in% NATURAL_INTRON_MOTIFS_RC] <- '-'
  return(motifStrand)
}


#'@title EVALANNOTATIONOVERLAP
#'@param intronRanges
#'@param intronsByTx
#'@param ignore.strand
evalAnnotationOverlap <- function(intronRanges, intronsByTx, ignore.strand=FALSE)
{
  return(table(!is.na(GenomicRanges::match(intronRanges, unique(unlist(intronsByTx)),ignore.strand=ignore.strand))))
}
#'@title TRIMFIRSTLASTEXONS
#'@description function returns exon granges list where the first and last exons are trimmed to width 2, all introns are preserved, requires exon_rank and exon_endRank
#'@param grangesListWithExonRanks
trimFirstLastExons <- function(grangesListWithExonRanks) {
  unlistedGrangesList <- unlist(grangesListWithExonRanks, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(elementNROWS(grangesListWithExonRanks)), names=NULL)
  startExonsSet <- (which((unlistedGrangesList$exon_rank==1 & as.character(strand(unlistedGrangesList))!='-')|(unlistedGrangesList$exon_endRank==1 & as.character(strand(unlistedGrangesList))=='-')))
  endExonsSet <- (which((unlistedGrangesList$exon_rank==1 & as.character(strand(unlistedGrangesList))=='-')|(unlistedGrangesList$exon_endRank==1 & as.character(strand(unlistedGrangesList))!='-')))
  start(unlistedGrangesList[startExonsSet]) <-  end(unlistedGrangesList[startExonsSet])-1
  end(unlistedGrangesList[endExonsSet]) <-  start(unlistedGrangesList[endExonsSet])+1
  return(relist(unlistedGrangesList, partitioning))
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

#'@title PREDICTSPLICEJUNCTIONS
#'@description Function to predict splice site as true or false positive based on annotations, requires annotated junctions object, optional list of models learned on the data
#'@param annotatedJunctions
#'@param junctionModel
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

#'@title MYPERFORMANCE
#'@param labels
#'@param scores
#'@param descreasing
myPerformance <- function(labels, scores, decreasing = TRUE){
  labels <- labels[order(scores, decreasing = decreasing)]
  results <- list()
  results[['TPR']] <- cumsum(labels)/sum(labels)  # TP/(TP+FP); True Positive Rate;Sensitivity; recall
  results[['FPR']] <- cumsum(!labels)/sum(!labels)  # FP/(FP+TN); False Positive Rate;1-Specificity
  results[['precision']] <- cumsum(labels)/(1:length(labels))  # TP/(TP+FP); positive predictive value;precision
  results[['AUC']] <- sum(results[['TPR']][!duplicated( results[['FPR']],fromLast=TRUE)]/sum(!duplicated( results[['FPR']],fromLast=TRUE)))
  return(results)
}


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
grangesListToBed<-function(x) # note: 0 based coordinates
{
  # useful: generate UCSC custom track with the following steps


  xUnlist <- unlist(range(x))

  bedData = data.frame(as.character(seqnames(xUnlist)),start(xUnlist)-1,end(xUnlist),names(xUnlist),rep(1000,length(xUnlist)),as.character(strand(xUnlist)),start(xUnlist)-1,end(xUnlist),rep(0,length(xUnlist)),elementNROWS(x),unlist(lapply(width(x),paste,collapse=',')),unlist(lapply((start(x)-min(start(x))),paste,collapse=','))) #, stringsAsFactors=FALSE
  colnames(bedData) <-  c('chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts')
  return(bedData)
}
