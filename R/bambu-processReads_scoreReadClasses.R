#' Assigns each read class geneScore and txScore
#' @param se summerized experiment object with read classes/ranges
#' @param genomeSequence genomeSequence
#' @param annotations GRangesList of annotations
scoreReadClasses = function(se, genomeSequence, annotations, 
                                    min.readCount = 2){
    options(scipen = 999) #maintain numeric basepair locations not sci.notfi.
    rowData(se)$GENEID = assignGeneIds(rowRanges(se), annotations)
    countsTBL = calculateGeneProportion(counts=mcols(se)$readCount,
                                        geneIds=mcols(se)$GENEID)
    rowData(se)$geneReadProp = countsTBL$geneReadProp
    rowData(se)$geneReadCount = countsTBL$geneReadCount
    thresholdIndex = which(rowData(se)$readCount
                         >=min.readCount)
    newRowData = addRowData(se[thresholdIndex,] , genomeSequence, annotations)
    rowData(se)[names(newRowData)] = NA
    rowData(se)[thresholdIndex,names(newRowData)] = newRowData
    geneScore = getGeneScore(se[thresholdIndex,])
    rowData(se)$geneScore = rep(0,nrow(se))
    rowData(se)$geneFDR = rep(0,nrow(se))
    rowData(se)$geneScore[thresholdIndex] = geneScore$geneScore
    rowData(se)$geneFDR[thresholdIndex] = geneScore$geneFDR
    txIndex = which(rowData(se)$readCount
                    >=min.readCount & !rowData(se)$novel)
    txScore = getTranscriptScore(se[txIndex,])
    rowData(se)$txScore = rep(0,nrow(se))
    rowData(se)$txFDR = rep(0,nrow(se))
    rowData(se)$txScore[thresholdIndex] = txScore$txScore
    rowData(se)$txFDR[thresholdIndex] = txScore$txFDR
    
    return(se)
}

#' calculates labels and features used in model generation
addRowData = function(se, genomeSequence, annotations){
  compTable <- isReadClassCompatible(rowRanges(se), annotations)
  polyATerminals = countPolyATerminals(rowRanges(se), genomeSequence)

  rowData = data.frame(numExons = elementNROWS(rowRanges(se)),
                          equal = compTable$equal,
                          compatible = compTable$compatible,
                          novel = grepl("gene.", rowData(se)$GENEID),
                          numAstart = polyATerminals$numAstart,
                          numAend = polyATerminals$numAend,
                          numTstart = polyATerminals$numTstart,
                          numTend = polyATerminals$numTend)
  return(rowData)
}

#' % of a genes read counts assigned to each read class
calculateGeneProportion = function(counts, geneIds){
  countsTBL <- tibble(counts, geneIds) %>%
    group_by(geneIds) %>% mutate(geneReadCount = sum(counts),
                                 geneReadProp = counts/geneReadCount)
  return(countsTBL)

}

#' returns number of ref anno each read class is a subset of
isReadClassCompatible =  function(query, subject){
  outData <- data.frame(compatible=rep(0, length(query)), equal = rep(FALSE, length(query)))
  query <- cutStartEndFromGrangesList(query)
  subject <- cutStartEndFromGrangesList(subject)
  
  olap = findOverlaps(query, subject, ignore.strand = F, type = 'within')
  query <- query[queryHits(olap)]
  subject <- subject[subjectHits(olap)]
  splice <- myGaps(query)

  comp <- myCompatibleTranscription(query = query, subject = subject, splice = splice)
  equal <- elementNROWS(query)==elementNROWS(subject) & comp
  
  outData$compatible <- countQueryHits(olap[comp])
  outData$equal <- countQueryHits(olap[equal])>0
  
  return(outData)
}

#' returns number of A/T's each read class aligned 5' and 3' end
countPolyATerminals = function(grl, genomeSequence){
  start <- resize(granges(unlist(selectStartExonsFromGrangesList(grl, exonNumber = 1), 
                                 use.names = F)), width = 10, fix = 'start', ignore.strand=F)
  end <- resize(granges(unlist(selectEndExonsFromGrangesList(grl, exonNumber = 1), 
                               use.names = F)), width = 10, fix = 'end', ignore.strand=F)
  startTemp = start
  start[which(unlist(unique(strand(grl))) == '-')] = end[which(unlist(unique(strand(grl))) == '-')]
  end[which(unlist(unique(strand(grl))) == '-')] = startTemp[which(unlist(unique(strand(grl))) == '-')]
  startSeqs = BSgenome::getSeq(genomeSequence,start)
  endSeqs = BSgenome::getSeq(genomeSequence,end)
  numATstart = letterFrequency(startSeqs, c("A","T"))
  numATend= letterFrequency(endSeqs, c("A","T"))
  return(data.frame(numAstart=numATstart[,"A"], numAend= numATend[,"A"],  
                    numTstart=numATstart[,"T"], numTend=numATend[,"T"]))
}

#' calculates a score based on how likely the read class is associated with a 
#' real gene
getGeneScore = function(se){
  geneFeatures = prepareGeneModelFeatures(rowData(se))
  if(checkFeatures(geneFeatures)){
    geneModel = fit_xgb(dplyr::select(geneFeatures,!c(labels, names)),geneFeatures$labels)
    geneScore = as.numeric(predict(geneModel, as.matrix(features), 
                                   s = "lambda.min",type="response"))
    geneFDR = calculateFDR(geneScore, geneFeatures$labels)
    geneRCMap = match(rowData(se)$GENEID, geneFeatures$names)
    geneScore = geneScore[geneRCMap]
    geneFDR = geneFDR[geneRCMap]

  } else {
    message("Gene Score not calculated")
    geneScore = rep(1,nrow(se))
    geneFDR = rep(1,nrow(se))
  }
  return(data.frame(geneScore = geneScore, geneFDR = geneFDR))   
}

#' calculates the minimum FDR for each score 
calculateFDR = function(score, labels){
  scoreOrder = order(score, decreasing = T)
  orderSave = (1:length(score))[scoreOrder]
  labels = labels[scoreOrder]
  score = score[scoreOrder]
  FDR = cumsum(labels)/(1:length(score))
  FDR = rev(cummin(rev(FDR)))
  FDR = FDR[order(orderSave)]
  return(FDR)
}

#' calculate and format features by gene for model
prepareGeneModelFeatures = function(rowData){
  outData <- as_tibble(rowData) %>% group_by(GENEID) %>% 
    summarise(numReads = geneReadCount[1],
              #numReads = sum(readCount, na.rm=T), 
              names = GENEID[1],
              labels = !novel[1], 
              strand_bias = 1-abs(0.5-(sum(readCount.posStrand, na.rm=T)/numReads)), 
              numRCs=n(), 
              numExons = max(numExons, na.rm=T), 
              isSpliced = numExons>1, 
              # subsets are compatible with 2 or more RCs
              # or compat with only 1 but are not equal
              numNonSubsetRCs = numRCs - sum(compatible >= 2 |(compatible == 1 & !equal)),
              highConfidence=any(confidenceType=='highConfidenceJunctionReads')) %>%
    mutate(numReads = log2(pmax(1,numReads)))
  return(outData)
}

#' ensures that the data is trainable after filtering
checkFeatures = function(features){
  labels = features$labels
  if(sum(labels)==length(labels) | sum(labels)==0){
    message("Missing presence of both TRUE and FALSE labels.")
    return(F)
  }
  if(length(labels)<50){
    message("Not enough data points")
    return(F)
  }
  return(T)
}

#' calculates a score based on how likely a read class is full length
getTranscriptScore = function(se){
  txFeatures = prepareTranscriptModelFeatures(rowData(se))
  features = cbind(txFeatures$numReads,txFeatures$startSD, txFeatures$endSD, 
                   txFeatures$geneReadProp, txFeatures$tx_strand_bias,
                   txFeatures$numAstart, txFeatures$numAend, 
                   txFeatures$numTstart, txFeatures$numTend)
  if(checkFeatures(txFeatures)){
    transcriptModel = fit_xgb(features, 
      txFeatures$labels)
    txScore = predict(transcriptModel, as.matrix(features))

    #calculates the FDR for filtering RCs based on wanted precision
    txFDR = calculateFDR(txScore, !rowData(se)$equal)

  } else {
    message("Transcript Score not calculated")
    txScore = rep(1,nrow(se))
    txFDR = rep(1,nrow(se))
  }
  return(data.frame(txScore = txScore, txFDR = txFDR))
}

#' calculate and format read class features for model training
prepareTranscriptModelFeatures = function(rowData){
  outData <- as_tibble(rowData) %>%  
    dplyr::select(numReads = readCount, geneReadProp, startSD, endSD, numAstart, numAend, 
           numTstart,numTend, tx_strand_bias = readCount.posStrand, labels = equal) %>%
    mutate(numReads = log2(pmax(1,numReads)), 
           tx_strand_bias=(1-abs(0.5-(tx_strand_bias/numReads))))
  return(outData)
}

#' helper function to get weightedMean
applyWeightedMean = function(x){
  weighted.mean(x[1:(length(x)/2)],x[((length(x)/2)+1):length(x)], na.rm = T)
}

#train a model using xgboost, using part of the features as training set
fit_xgb = function(features, labels) {
  # Fit the xgb model
  xgb_model = xgboost(data = as.matrix(features), label = labels,
  nthread=1, nround= 50, objective = "binary:logistic", 
  eval_metric='error', verbose = 0)

  return(xgb_model)
  
}