#' Assigns each read class geneScore and txScore
#' @param se summerized experiment object with read classes/ranges
#' @param genomeSequence genomeSequence
#' @param annotations GRangesList of annotations
txrange.scoreReadClasses = function(se, genomeSequence, annotations, 
                                    min.readCount = 2){
    options(scipen = 999)
    se = addRowData(se, genomeSequence, annotations)
    thresholdIndex = which(rowData(se)$readCount
        >=min.readCount)
    geneScore = getGeneScore(se, thresholdIndex)
    rowData(se)$geneScore = geneScore$geneScore
    rowData(se)$geneFDR = geneScore$geneFDR
    txScore = getTranscriptScore(se, thresholdIndex)
    rowData(se)$txScore = txScore$txScore
    rowData(se)$txFDR = txScore$txFDR
    
    se = se[order(unlist(unique(seqnames(rowRanges(se)))), 
      min(start(rowRanges(se)))),]
    
    return(se)
}

#' calculates labels and features used in model generation
addRowData = function(se, genomeSequence, annotations){
  rowData(se)$numExons <- elementNROWS(rowRanges(se))
  rowData(se)$equal = isReadClassEqual(rowRanges(se), annotations)
  compTable <- isReadClassCompatible(rowRanges(se), annotations)
  rowData(se)$equal = compTable$equal
  rowData(se)$compatible = compTable$compatible
  rowData(se)$GENEID = assignGeneIds(rowRanges(se), annotations)
  rowData(se)$novel = grepl("gene.", 
      rowData(se)$GENEID)
  countsTBL = calculateGeneProportion(counts=mcols(se)$readCount,
                          geneIds=mcols(se)$GENEID)
  rowData(se)$geneReadProp = countsTBL$geneReadProp
  rowData(se)$geneReadCount = countsTBL$geneReadCount
  polyATerminals = countPolyATerminals(rowRanges(se), genomeSequence)
  rowData(se)$numAstart = polyATerminals$numAstart
  rowData(se)$numAend = polyATerminals$numAend
  rowData(se)$numTstart = polyATerminals$numTstart
  rowData(se)$numTend = polyATerminals$numTend
  return(se)
}

#' % of a genes read counts assigned to each read class
calculateGeneProportion = function(counts, geneIds){
  countsTBL <- tibble(counts, geneIds) %>%
    group_by(geneIds) %>% mutate(geneReadProp = counts/sum(counts), 
                                 geneReadCount = sum(counts))
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

  counts <- countQueryHits(olap[comp])

  outData$compatible <- counts
  outData$equal <- countQueryHits(olap[equal])>0
  
  return(outData)
}

#' returns number of A/T's each read class aligned 5' and 3' end
countPolyATerminals = function(grl, genomeSequence){
  start <- resize(granges(unlist(selectStartExonsFromGrangesList(grl, exonNumber = 1), use.names = F)), width = 10, fix = 'start', ignore.strand=F)
  end <- resize(granges(unlist(selectEndExonsFromGrangesList(grl, exonNumber = 1), use.names = F)), width = 10, fix = 'end', ignore.strand=F)
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
getGeneScore = function(se, thresholdIndex){
  geneFeatures = prepareGeneModelFeatures(rowData(se)[thresholdIndex,])
  labels = geneFeatures$labels
  features = cbind(geneFeatures$numReads, geneFeatures$strand_bias, 
                   geneFeatures$numRCs, geneFeatures$numExons, 
                   geneFeatures$isSpliced , geneFeatures$highConfidence)
  if(checkFeatures(geneFeatures)){
  if(sum(labels)==length(labels) | 
    sum(labels)==0){geneModel = NULL}
  geneModel = fit_xgb(features,labels)
  geneScore = as.numeric(predict(geneModel, as.matrix(features), 
      s = "lambda.min",type="response"))
  names(geneScore) = geneFeatures$names
  
  labels = labels[order(geneScore, decreasing = T)]  
  geneScore = geneScore[order(geneScore, decreasing = T)]  
  geneFDR = cumsum(labels)/(1:length(geneScore))
  geneFDR = rev(cummin(rev(geneFDR)))
  names(geneFDR) = names(geneScore)
  
  geneFDR = calculateFDR(geneScore, labels)
  names(geneFDR) = names(geneScore)

  geneScore = geneScore[rowData(se)$GENEID]
  geneFDR = geneFDR[rowData(se)$GENEID]

  } else {
    message("Gene Score not calculated")
    geneScore = rep(1,nrow(se))
    geneFDR = rep(1,nrow(se))
  }
  return(list(geneScore = geneScore, geneFDR = geneFDR))  
}

#' calculates the minimum FDR for each score 
calculateFDR = function(score, labels){
  scoreOrder = order(score, decreasing = T)
  orderSave = (1:length(score))[scoreOrder]
  labels = labels[scoreOrder]
  score = score[scoreOrder]
  FDR = cumsum(labels)/(1:length(score))
  FDR = rev(cummin(rev(FDR)))
  FDR = FDR[orderSave]
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
              numNonSubsetRCs = numRCs - sum(compatible >= 2 |(compatible == 1 & equal)),
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
getTranscriptScore = function(se, thresholdIndex){
  txFeatures = prepareTranscriptModelFeatures(se)
  if(checkFeatures(txFeatures)){
    txIndex = thresholdIndex[thresholdIndex %in% 
      which(!rowData(se)$novel)]
    transcriptModel = fit_xgb(txFeatures$features[txIndex,], 
      txFeatures$labels[txIndex])
    txScore = predict(transcriptModel, as.matrix(txFeatures$features),
      s = "lambda.min", type="response")

    #calculates the FDR for filtering RCs based on wanted precision
    txFDR = calculateFDR(txScore, !rowData(se)$equal)

  } else {
    message("Transcript Score not calculated")
    txScore = rep(1,nrow(se))
    txFDR = rep(1,nrow(se))
  }
  return(list(txScore = txScore, txFDR = txFDR))
}

#' calculate and format read class features for model training
prepareTranscriptModelFeatures = function(input){
  labels = rowData(input)$equal
  
  numReads = rowData(input)$readCount
  logNumReads = log(numReads, 2)
  geneReadProp=rowData(input)$geneReadProp
  geneReadProp[is.na(geneReadProp)]=0
  tx_strand_bias=(1-abs(0.5-(rowData(input)$readCount.posStrand/numReads)))
  SD = apply(cbind(rowData(input)$startSD, rowData(input)$readCount),MARGIN = 1, 
    FUN = applyWeightedMean)*-1
  SD[which(is.na(SD))] = 1
  SDend = apply(cbind(rowData(input)$endSD, rowData(input)$readCount),1,
    FUN = applyWeightedMean)*-1
  SDend[which(is.na(SDend))] = 1
  numAstart = rowData(input)$numAstart
  numAstart[is.infinite(numAstart)]=0
  numAend = rowData(input)$numAend
  numAend[is.infinite(numAend)]=0
  numTstart = rowData(input)$numTstart
  numTstart[is.infinite(numTstart)]=0
  numTend = rowData(input)$numTend
  numTend[is.infinite(numTend)]=0
  
  features = cbind(numReads, SD, SDend, geneReadProp, tx_strand_bias,
    numAstart, numAend, numTstart, numTend)
  return(list(features = features, labels = labels))
}

#' helper function to get weightedMean
applyWeightedMean = function(x){
  weighted.mean(x[1:(length(x)/2)],x[((length(x)/2)+1):length(x)], na.rm = T)
}

#train a model using xgboost, using part of the features as training set
fit_xgb = function(features, labels) {

  # Sample the data into train, val and test sets
  train_idx = sample(nrow(features), floor(0.9*nrow(features)))
  val_idx = setdiff(seq_len(nrow(features)),train_idx)
  
  train_data= features[train_idx,]
  val_data = features[val_idx,]
  train_labels = labels[train_idx]
  val_labels = labels[val_idx]
  
  x_mat_train = as.matrix(train_data)
  x_mat_val = as.matrix(val_data)

  # Fit the xgb model
  xgb_model = xgboost(data = x_mat_train, label = train_labels,
  nthread=1, nround= 50, objective = "binary:logistic", 
  eval_metric='error', verbose = 0)
  xgb_probs = predict(xgb_model, x_mat_val)

  return(xgb_model)
  
}