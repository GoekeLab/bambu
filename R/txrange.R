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
  rowData(se)$compatible = isReadClassCompatible(rowRanges(se), annotations)
  rowData(se)$GENEID = assignGeneIds(rowRanges(se), annotations)
  rowData(se)$novel = grepl("gene.", 
      rowData(se)$GENEID)
  se = calculateGeneProportion(se)
  polyATerminals = countPolyATerminals(rowRanges(se), genomeSequence)
  rowData(se)$numAstart = polyATerminals$numAstart
  rowData(se)$numAend = polyATerminals$numAend
  rowData(se)$numTstart = polyATerminals$numTstart
  rowData(se)$numTend = polyATerminals$numTend
  return(se)
}

#' % of a genes read counts assigned to each read class
calculateGeneProportion = function(resultOutput){
  countsTBL <- as_tibble(rowData(resultOutput)$readCount) %>%
    mutate(geneId = rowData(resultOutput)$GENEID) %>%
    group_by(geneId) %>%
    mutate_at(vars(-geneId), .funs = sum) %>%
    ungroup() %>%
    dplyr::select(-geneId)
  geneReadProp <- rowData(resultOutput)$readCount / countsTBL
  rowData(resultOutput)$geneReadProp = unlist(geneReadProp)
  return(resultOutput)
}

#' checks to see if a read classes intron junctions fully matches an annotation
isReadClassEqual = function(query, subject){
    olapEqual = findOverlaps(cutStartEndFromGrangesList(query),
      cutStartEndFromGrangesList(subject), ignore.strand = F, type = 'equal',
      select= 'first')
    equal = !is.na(olapEqual)
    return(equal)
}

#' returns number of ref anno each read class is a subset of
isReadClassCompatible =  function(query, subject){
  compatible = rep(0, length(query))
  
  olap = findOverlaps(cutStartEndFromGrangesList(query),
                      cutStartEndFromGrangesList(subject), ignore.strand = F, type = 'within')

  query <- query[queryHits(olap)]
  subject <- subject[subjectHits(olap)]
  splice <- myGaps(query)

  qrng <- ranges(query)
  srng <- ranges(subject)
  sprng <- ranges(splice)
  
  #calculates if query is a subset of subject
  bnds <- elementNROWS(GenomicRanges::setdiff(qrng, srng)) == 0L
  splc <- elementNROWS(GenomicRanges::intersect(srng, sprng)) == 0L
  
  #count number of compatible matches
  counts = by(bnds & splc, queryHits(olap), sum)
  compatible[as.numeric(names(counts))] = counts
  
  return(compatible)
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
  return(data.frame(numAstart=numATstart[,"A"], numAend= numATend[,"A"],  numTstart=numATstart[,"T"], numTend=numATend[,"T"]))
}

#' calculates a score based on how likely the read class is associated with a 
#' real gene
getGeneScore = function(se, thresholdIndex){
  geneFeatures = prepareGeneModelFeatures(se[thresholdIndex,])
  if(checkFeatures(geneFeatures)){
  if(sum(geneFeatures$labels)==length(geneFeatures$labels) | 
    sum(geneFeatures$labels)==0){geneModel = NULL}
  geneModel = fit_xgb(geneFeatures$features,geneFeatures$labels)$cvfit
  geneScore = as.numeric(predict(geneModel, as.matrix(geneFeatures$features), 
      s = "lambda.min",type="response"))
  names(geneScore) = geneFeatures$names
  
  labels = geneFeatures$labels[order(geneScore, decreasing = T)]  
  geneScore = geneScore[order(geneScore, decreasing = T)]  
  geneFDR = cumsum(labels)/(1:length(geneScore))
  geneFDR = rev(cummin(rev(geneFDR)))
  names(geneFDR) = names(geneScore)
  
  geneFDR = calculateFDR(geneScore, geneFeatures$labels)
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
prepareGeneModelFeatures = function(se){
  temp = by(rowData(se)$readCount, rowData(se)$GENEID, sum)
  geneIDs = names(temp)
  labels = !grepl("gene.",geneIDs)
  numReads = as.numeric(temp)
  numReadsLog = log(numReads,2)
  numReadsLog[is.infinite(numReadsLog)]=0

  strand_bias = 1-abs(0.5-(as.numeric(by(rowData(se)$readCount.posStrand, 
    rowData(se)$GENEID, sum))/numReads))
  strand_bias[is.na(strand_bias)]=0
  #how many read classes does a gene have
  numRCs = table(rowData(se)$GENEID)
  
  #number of non-subset read classes
  subset = rowData(se)$compatible >= 2 | 
              (rowData(se)$compatible == 1 & 
              !rowData(se)$equal)
  subsetRCs = table(rowData(se)$GENEID[subset])
  temp = numRCs
  temp[names(subsetRCs)] = temp[names(subsetRCs)]-subsetRCs
  numNonSubsetRCs = temp
  #max number of exons of the read classes for each gene
  numExons = as.numeric(by(rowData(se)$numExons, rowData(se)$GENEID, max))
  
  #is spliced?
  isSpliced = ifelse(numExons>1,1,0)
  
  #highCOnfidence
  highConfidence = by(rowData(se)$confidenceType, rowData(se)$GENEID, 
    function(x){
      ifelse("highConfidenceJunctionReads" %in% x, 1,0)
    })
  highConfidence = as.numeric(highConfidence)

  features = cbind(numReadsLog,strand_bias, numRCs,
    numExons, isSpliced, highConfidence)

  return(list(features=features, labels = labels, names = geneIDs))
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
      txFeatures$labels[txIndex])$cvfit
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

  return (list(score = xgb_probs, cvfit = xgb_model, testData = val_data, 
    testLabels = val_labels, trainCount = nrow(train_data), 
    indicesTest = val_idx))
  
}