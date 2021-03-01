txrange.scoreReadClasses = function(se, genomeSequence, annotations, 
  withAdapters = FALSE, min.readCount = 2){
    options(scipen = 999)
    se = addRowData(se, genomeSequence, annotations)
    thresholdIndex = which(rowSums(rowData(se)$counts)
        >=min.readCount)
    rowData(se)$geneScore = getGeneScore(se, thresholdIndex, 
      method = "xgboost")
    rowData(se)$txScore = getTranscriptScore(se, thresholdIndex, 
      method = "xgboost")

    return(se)
}

addRowData = function(se, genomeSequence, annotations){
  exons = str_split(rowData(se)$intronStarts,",")
  rowData(se)$numExons = sapply(exons, FUN = length)+1
  rowData(se)$numExons[is.na(exons)] = 1
  rowData(se)$equal = isReadClassEqual(rowRanges(se), annotations)
  rowData(se)$GENEID = assignGeneIds(rowRanges(se), annotations)
  rowData(se)$novel = grepl("gene.", 
      rowData(se)$GENEID)
  se = calculateGeneProportion(se)
  se = countPolyATerminals(se, genomeSequence)
 
  return(se)
}

calculateGeneProportion = function(resultOutput){
  countsTBL <- as_tibble(rowData(resultOutput)$counts) %>%
    mutate(geneId = rowData(resultOutput)$GENEID) %>%
    group_by(geneId) %>%
    mutate_at(vars(-geneId), .funs = sum) %>%
    ungroup() %>%
    dplyr::select(-geneId)
  geneReadProp <- rowData(resultOutput)$counts / countsTBL
  rowData(resultOutput, withDimnames = F)$geneReadProp = geneReadProp
  rowData(resultOutput)$totalGeneReadProp = 
    rowSums(rowData(resultOutput)$counts / rowSums(countsTBL))
  return(resultOutput)
}

isReadClassEqual = function(query, subject, maxDist = 5){
    olapEqual = findOverlaps(cutStartEndFromGrangesList(query),
      cutStartEndFromGrangesList(subject), ignore.strand = F, type = 'equal',
      select= 'first')
    equal = !is.na(olapEqual)
    return(equal)
}

countPolyATerminals = function(se, genomeSequence){
  #counts A/T's at 5' and 3' of RCs on genome
  #get all first/last exons

  start = min(start(rowRanges(se)))
  end = max(end(rowRanges(se)))
  strand = unlist(runValue(strand(rowRanges(se))))
  seqname = unlist(runValue(seqnames(rowRanges(se))))
  df = data.frame(start,end,strand,seqname)
  RCranges = makeGRangesFromDataFrame(df)
  strand(RCranges)[which(as.character(strand(RCranges))=='*')]='+'
  genomeSequence = checkInputSequence(genomeSequence)
  
  startSeqs = as.character(BSgenome::getSeq(genomeSequence,resize(RCranges,10,
    fix="start")[which(as.character(seqnames(RCranges)) %in%
    names(genomeSequence))]))
  endSeqs = as.character(BSgenome::getSeq(genomeSequence,resize(RCranges,10,
    fix="end")[which(as.character(seqnames(RCranges)) 
    %in% names(genomeSequence))]))
  #count number of A's in the first/last 10 bp
  index = which(as.character(seqnames(RCranges)) %in% names(genomeSequence))
  numAstart = rep(NA,nrow(se))  
  numAstart[index] = str_count(startSeqs, "A")
  numAend = rep(NA,nrow(se)) 
  numAend[index] = str_count(endSeqs, "A")
  
  numTstart = rep(NA,nrow(se))  
  numTstart[index] = str_count(startSeqs, "T")
  numTend = rep(NA,nrow(se)) 
  numTend[index] = str_count(endSeqs, "T")
  
  rowData(se)$numAstart = numAstart
  rowData(se)$numAend = numAend
  rowData(se)$numTstart = numTstart
  rowData(se)$numTend = numTend
  return(se)
}

prepareTranscriptModelFeatures = function(input, withAdapters = F){

  labels = rowData(input)$equal
  
  features = getAgnosticFeatures(input)
  
  return(list(features = features, labels = labels))
}

applyWeightedMean = function(x){
  weighted.mean(x[1:(length(x)/2)],x[((length(x)/2)+1):length(x)], na.rm = T)
}

getAgnosticFeatures = function(input){
  numReads = rowSums(rowData(input)$counts)
  logNumReads = log(numReads, 2)
  geneReadProp=rowData(input)$totalGeneReadProp
  geneReadProp[is.na(geneReadProp)]=0
  tx_strand_bias=(1-abs(0.5-(rowSums(rowData(input)$strand_bias)/numReads)))
  SD = apply(cbind(rowData(input)$startSD, rowData(input)$counts),MARGIN = 1, 
    FUN = applyWeightedMean)*-1
  SD[which(is.na(SD))] = 1
  SDend = apply(cbind(rowData(input)$endSD, rowData(input)$counts),1,
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
  transcriptProp = as.numeric(rowData(input)$transcriptProp)
  transcriptProp[is.na(transcriptProp)] = 0
  transcriptPropMin = as.numeric(rowData(input)$transcriptPropMin)
  transcriptPropMin[is.na(transcriptPropMin)] = 0

  features = cbind(numReads, SD, SDend, geneReadProp, tx_strand_bias,
    numAstart, numAend, numTstart, numTend, transcriptProp,
    transcriptPropMin)
  
  return(features)
}

prepareGeneModelFeatures = function(se){
  #group read classes by gene
  #summerize the features
  temp = by(rowData(se)$counts, rowData(se)$GENEID, sum)
  geneIDs = names(temp)
  labels = !grepl("gene.",geneIDs)
  numReads = as.numeric(temp)
  numReadsLog = log(numReads,2)
  numReadsLog[is.infinite(numReadsLog)]=0

  strand_bias = 1-abs(0.5-(as.numeric(by(rowData(se)$strand_bias, 
    rowData(se)$GENEID, sum))/numReads))
  strand_bias[is.na(strand_bias)]=0
  #how many read classes does a gene have
  numRCs = table(rowData(se)$GENEID)
  
  #number of non-subset read classes
  subset = (sapply(str_split(rowData(se)$resultOutput.compatible, ';'),
    length)>=2 | (sapply(str_split(rowData(se)$resultOutput.compatible, ';'),
    length)==1 & rowData(se)$resultOutput.equal == "NOVEL" & 
    rowData(se)$resultOutput.compatible != "non_compatible_match"))
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

  features = cbind(numReadsLog,strand_bias, numRCs, numNonSubsetRCs,
    numExons, isSpliced, highConfidence)

  return(list(features=features, labels = labels, names = geneIDs))
}

getGeneScore = function(se, thresholdIndex, method = "sgboost"){
  geneFeatures = prepareGeneModelFeatures(se[thresholdIndex,])
    if(checkFeatures(geneFeatures)){
    geneModel = trainGeneModel(geneFeatures$features, 
      geneFeatures$labels, geneFeatures$names, method)
    geneScore = calculateGeneScore(geneFeatures$features, geneFeatures$labels, 
      geneFeatures$names, model = geneModel, method = method)$score
    rowData(se)$geneScore = geneScore[rowData(se)$GENEID]
    # seSorted = seTrimmed[order(geneScore, decreasing = T),]  
    # rowData(seSorted)$FDR = cumsum(!grepl("gene.",
    #   rowData(seSorted)$GENEID))/(1:nrow(rowData(seSorted)))
    } else {
      message("Gene Score not calculated")
      rowData(se)$geneScore = rep(1,nrow(se))
    }
}

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

trainGeneModel = function(features, labels, names, 
  method = "sgboost"){
  if(sum(labels)==length(labels) | sum(labels)==0){return(NULL)}
  geneScore = calculateGeneScore(features, labels, names, method = method)
  return(geneScore$model$cvfit)
}

calculateGeneScore = function(features, labels, names, model = NULL, 
  method = "xgboost"){
  if(!is.null(model)){
    score = as.numeric(predict(model, as.matrix(features), 
      s = "lambda.min",type="response"))
  } else{
    result = calculateScore(features,labels, method = method)
    model = result
    score = as.numeric(predict(result$cvfit, as.matrix(features), 
      s = "lambda.min",type="response"))
  }
  names(score) = names
  
  return(list(score=score, model = model))
}

getTranscriptScore = function(se, thresholdIndex, 
  method = "xgboost"){
  txFeatures = prepareTranscriptModelFeatures(se,
    withAdapters = withAdapters)
  if(checkFeatures(txFeatures)){
    txIndex = thresholdIndex[thresholdIndex %in% 
      which(!rowData(se)$novel)]
    #txIndex = thresholdIndex
    transcriptModel = trainTranscriptModel(txFeatures$features[txIndex,], 
      txFeatures$labels[txIndex], method = method)
    transcriptScore = predict(transcriptModel, txFeatures$features,
      s = "lambda.min", type="response")
    rowData(se)$transcriptScore = transcriptScore
    # seSorted = seFrac[order(transcriptScore, decreasing = T),]
    # rowData(seSorted)$FDR = cumsum(rowData(seSorted)$equal != "NOVEL")/
    # (1:nrow(rowData(seSorted)))
  } else {
    message("Transcript Score not calculated")
    rowData(se)$transcriptScore = rep(1,nrow(se))
  }
}

trainTranscriptModel = function(features, labels, 
  method = "xgboost"){
  result=calculateScore(features, labels, method = method)
  return(result$cvfit)
}

calculateScore = function(input, labels, model = NULL, log = F, 
  method = "xgboost"){
  if(method == "xgboost"){
    return(fit_xgb(input, labels))
  } else{
    return(fit_glmnet(input,labels, log = log))
  }
}

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
  negative_labels = sum(train_labels == 0)
  positive_labels = sum(train_labels == 1)
  xgb_time = system.time({xgb_model = xgboost(data = x_mat_train, 
  label = train_labels, nthread=1, max.depth = 3, nround= 100, 
  objective = "binary:logistic", 
  scale_pos_weight=negative_labels/positive_labels, verbose = 0)})
  xgb_probs = predict(xgb_model, x_mat_val)
  return (list(score = xgb_probs, cvfit = xgb_model, testData = val_data, 
    testLabels = val_labels, trainCount = nrow(train_data), 
    indicesTest = val_idx))
  
}

fit_glmnet = function(input, labels, model = NULL, log = F) {
  input = cbind(input,labels)
  finalFeatureIndex = ncol(input)-1
  sampleSizeTraining=0
  sampleSizeValidation=0
  if(is.null(model)){     #split data for training
    fractionTraining   <- 0.60
    fractionValidation <- 0.20
    fractionTest       <- 0.20
    sampleSizeTraining   <- floor(fractionTraining   * nrow(input))
    sampleSizeValidation <- floor(fractionValidation * nrow(input))
    sampleSizeTest       <- floor(fractionTest       * nrow(input))
    indicesTraining    <- sort(sample(seq_len(nrow(input)), 
      size=sampleSizeTraining))
    indicesNotTraining <- setdiff(seq_len(nrow(input)), indicesTraining)
    indicesValidation  <- sort(sample(indicesNotTraining, 
      size=sampleSizeValidation))
    indicesTest        <- setdiff(indicesNotTraining, indicesValidation)
    dfTraining   <- input[indicesTraining, ]
    dfValidation <- input[indicesValidation, ]
    dfTest       <- input[indicesTest, ]
    if(!log){
      cvfit = cv.glmnet(as.matrix(rbind(dfTraining[,1:finalFeatureIndex],
        dfValidation[,1:finalFeatureIndex])), 
        c(dfTraining[,finalFeatureIndex+1],dfValidation[,finalFeatureIndex+1]),
        family = "binomial", type.measure = "class")
    } else {
      cvfit = glmnet(as.matrix(rbind(dfTraining[,1:finalFeatureIndex],
        dfValidation[,1:finalFeatureIndex])), c(dfTraining[,finalFeatureIndex+1],
        dfValidation[,finalFeatureIndex+1]),family = "binomial", 
        type.measure = "class")
      }
    } else{
    dfTest = input
    cvfit = model
  }
  if(!log){
    score = as.numeric(predict(cvfit, 
      newx = as.matrix(dfTest[,1:finalFeatureIndex]), 
      s="lambda.min",type="response"))
  } else {
    score = as.numeric(predict(cvfit, 
      newx = as.matrix(dfTest[,1:finalFeatureIndex]),
      s = 0.01, type="response"))
  }
  testLabels = dfTest[,finalFeatureIndex+1]
  testData = dfTest[,1:finalFeatureIndex]
  return(list(score=score, cvfit=cvfit, testData = testData, 
    testLabels=testLabels, 
    trainCount = sampleSizeTraining + sampleSizeValidation, 
    indicesTest = indicesTest))
}
