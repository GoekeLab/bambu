suppressMessages(require(stringr))
library(BSgenome)
library(ROCit)
library(glmnet)
library(xgboost)

if(F){
  load("se_bambuTest.Rdata")
  load("readGrgList_bambuTest.Rdata")
  min.readCount = 2
  annotations = bambuAnnotations
  genomeSequence = fa.file
}

txrange.filterReadClasses = function(se, readGrgList, genomeSequence,
      annotations, withAdapters = FALSE, min.readCount = 2){

    options(scipen = 999)
    #alignData = createAlignData(readGrgList)
    #alignData = annotateReadStartsAndEnds(alignData, se)
    rm(readGrgList)
    se = combineSEs(list(se), annotations)
    se = addRowData(se, genomeSequence, annotations)
    thresholdIndex = which(rowSums(assays(se)$counts)
        >=min.readCount)
    rowData(se)$geneScore = getGeneScore(se, thresholdIndex, 
      plot = NULL, method = "xgboost")
    rowData(se)$txScore = getTranscriptScore(se, thresholdIndex, 
      plot = NULL, method = "xgboost")
    
    # test model on subset and nonsubset RCs
    # subset = rowData(se)$compatibleCount >= 2 | 
    #  (!rowData(se)$equal & rowData(se)$compatibleCount == 1)
    # getTranscriptScore(se[which(!subset),], thresholdIndex, 
    #                    plot = "NULL", method = "xgboost")
    # getTranscriptScore(se[which(subset),], thresholdIndex, 
    #                    plot = "NULL", method = "xgboost")
    return(se)
}

#deprecated
createAlignData = function(readGrgList){
  #converts the readGrgList into a simple matrix
  #todo, instead of this function just use mcols
  seqnames = unlist(runValue(seqnames(readGrgList)))
  strand = unlist(runValue(strand(readGrgList)))
  start = min(start(readGrgList))
  end = max(end(readGrgList))
  qwidth = sum(width(readGrgList))
  readClass = mcols(readGrgList)$readClass
  readClassIndex = mcols(readGrgList)$readClassIndex
  alignData = as.data.frame(cbind(seqnames,strand,qwidth, start,end, 
      readClass, readClassIndex))
  return(alignData)
}

#deprecated
annotateReadStartsAndEnds = function(alignData, se){
  #assign labels to reads
  #rownames(allReadClasses)=allReadClasses$`readClasses$readClassId`
  #RCstrand = allReadClasses[alignData$readClass,]$strand
  RCstrand = rowData(se)[alignData$readClass,]$strand.rc
  starts = alignData$start
  starts[which(RCstrand == '-' )]=alignData[which(RCstrand == '-' ),]$end
  #TODO deal with single exon reads which are assigned *
  starts = as.numeric(starts)
  ends = alignData$end
  ends[which(RCstrand == '-' )]=alignData[which(RCstrand == '-' ),]$start
  ends = as.numeric(ends)
  
  txstrand = alignData$strand
  revindex = which(RCstrand == '-' & txstrand == '-')
  txstrand[which(RCstrand == '-' & txstrand == '+')] = '-'
  txstrand[revindex] = '+'
  txstrand[which(RCstrand == "*")] = "*"
  
  alignData$tssStart=as.numeric(starts)
  alignData$tesEnd=as.numeric(ends)
  alignData$txstrand = txstrand
  return(alignData)
}

#deprecated
simplifyAdapterData = function(alignData){
  # NOT USED WHEN NO ADAPTERS IN DATA
  alignData$adap5 = grepl('5', alignData$predictedStrand, fixed=TRUE)
  alignData$adap3 = grepl('3', alignData$predictedStrand, fixed=TRUE)

  alignData$adap5pos = grepl('5+', alignData$predictedStrand, fixed=TRUE)
  alignData$adap3pos = grepl('3+', alignData$predictedStrand, fixed=TRUE)
  alignData$adap5neg = grepl('5-', alignData$predictedStrand, fixed=TRUE)
  alignData$adap3neg = grepl('3-', alignData$predictedStrand, fixed=TRUE)

  alignData$adapPos = (alignData$adap5pos & alignData$adap3pos & 
    alignData$strand=="+") | 
    (alignData$adap5neg & alignData$adap3neg & alignData$strand=="-")
  alignData$adapNeg = (alignData$adap5neg & alignData$adap3neg & 
    alignData$strand=="+") | 
    (alignData$adap5pos & alignData$adap3pos & alignData$strand=="-")
  
  alignData$adapDist5 = as.numeric(alignData$barcode5Prime.minus.start)
  pos5 = which(grepl('5+',alignData$predictedStrand, fixed=TRUE))
  alignData$adapDist5[pos5] = as.numeric(alignData$widthRead5prime[pos5])-
    as.numeric(alignData$barcode5Prime.plus.end[pos5])
  alignData$adapDist5[!alignData$adap5] = NA
  
  alignData$adapDist3 = as.numeric(alignData$barcode3Prime.minus.start)
  pos3 = which(grepl('3-',alignData$predictedStrand, fixed=TRUE))
  alignData$adapDist3[pos3] = as.numeric(alignData$widthRead3prime[pos3])-
    as.numeric(alignData$barcode3Prime.plus.end[pos3])
  alignData$adapDist3[!alignData$adap3] = NA
  
  #polyA data
  #alignData$polyApos = rep(0,nrow(alignData))
  alignData$polyApos = as.numeric(alignData$polyA3Prime.plus.score) *
    100/as.numeric(alignData$polyA3Prime.plus.pid)
  alignData$polyApos[which(alignData$strand=="-")] = 
    as.numeric(alignData$polyA3Prime.minus.score)[which(alignData$strand=="-")] * 
    100/as.numeric(alignData$polyA3Prime.minus.pid)[which(alignData$strand=="-")]
  alignData$polyAneg = as.numeric(alignData$polyA3Prime.minus.score) * 
    100/as.numeric(alignData$polyA3Prime.minus.pid)
  alignData$polyAneg[which(alignData$strand=="-")] = 
    as.numeric(alignData$polyA3Prime.plus.score)[which(alignData$strand=="-")] * 
    100/as.numeric(alignData$polyA3Prime.plus.pid)[which(alignData$strand=="-")]
  
  #if polyA is present opposide side of sequencing start site
  alignData$polyAFull = alignData$txstrand=="+" & alignData$polyApos > 10 

   return(alignData)
}

#deprecated
summeriseReadsByGroup <- function(alignData, withAdapters = FALSE){
  
  alignData=alignData[which(!is.na(alignData$readClass)),]

  if(withAdapters){
    result <- alignData %>%
      dplyr::select(readClassIndex, readClass, strand, txstrand, uniqueRead, 
        tssStart, tesEnd, adap5, adap3, adapPos, adapNeg, adapDist5, 
        adapDist3, polyApos, polyAneg) %>%
      group_by(readClassIndex) %>%
      summarise(readClass = readClass[1],
                counts = n(),
                strand_bias = sum(txstrand == '+'),
                startSD = sd(tssStart),
                endSD = sd(tesEnd),
                uniqueReads = sum(uniqueRead),
                adapNumReads = sum(adap5),
                adapEndNumReads = sum(adap3),
                adaptDist = median(adapDist5,na.rm = T),
                adapEndDist = median(adapDist3,na.rm = T),
                bothAdapters = sum(adapPos | adapNeg),
                polyA = max(median(polyApos), median(polyAneg)),
                polyAEnd = sum(polyApos >= 10 & txstrand == '+', na.rm = T)
      )
  } else{
    result <- alignData %>%
      dplyr::select(readClassIndex, readClass, strand, txstrand, 
        tssStart, tesEnd) %>%
      group_by(readClassIndex) %>%
      summarise(readClass = readClass[1],
                counts = n(),
                strand_bias = sum(txstrand == '+'),
                startSD = sd(tssStart),
                endSD = sd(tesEnd),
      )
  }
  rownames(result)=result$readClass
  return(result)
}

combineSEs = function(combinedOutputs, annotations){
  combinedTxCandidates = NULL
  for(i in 1:length(combinedOutputs)){
    combinedTxCandidates <- 
      isore.combineTranscriptCandidates(combinedOutputs[[i]], readClassSeRef =
      combinedTxCandidates)
  }
  colnames(rowData(combinedTxCandidates)) = c(c("chr.rc", "start.rc",
     "end.rc", "strand.rc"), 
    colnames(rowData(combinedOutputs[[1]])[-1:-4]))
  return(combinedTxCandidates)
}


addRowData = function(se, genomeSequence, annotations){
  exons = str_split(rowData(se)$intronStarts,",")
  rowData(se)$numExons = sapply(exons, FUN = length)+1
  rowData(se)$numExons[is.na(exons)] = 1
  readClassesList = convertSEtoGRangesList(se)
  classifications = getReadClassClassifications(readClassesList, annotations)
  rowData(se)$equal = classifications$equal
  rowData(se)$compatible = classifications$compatible
  rowData(se)$compatibleCount = classifications$compatibleCount
  rowData(se)$GENEID = assignGeneIds(readClassesList, annotations)
  rowData(se)$novel = grepl("gene.", 
      rowData(se)$GENEID)
  se = calculateGeneProportion(se)
  #se = getTranscriptProp(se, readClassesList)
  se = countPolyATerminals(se, genomeSequence)
 
  return(se)
}

calculateGeneProportion = function(resultOutput){
  countsTBL <- as_tibble(assays(resultOutput)$counts) %>%
    mutate(geneId = rowData(resultOutput)$GENEID) %>%
    group_by(geneId) %>%
    mutate_at(vars(-geneId), .funs = sum) %>%
    ungroup() %>%
    dplyr::select(-geneId)
  geneReadProp <- assays(resultOutput)$counts / countsTBL
  assays(resultOutput, withDimnames = F)$geneReadProp = geneReadProp
  rowData(resultOutput)$totalGeneReadProp = 
    rowSums(assays(resultOutput)$counts / rowSums(countsTBL))
  return(resultOutput)
}

convertSEtoGRangesList = function(combinedOutputs){
  
  exonEndsShifted <- paste(rowData(combinedOutputs)$intronStarts, 
    rowData(combinedOutputs)$end.rc + 1, sep = ',')
  exonStartsShifted <- paste(rowData(combinedOutputs)$start.rc - 1, 
    rowData(combinedOutputs)$intronEnds, sep = ',')
  #deal with single exon RCs
  singleExonIndex = which(is.na(rowData(combinedOutputs)$intronStarts))
  exonEndsShifted[singleExonIndex] = 
    rowData(combinedOutputs)$end.rc[singleExonIndex] + 1

  exonStartsShifted[singleExonIndex] = 
    rowData(combinedOutputs)$start.rc[singleExonIndex] - 1

  readClassList= makeGRangesListFromFeatureFragments(
    seqnames=rowData(combinedOutputs)$chr.rc,
    fragmentStarts=as.character(exonStartsShifted),
    fragmentEnds=as.character(exonEndsShifted),
    strand=rowData(combinedOutputs)$strand.rc)

   # correct junction to exon differences in coordinates
  readClassList <- narrow(readClassList, start = 2, end = -2) 
  
  unlistData <- unlist(readClassList, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(elementNROWS(readClassList)), 
    names = NULL)
  
  exon_rank <- sapply(width((partitioning)), seq, from = 1)
  exon_rank[which(rowData(combinedOutputs)$strand.rc == '-')] <- 
  # * assumes positive for exon ranking
  lapply(exon_rank[which(rowData(combinedOutputs)$strand.rc == '-')], rev)  
  exon_endRank <- lapply(exon_rank, rev)
  unlistData$exon_rank <- unlist(exon_rank)
  unlistData$exon_endRank <- unlist(exon_endRank)
  
  readClassList <- relist(unlistData, partitioning)
  names(readClassList) = 1:length(readClassList)
  
  return(readClassList)
}

getReadClassClassifications = function(query, subject, maxDist = 5){

    subjectExtend <- extendGrangesListElements(subject, by = maxDist)
    queryForOverlap <- dropGrangesListElementsByWidth(query,
        minWidth = maxDist, cutStartEnd = T)
    olap = findOverlaps(queryForOverlap, subjectExtend, ignore.strand = F,
        type = 'within')
    compatible = rep(F, length(query))
    compatible[queryHits(olap)] = T
    compatibleCount = rep(0, length(query))
    compatTable = table(queryHits(olap))
    compatibleCount[as.numeric(names(compatTable))] = compatTable
    
    olapEqual = findOverlaps(queryForOverlap, 
      cutStartEndFromGrangesList(subject), ignore.strand = F, type = 'equal')
    equal = rep(F, length(query))
    equal[queryHits(olapEqual)] = T
    
    return(list(equal = equal, compatible = compatible, 
      compatibleCount = compatibleCount))
}

getTranscriptProp = function(se, readClassesList){
  # calculates the min and max contributions a read class makes to all
  i= 0
  allOverlaps = list()
  readClassesListToExtract = readClassesList
  for(x in 1:length(seqlevels(readClassesList))){
    chr = seqlevels(readClassesList)[x]
    selection = which(unlist(runValue(seqnames(readClassesListToExtract))==chr))
    readClassListTemp = readClassesListToExtract[selection]
    readClassesListToExtract= readClassesListToExtract[-selection]
    query <- cutStartEndFromGrangesList(readClassListTemp)
    overlaps = countOverlaps(query, readClassListTemp, type = 'within')
    megaRCs = readClassesList[which(overlaps == 1),]
    overlaps = findOverlaps(query, megaRCs, type = 'within')
    #account for the index change due to seperating the readclasses by chr
    overlaps = as.matrix(overlaps)+i
    i = i + length(readClassListTemp)
    #allOverlaps = rbind(allOverlaps, overlaps)
    allOverlaps[[x]]=overlaps
  }
  allOverlaps = do.call(rbind, allOverlaps)
  allCompats = allOverlaps[,'subjectHits']
  compatIndex = allOverlaps[,'queryHits']
  compatCounts = rowSums(assays(se)$counts)[allOverlaps[,'queryHits']]
  
  #inverse the data so read classes are grouped now by transcript
  annotationSubsets = by(compatIndex,allCompats, FUN = c)
  transcriptSums = by(cbind(compatCounts,compatIndex),allCompats, 
    FUN = function(x){
      counts = as.numeric(x[,1])
      transcriptProp = counts/sum(counts)
      names(transcriptProp) =  x[,2]
      return(transcriptProp)
  })
  #align the proportions with the read class names
  transcriptProp = unlist(transcriptSums)
  indexes = unlist(sapply(transcriptSums, FUN = names))
  x = as.data.frame(cbind(transcriptProp, indexes))

  topTranscriptProp_RC = x %>% group_by(indexes) %>% 
    top_n(1, transcriptProp) %>% distinct(indexes, transcriptProp, 
    .keep_all = TRUE) 
  bottomTranscriptProp_RC = x %>% group_by(indexes) %>% 
    top_n(-1, transcriptProp) %>% 
    distinct(indexes, transcriptProp, .keep_all = TRUE) 

  rowData(se)$transcriptProp = rep(0, nrow(se))
  rowData(se)$transcriptPropMin = rep(0, nrow(se))
  rowData(se)[as.numeric(topTranscriptProp_RC$indexes),]$transcriptProp = 
    topTranscriptProp_RC$transcriptProp
  rowData(se)[as.numeric(bottomTranscriptProp_RC$indexes),]$transcriptPropMin =
    bottomTranscriptProp_RC$transcriptProp

  return(se)
}

countPolyATerminals = function(se, genomeSequence){
  #counts A/T's at 5' and 3' of RCs on genome
  #get all first/last exons
  RCranges = makeGRangesFromDataFrame(rowData(se),start.field = "start.rc", 
    end.field = "end.rc", seqnames.field = "chr.rc",
    strand.field = "strand.rc")
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
  
  if("adapNumReads" %in% names(assays(input))){
    features(cbind(features,getAdapterFeatures(input)))
  }
  
  return(list(features = features, labels = labels))
}

applyWeightedMean = function(x){
  weighted.mean(x[1:(length(x)/2)],x[((length(x)/2)+1):length(x)], na.rm = T)
}

getAgnosticFeatures = function(input){
  numReads = rowSums(assays(input)$counts)
  logNumReads = log(numReads, 2)
  geneReadProp=rowData(input)$totalGeneReadProp
  geneReadProp[is.na(geneReadProp)]=0
  tx_strand_bias=(1-abs(0.5-(rowSums(assays(input)$strand_bias)/numReads)))
  SD = apply(cbind(assays(input)$startSD, assays(input)$counts),MARGIN = 1, 
    FUN = applyWeightedMean)*-1
  SD[which(is.na(SD))] = 1
  SDend = apply(cbind(assays(input)$endSD, assays(input)$counts),1,
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

getAdapterFeatures = function(input){
adapter_reads=log(rowSums(assays(input)$adapNumReads),2)
  adapter_reads[is.na(adapter_reads)]=0
  
  adapterEnd_reads=log(rowSums(assays(input)$adapEndNumReads),2)
  adapterEnd_reads[is.na(adapterEnd_reads)]=0
  
  adapter_prop = adapter_reads/numReads
  adapter_prop[is.na(adapter_prop)]=0
  adapter_prop[is.infinite(adapter_prop)]=0
  
  adapterend_prop = adapterEnd_reads/numReads
  adapterend_prop[is.na(adapterend_prop)]=0
  adapterend_prop[is.infinite(adapterend_prop)]=0
  
  adapter_distance = apply(cbind(assays(input)$adaptDist, 
    assays(input)$counts),1,FUN = applyWeightedMean)
  adapter_distance[is.na(adapter_distance)]=100
  adapter_distance[adapter_distance<=40]=1
  adapter_distance[adapter_distance>40]=0
  
  adapterend_distance = apply(cbind(assays(input)$adaptEndDist, 
    assays(input)$counts),1,FUN = applyWeightedMean)
  adapterend_distance[is.na(adapterend_distance)]=100
  adapterend_distance[adapterend_distance<=30]=1
  adapterend_distance[adapterend_distance>30]=0
  
  bothAdapters = rowSums(assays(input)$bothAdapters)
  
  polyAEnd = rowSums(assays(input)$polyAEnd)
  
  bothAdaptersProp = rowSums(assays(input)$bothAdapters)/numReads
  
  polyAEndProp = rowSums(assays(input)$polyAEnd)/numReads
  
  features = cbind(adapter_prop, adapterend_prop, adapter_distance, 
    adapterend_distance, bothAdapters, polyAEnd, bothAdaptersProp, 
    polyAEndProp)

  return(features)
}

prepareGeneModelFeatures = function(se){
  #group read classes by gene
  #summerize the features
  temp = by(assays(se)$counts, rowData(se)$GENEID, sum)
  geneIDs = names(temp)
  labels = !grepl("gene.",geneIDs)
  numReads = as.numeric(temp)
  numReadsLog = log(numReads,2)
  numReadsLog[is.infinite(numReadsLog)]=0

  strand_bias = 1-abs(0.5-(as.numeric(by(assays(se)$strand_bias, 
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

getGeneScore = function(se, thresholdIndex, plot = NULL, method = "sgboost"){
  geneFeatures = prepareGeneModelFeatures(se[thresholdIndex,])
    if(checkFeatures(geneFeatures)){
    geneModel = trainGeneModel(geneFeatures$features, 
      geneFeatures$labels, geneFeatures$names, 
      plot = plot, saveFig = F, method)
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

trainGeneModel = function(features, labels, names, plot = NULL, saveFig = T, 
  method = "sgboost"){
  if(sum(labels)==length(labels) | sum(labels)==0){return(NULL)}
  geneScore = calculateGeneScore(features, labels, names, method = method)
  if(!is.null(plot)){
    if(saveFig){
      svg(paste0(savePath,'/',plot,"_ROC.svg"), width = 7, height = 7)
      plotGeneModel(features, labels, method = method)
      dev.off()
    } else{
      plotGeneModel(features, labels, method = method)
    }
    if(saveFig){
      svg(paste0(savePath,plot, "_Pres_Sens.svg"), width = 7, height = 7)
    }
    measure = measureit(score = geneScore$model$score, 
      class = geneScore$model$testLabels, measure = c("SENS", "PREC"))
    measure$PREC[1]=1
    base::plot(measure$PREC~measure$SENS, type = "l")
    colours = c("red","orange","yellow","green", "turquoise", "blue",
      "purple", "plum", "gray","black")
    thresholds = c(.9, .8, .7, .6, .5, .4, .3, .2, .1, 0)
    for(i in 1:10){
      v = measure$SENS[sum(unique(geneScore$model$score)>thresholds[i])+1]
      abline(v = v, lwd=1, lty = 2, col = colours[i])
    }
    if(saveFig){
      dev.off()
    }
  }
  
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

getTranscriptScore = function(se, thresholdIndex, plot = NULL, 
  method = "xgboost"){
  txFeatures = prepareTranscriptModelFeatures(se,
    withAdapters = withAdapters)
  if(checkFeatures(txFeatures)){
    txIndex = thresholdIndex[thresholdIndex %in% 
      which(!rowData(se)$novel)]
    #txIndex = thresholdIndex
    transcriptModel = trainTranscriptModel(txFeatures$features[txIndex,], 
      txFeatures$labels[txIndex], plot = plot, saveFig = F, method = method)
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
trainTranscriptModel = function(features, labels, plot = NULL, saveFig = T, 
  method = "xgboost"){
  print(method)
  result=calculateScore(features, labels, method = method)
  if(!is.null(plot)){
    if(saveFig){
      svg(paste0('/',plot,"_ROC.svg"), width = 7, height = 7)
      plotTranscriptModel(features, labels, plot, method = method)
      dev.off()
    } else{
      plotTranscriptModel(features, labels, plot, method = method)
    }

    if(saveFig){
      svg(paste0(savePath,'/',plot, "_Pres_Sens.svg"), width = 7, height = 7)
    }
    measure = measureit(score = result$score, class = result$testLabels,
      measure = c("SENS", "PREC"))
    base::plot(measure$PREC~measure$SENS, type = "l")
    #draw lines to where the thresholds are at the moment
    colours = c("red","orange","yellow","green", "turquoise", "blue","purple", 
      "plum", "gray","black")
    thresholds = c(.9, .8, .7, .6, .5, .4, .3, .2, .1, 0)
    for(i in 1:10){
      abline(v = measure$SENS[sum(unique(result$score)>thresholds[i])+1], 
        lwd=1, lty = 2, col = colours[i])
    }
    if(saveFig){
      dev.off()
    }
  }
  return(result$cvfit)
}

plotGeneModel = function(features, labels, method = "xgboost"){
  colours = c('black','gray','red','blue', 'green', 'purple', 'orange', 
    'yellow', 'brown', 'pink', 'teal', 'darksalmon', 'darkslategray4', 
    'deeppink2', 'goldenrod', 'burlywood1')
  legend = NULL
  i=1
  result = calculateScore(features,labels, method = method)
  x=rocit(score=result$score, class=result$testLabels)
  base::plot(x$FPR, x$TPR, col=colours[i],  main="Gene Score ROC", type ='l')
  text(0.8,0.6,c(paste0("All-AUC:",as.character(signif(x$AUC,3)))), cex=1.5)
  i = i+1
  legend = c(legend,"all")
  
  #baseline
  x=rocit(score=result$testData[,"numReadsLog"], class=result$testLabels)
  lines(x$FPR, x$TPR, col=colours[i],  main="reads")
  i = i+1
  legend = c(legend,"numReads")
  #random baseline
  lines(c(0,1),c(0,1))

  for(feature in colnames(features)[-which(colnames(features)=="numReadsLog")]){
    plotLine(features[,c("numReadsLog",feature)],labels, colours[i])
    legend = c(legend,feature)
    i = i+1
  }
  legend("bottomright", legend=legend, col = colours[1:i],lty=1, cex=1)
}

plotTranscriptModel = function(features, labels, plot, method = "xgboost"){
  legend = NULL
  i= 1
  colours = c('black','gray','red','blue', 'green', 'purple', 'orange', 
    'yellow', 'brown', 'pink', 'darksalmon', 'darkslategray4', 'deeppink2', 
    'goldenrod', 'burlywood1')
  
  result = calculateScore(features, labels, method = method)
  x=rocit(score=result$score, class=result$testLabels)
  base::plot(x$FPR, x$TPR, col=colours[i],  main="Transcript Score ROC", 
    type ='l')
  text(0.4,0.2,c(paste0("Trained on:",as.character(result$trainCount))), 
    cex=1.5)
  text(0.4,0.1,c(paste0("Tested on:",as.character(length(result$testLabels)))),
    cex=1.5)
  text(0.8,0.6,c(paste0("All-AUC:",as.character(signif(x$AUC,3)))), cex=1.5)
  legend = c(legend, c("All"))
  i = i+1
  
  x=rocit(score=result$testData[,1], class=result$testLabels)
  lines(x$FPR, x$TPR, col=colours[i],  main="reads")
  text(0.8,0.5,c(paste0("Read-AUC:",as.character(signif(x$AUC,3)))), cex=1.5)
  legend = c(legend, c("Number of Reads"))
  i = i+1
  
  lines(c(0,1),c(0,1))
  
  # for(feature in colnames(features)[-which(colnames(features)=="numReads")]){
  #   plotLine(features[,c("numReads",feature)],labels, colours[i])
  #   legend = c(legend,features)
  #   i = i+1
  # }

  for(feature in colnames(features)[-which(colnames(features)=="numReads")]){
    plotLine(features[,feature],labels, colours[i])
    legend = c(legend,feature)
    i = i+1
  }
  legend("bottomright", legend=legend, col = colours[1:i],lty=1, cex=1)
}

plotLine = function(features, labels, colour, log = F, method = "xgboost"){
  if(!is.matrix(features)){
    score = features
  } else {
    result = calculateScore(features,labels, log = log, method = method)
    score = result$score
    labels = result$testLabels
  }
  x=rocit(score=score, class=labels)
  lines(x$FPR, x$TPR, col=colour,  main="reads")
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
  label = train_labels, nthread=2, max.depth = 3, nround= 100, 
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
