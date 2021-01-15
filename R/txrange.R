suppressMessages(require(stringr))
library(BSgenome)
library(ROCit)
library(glmnet)

txrange.filterReadClasses = function(se, readGrgList, genomeSequence, annotations, withAdapters = FALSE, min.readCount = 5){
  
    save(se, file="se_bambuTest.Rdata")
    save(readGrgList, file="readGrgList_bambuTest.Rdata")
    
    if(F){
      load("se_bambuTest.Rdata")
      load("readGrgList_bambuTest.Rdata")
      annotations = bambuAnnotations
      genomeSequence = fa.file
      withAdapters = FALSE
      
      load("se_bambuTest_13_1_21.Rdata")
      load("alignData_bambuTest_13_1_21.Rdata")
    }
    
    if(F){
      bamFile <- Rsamtools::BamFile(test.bam)
      bf <- open(bamFile)
      readGrgList <- GenomicAlignments::grglist(GenomicAlignments::readGAlignments(bf,
        param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isSecondaryAlignment = FALSE)),
        use.names = FALSE))
      readGrgList <- readGrgList[GenomicRanges::width(readGrgList) > 1]
      mcols(readGrgList)$id <- seq_along(readGrgList)
      load("C:/Users/simandred/Desktop/txrange-dev/results/GIS_Hct116_directcDNA_Rep3-Run1/singleOutput_GIS_Hct116_directcDNA_Rep3-Run1.bamdev.Rdata")  
      se = singleOutput[[1]]$se
      colnames(rowData(se)) = c("chr.rc", "strand.rc", "intronStarts", "intronEnds", "confidenceType", "start.rc", "end.rc")
      rowData(combinedOutputs)$start.rc = as.numeric(rowData(combinedOutputs)$start.rc)
      rowData(combinedOutputs)$end.rc = as.numeric(rowData(combinedOutputs)$end.rc)
      rm(singleOutput)
      load("C:/Users/simandred/Desktop/txrange-dev/results/GIS_Hct116_directcDNA_Rep3-Run1/GIS_Hct116_directcDNA_Rep3-Run1.bam_Barcoded_assignedRCdev.Rdata")
    }
    
    print("txRange start")
    options(scipen = 999)
    alignData = createAlignData(readGrgList)
    #rm(readGrgList)
    
    alignData = annotateReadStartsAndEnds(alignData, se)
    print("annotated read start and ends")
    resultOutput = summeriseReadsByGroup(alignData, withAdapters = FALSE)
    resultOutput = resultOutput[order(as.numeric(resultOutput$readClassIndex)),] #tidy up outputs
    #counts = assays(se)$counts
    assays(se, withDimnames = F) = lapply(resultOutput, as.matrix)
    #assays(se)$counts = counts
    rownames(se)=resultOutput$readClass
    combinedOutputs = combineSEs(list(se), annotations)
    combinedOutputs = addRowData(combinedOutputs, genomeSequence, annotations)

    rowData(combinedOutputs)$novel = grepl("gene.", rowData(combinedOutputs)$GENEID)
    seFiltered = combinedOutputs[which(!rowData(combinedOutputs)$novel),]
    #seFiltered = combinedOutputs[which(rowData(combinedOutputs)$non_compatible_match != "NOVEL" | rowData(combinedOutputs)$compatible != "NOVEL"),]
    #seFiltered = combinedOutputs[which(rowData(combinedOutputs)$non_compatible_match != "NOVEL"),]
    seTrimmed = seFiltered[which(rowSums(assays(seFiltered)$counts)>=min.readCount)]
    seWithNovel = combinedOutputs
    seWithNovel = seWithNovel[which(rowSums(assays(seWithNovel)$counts)>=min.readCount)]

    # seTrimmed = seTrimmed[which(rowData(seTrimmed)$chr!="chrIS"),]
    # seWithNovel = seWithNovel[which(rowData(seWithNovel)$chr!="chrIS"),]
    
    geneFeatures = prepareGeneModelFeatures(seWithNovel)
    geneModel = trainGeneModel(geneFeatures$features, geneFeatures$labels, geneFeatures$names, plot = "", saveFig = F)
    
    start_time = Sys.time()
    transcriptFeatures = prepareTranscriptModelFeatures(seTrimmed)
    
    #temp line to deal with NA features (there should be none when working fine)
    transcriptFeatures$features[which(is.na(transcriptFeatures$features))]=0
    
    transcriptModel = trainTranscriptModel(transcriptFeatures$features, transcriptFeatures$labels, plot = "", saveFig = F)
    print(start_time - Sys.time())
    
    # subset = (sapply(str_split(rowData(seTrimmed)$compatible, ';'),length)>=2 | (sapply(str_split(rowData(seTrimmed)$compatible, ';'),length)==1 & rowData(seTrimmed)$equal == "NOVEL" & rowData(seTrimmed)$compatible != "non_compatible_match"))
    # nonSubsetTranscriptModel = trainTranscriptModel(transcriptFeatures$features[!subset,], transcriptFeatures$labels[!subset], plot = "subscriptScore/nonSubsetModel/")
    # #subsetTranscriptModel = trainTranscriptModel(transcriptFeatures$features[subset,], transcriptFeatures$labels[subset], plot = "subscriptScore/subsetModel/", saveFig = F)
    # subsetTranscriptModel = trainTranscriptModel(transcriptFeatures$features[subset,], transcriptFeatures$labels[subset], plot = "subscriptScore/subsetModel/")

    #filter the se by removing unlikely novel genes and read classes
    #geneScore
    geneScore = getGeneScore(geneFeatures$features, geneFeatures$labels, geneFeatures$names, model = geneModel)$score[rowData(seWithNovel)$GENEID]
    seSorted = seWithNovel[order(geneScore, decreasing = T),]  
    rowData(seSorted)$FDR = cumsum(!grepl("gene.",rowData(seSorted)$GENEID))/(1:nrow(rowData(seSorted)))
    threshold = .98
    thresholdIndex = which.max(cumsum(rowData(seSorted)$FDR > threshold))
    seFrac = seSorted[1:thresholdIndex, ]

    #transcriptScore
    features = prepareTranscriptModelFeatures(seFrac, withAdapters = withAdapters)$features
    transcriptScore = predict(transcriptModel, newx = features, s = "lambda.min", type="response")
    seSorted = seFrac[order(transcriptScore, decreasing = T),]
    rowData(seSorted)$FDR = cumsum(rowData(seSorted)$equal != "NOVEL")/(1:nrow(rowData(seSorted)))
    threshold = .9
    thresholdIndex = which.max(cumsum(rowData(seSorted)$FDR > threshold))
    seFrac = seSorted[1:thresholdIndex, ]

    print(seFrac)
    
    return(seFrac)
}

createAlignData = function(readGrgList){
  seqnames = unlist(runValue(seqnames(readGrgList)))
  strand = unlist(runValue(strand(readGrgList)))
  start = min(start(readGrgList))
  end = max(end(readGrgList))
  qwidth = sum(width(readGrgList))
  readClass = mcols(readGrgList)$readClass
  readClassIndex = mcols(readGrgList)$readClassIndex
  alignData = as.data.frame(cbind(seqnames,strand,qwidth, start,end, readClass, readClassIndex))
  return(alignData)
}

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

simplifyAdapterData = function(alignData){
  alignData$adap5 = grepl('5', alignData$predictedStrand, fixed=TRUE)
  alignData$adap3 = grepl('3', alignData$predictedStrand, fixed=TRUE)

  alignData$adap5pos = grepl('5+', alignData$predictedStrand, fixed=TRUE)
  alignData$adap3pos = grepl('3+', alignData$predictedStrand, fixed=TRUE)
  alignData$adap5neg = grepl('5-', alignData$predictedStrand, fixed=TRUE)
  alignData$adap3neg = grepl('3-', alignData$predictedStrand, fixed=TRUE)

  #alignData$adapPos = (alignData$adap5pos & alignData$adap3pos & alignData$txstrand=="+") | (alignData$adap5neg & alignData$adap3neg & alignData$txstrand=="+")
  #alignData$adapNeg = (alignData$adap5neg & alignData$adap3neg & alignData$txstrand=="-") | (alignData$adap5pos & alignData$adap3pos & alignData$txstrand=="-")
  
  alignData$adapPos = (alignData$adap5pos & alignData$adap3pos & alignData$strand=="+") | (alignData$adap5neg & alignData$adap3neg & alignData$strand=="-")
  alignData$adapNeg = (alignData$adap5neg & alignData$adap3neg & alignData$strand=="+") | (alignData$adap5pos & alignData$adap3pos & alignData$strand=="-")
  
  alignData$adapDist5 = as.numeric(alignData$barcode5Prime.minus.start)
  pos5 = which(grepl('5+',alignData$predictedStrand, fixed=TRUE))
  alignData$adapDist5[pos5] = as.numeric(alignData$widthRead5prime[pos5])-as.numeric(alignData$barcode5Prime.plus.end[pos5])
  alignData$adapDist5[!alignData$adap5] = NA
  
  alignData$adapDist3 = as.numeric(alignData$barcode3Prime.minus.start)
  pos3 = which(grepl('3-',alignData$predictedStrand, fixed=TRUE))
  alignData$adapDist3[pos3] = as.numeric(alignData$widthRead3prime[pos3])-as.numeric(alignData$barcode3Prime.plus.end[pos3])
  alignData$adapDist3[!alignData$adap3] = NA

  print("txstrand +, both + adapters")  
  print(sum(alignData$txstrand == "+" & alignData$adapPos))
  print("txstrand +, both - adapters")
  print(sum(alignData$txstrand == "+" & alignData$adapNeg))
  print("txstrand -, both + adapters")
  print(sum(alignData$txstrand == "-" & alignData$adapPos))
  print("txstrand -, both - adapters")
  print(sum(alignData$txstrand == "-" & alignData$adapNeg))
  
  #polyA data
  #alignData$polyApos = rep(0,nrow(alignData))
  alignData$polyApos = as.numeric(alignData$polyA3Prime.plus.score) * 100/as.numeric(alignData$polyA3Prime.plus.pid)
  alignData$polyApos[which(alignData$strand=="-")] = as.numeric(alignData$polyA3Prime.minus.score)[which(alignData$strand=="-")] * 100/as.numeric(alignData$polyA3Prime.minus.pid)[which(alignData$strand=="-")]
  alignData$polyAneg = as.numeric(alignData$polyA3Prime.minus.score) * 100/as.numeric(alignData$polyA3Prime.minus.pid)
  alignData$polyAneg[which(alignData$strand=="-")] = as.numeric(alignData$polyA3Prime.plus.score)[which(alignData$strand=="-")] * 100/as.numeric(alignData$polyA3Prime.plus.pid)[which(alignData$strand=="-")]
  
  #if polyA is present opposide side of sequencing start site
  alignData$polyAFull = alignData$txstrand=="+" & alignData$polyApos > 10 
  print("txstrand +, polyA+")  
  print(sum(alignData$txstrand=="+" & alignData$polyApos > 10,na.rm=T)) #377282
  print("txstrand +, polyA-")  
  print(sum(alignData$txstrand=="+" & alignData$polyAneg > 10,na.rm=T)) #76483
  print("txstrand -, polyA+")  
  print(sum(alignData$txstrand=="-" & alignData$polyApos > 10,na.rm=T)) #129012
  print("txstrand -, polyA-")  
  print(sum(alignData$txstrand=="-" & alignData$polyAneg > 10,na.rm=T)) #638649
  
  #alignData = alignData[,!names(alignData) %in% c("predictedStrand","barcode5Prime.minus.start","barcode3Prime.minus.start","barcode5Prime.plus.end","barcode3Prime.plus.end","widthRead5prime","widthRead3prime")]
  return(alignData)
}

summeriseReadsByGroup <- function(alignData, withAdapters = FALSE){
  
  alignData=alignData[which(!is.na(alignData$readClass)),]

  if(withAdapters){
    result <- alignData %>%
      dplyr::select(readClassIndex, readClass, strand, txstrand, uniqueRead, tssStart, tesEnd, adap5, adap3, adapPos, adapNeg, adapDist5, adapDist3, polyApos, polyAneg) %>%
      group_by(readClassIndex) %>%
      summarise(readClass = readClass[1],
                counts = n(),
                #strand_bias=sum(strand == '+')/counts, 
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
                #tssStart[which.max(cumSum(tssStart==lag(tssStart))))
                #peakPosition = names(which.max(table(tssStart))),
                #peakRelStrength = max(table(tssStart))/counts,
                #endPeakPosition = names(which.max(table(tesEnd))),
                #endPeakRelStrength = max(table(tesEnd))/counts,
      )
  } else{
    result <- alignData %>%
      dplyr::select(readClassIndex, readClass, strand, txstrand, tssStart, tesEnd) %>%
      group_by(readClassIndex) %>%
      summarise(readClass = readClass[1],
                counts = n(),
                #strand_bias=sum(strand == '+')/counts, 
                strand_bias = sum(txstrand == '+'),
                startSD = sd(tssStart),
                endSD = sd(tesEnd),
                #uniqueReads = sum(uniqueRead)
      )
  }
  rownames(result)=result$readClass
  return(result)
}

combineSEs = function(combinedOutputs, annotations){
  combinedTxCandidates = NULL
  for(i in 1:length(combinedOutputs)){
    combinedTxCandidates <- isore.combineTranscriptCandidates(combinedOutputs[[i]], readClassSeRef = combinedTxCandidates)
  }
  return(se)
}


addRowData = function(combinedOutputs, genomeSequence, annotations){
  
  start_time = Sys.time()

  exons = str_split(rowData(combinedOutputs)$intronStarts,",")
  rowData(combinedOutputs)$numExons = sapply(exons, FUN = length)+1
  rowData(combinedOutputs)$numExons[is.na(exons)] = 1
  print("got numExons")
  print(Sys.time()-start_time)
  
  readClassesList = convertSEtoGRangesList(combinedOutputs)
  print("convertSEtoGRangesList")
  
  #######
  #alternative equal/compatible finder
  classifications = getReadClassClassifications(readClassesList, annotations)
  rowData(combinedOutputs)$equal = classifications$equal
  rowData(combinedOutputs)$compatible = classifications$compatible
  print("testing new equal/compat classifier")
  print(Sys.time()-start_time)
  
  #################################################################################
  #calculateDistToAnnotation to match classes to annotation
  # print(Sys.time()-start_time)
  # distNewTx <- calculateDistToAnnotation(exByTx=readClassesList,
  #                                        exByTxRef=annotations)
  # 
  # print("assigned annotations to read classes")
  # print(Sys.time()-start_time)
  # ###############################
  # # classify read classes as equal/compatible/overlap/novel
  # 
  # #readClasses = sort(unique(as.numeric(resultOutput$readClassIndex)))
  # rowData(combinedOutputs)$equal = getEqualReadClasses(distNewTx, nrow(combinedOutputs))
  # compatibleOutput = getCompatibleReadClasses(distNewTx, nrow(combinedOutputs))
  # rowData(combinedOutputs)$compatible = compatibleOutput$compatible
  # rowData(combinedOutputs)$non_compatible_match = compatibleOutput$non_compatible_match
  # print("got compatibles")
  # print(Sys.time()-start_time)
  ####
  
  combinedOutputs = assignGeneIDs(combinedOutputs, annotations)
  print("assigned geneID")
  print(Sys.time()-start_time)
  combinedOutputs = calculateGeneProportion(combinedOutputs)
  print("gotGeneProp")
  print(Sys.time()-start_time)
  
  #assign the max transcript proportions for each read class
  combinedOutputs = getTranscriptProp(combinedOutputs, readClassesList)
  print("gotTxProp")
  print(Sys.time()-start_time)
  #hist(alignedPolyAstart[which(resultOutput$equal=="NOVEL" & (resultOutput$non_compatible_match!="NOVEL" | resultOutput$compatible!="NOVEL"))])
  
  combinedOutputs = countPolyATerminals(combinedOutputs, genomeSequence)
  print("gotPolyAStats")
  print(Sys.time()-start_time)
  
  return(combinedOutputs)
}

assignGeneIDs <- function(se, annotationGrangesList, min.exonOverlap = 35){
  
  start_time = Sys.time()
  
  ## (1) Spliced Reads
  #se <- isore.combineTranscriptCandidates(se) #does this do anything for just one sample other than making start and end columns?
  seFilteredSpliced <- se[(rowData(se)$confidenceType == 'highConfidenceJunctionReads' | rowData(se)$confidenceType == 'lowConfidenceJunctionReads'),]
  
  mcols(seFilteredSpliced)$GENEID <- NA
  
  intronsByReadClass= makeGRangesListFromFeatureFragments(seqnames=rowData(seFilteredSpliced)$chr.rc,
                                                          fragmentStarts=rowData(seFilteredSpliced)$intronStarts,
                                                          fragmentEnds=rowData(seFilteredSpliced)$intronEnds,
                                                          strand=rowData(seFilteredSpliced)$strand.rc)
  
  names(intronsByReadClass) <- 1:length(intronsByReadClass)
  seqlevels(intronsByReadClass) <-  unique(c(seqlevels(intronsByReadClass), seqlevels(annotationGrangesList)))
  
  exonEndsShifted <-paste(rowData(seFilteredSpliced)$intronStarts,
                          rowData(seFilteredSpliced)$end.rc + 1,
                          sep=',')
  exonStartsShifted <- paste(rowData(seFilteredSpliced)$start.rc - 1,
                             rowData(seFilteredSpliced)$intronEnds,
                             sep=',')
  
  exonsByReadClass <- makeGRangesListFromFeatureFragments(seqnames=rowData(seFilteredSpliced)$chr.rc,
                                                          fragmentStarts=exonStartsShifted,
                                                          fragmentEnds=exonEndsShifted,
                                                          strand=rowData(seFilteredSpliced)$strand.rc)
  exonsByReadClass <- narrow(exonsByReadClass, start = 2, end = -2)  # correct junction to exon differences in coordinates
  names(exonsByReadClass) <- 1:length(exonsByReadClass)
  
  # add exon start and exon end rank
  unlistData <- unlist(exonsByReadClass, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(elementNROWS(exonsByReadClass)), names=NULL)
  
  exon_rank <- sapply(width((partitioning)), seq, from=1)
  exon_rank[which(rowData(seFilteredSpliced)$strand == '-')] <- lapply(exon_rank[which(rowData(seFilteredSpliced)$strand == '-')], rev)  # * assumes positive for exon ranking
  exon_endRank <- lapply(exon_rank, rev)
  unlistData$exon_rank <- unlist(exon_rank)
  unlistData$exon_endRank <- unlist(exon_endRank)
  
  exonsByReadClass <- relist(unlistData, partitioning)
  
  seqlevels(exonsByReadClass) <-  unique(c(seqlevels(exonsByReadClass), seqlevels(annotationGrangesList)))
  
  ovExon <- findSpliceOverlapsQuick(cutStartEndFromGrangesList(exonsByReadClass),
                                    cutStartEndFromGrangesList(annotationGrangesList))
  
  ## annotate with transcript and gene Ids
  mcols(seFilteredSpliced)$GENEID[queryHits(ovExon[mcols(ovExon)$compatible][!duplicated(queryHits(ovExon[mcols(ovExon)$compatible]))])] <- mcols(annotationGrangesList[subjectHits(ovExon[mcols(ovExon)$compatible])[!duplicated(queryHits(ovExon[mcols(ovExon)$compatible]))]])$GENEID # annotate with compatible gene id,
  mcols(seFilteredSpliced)$GENEID[queryHits(ovExon[mcols(ovExon)$equal][!duplicated(queryHits(ovExon[mcols(ovExon)$equal]))])] <- mcols(annotationGrangesList[subjectHits(ovExon[mcols(ovExon)$equal])[!duplicated(queryHits(ovExon[mcols(ovExon)$equal]))]])$GENEID # annotate as identical,

  ## using intron matches
  unlistedIntrons <- unlist(intronsByReadClass, use.names=TRUE)
  partitioning <- PartitioningByEnd(cumsum(elementNROWS(intronsByReadClass)), names=NULL)
  
  unlistedIntronsAnnotations <- unlist(myGaps(annotationGrangesList))
  mcols(unlistedIntronsAnnotations)$GENEID <- mcols(annotationGrangesList)$GENEID[match(names(unlistedIntronsAnnotations), mcols(annotationGrangesList)$TXNAME)]
  
  intronMatches <- GenomicRanges::match(unlistedIntrons, unique(unlistedIntronsAnnotations), nomatch=0) > 0
  intronMatchesList <- relist(intronMatches, partitioning)
  
  lastJunctionMatch <- unlist(endoapply(endoapply(intronMatchesList, rev), '[[', 1))
  firstJunctionMatch <- unlist(endoapply(intronMatchesList, '[[', 1))
  
  ## assign gene ids based on the maximum number of matching introns/splice junctions
  overlapsNewIntronsAnnotatedIntrons <- findOverlaps(unlistedIntrons,
                                                     unlistedIntronsAnnotations,
                                                     type='equal',
                                                     select='all',
                                                     ignore.strand=FALSE)
  
  maxGeneCountPerNewTx <- tbl_df(data.frame(txId=names(unlistedIntrons)[queryHits(overlapsNewIntronsAnnotatedIntrons)],
                                            geneId=mcols(unlistedIntronsAnnotations)$GENEID[subjectHits(overlapsNewIntronsAnnotatedIntrons)],
                                            stringsAsFactors=FALSE)) %>%
    group_by(txId, geneId) %>%
    summarise(geneCount=n()) %>%
    group_by(txId) %>%
    filter(geneCount == max(geneCount)) %>%
    filter(!duplicated(txId)) %>%
    ungroup()
  
  geneIdByIntron <- rep(NA,length(exonsByReadClass))
  geneIdByIntron <- maxGeneCountPerNewTx$geneId[match(names(exonsByReadClass), maxGeneCountPerNewTx$txId)]
  mcols(seFilteredSpliced)$GENEID[is.na(mcols(seFilteredSpliced)$GENEID)] <- geneIdByIntron[is.na(mcols(seFilteredSpliced)$GENEID)]
  
  ## unspliced transcripts
  if(any(rowData(se)$confidenceType == 'unsplicedNew')) {
    seFilteredUnspliced <- se[(rowData(se)$confidenceType == 'unsplicedNew' | rowData(se)$confidenceType == 'unsplicedWithin'), ]
    exonsByReadClassUnspliced= GRanges(seqnames=rowData(seFilteredUnspliced)$chr.rc,
                                       ranges=IRanges(start=rowData(seFilteredUnspliced)$start.rc,
                                                      end=rowData(seFilteredUnspliced)$end.rc),
                                       strand=rowData(seFilteredUnspliced)$strand.rc)
    
    partitioning <- PartitioningByEnd(1:length(exonsByReadClassUnspliced), names=NULL)
    exonsByReadClassUnspliced$exon_rank <- rep(1, length(exonsByReadClassUnspliced))
    exonsByReadClassUnspliced$exon_endRank <- rep(1, length(exonsByReadClassUnspliced))
    exonsByReadClassUnspliced <- relist(exonsByReadClassUnspliced, partitioning)
    seqlevels(exonsByReadClassUnspliced) <-  unique(c(seqlevels(exonsByReadClassUnspliced), seqlevels(annotationGrangesList)))
    mcols(seFilteredUnspliced)$GENEID <- NA
    #mcols(seFilteredUnspliced)$readClassType <- 'unsplicedNew'
    
    ## here: add filter to remove unspliced transcripts which overlap with known transcripts/high quality spliced transcripts
    # overlapUnspliced <- findOverlaps(exonsByReadClassUnspliced,
    #                                  annotationGrangesList,
    #                                  minoverlap=min.exonOverlap,
    #                                  select='first')
    # seFilteredUnspliced <- seFilteredUnspliced[is.na(overlapUnspliced)]
    # exonsByReadClassUnspliced <- exonsByReadClassUnspliced[is.na(overlapUnspliced)]
    seFilteredUnspliced <- seFilteredUnspliced
    exonsByReadClassUnspliced <- exonsByReadClassUnspliced
    
    ## combined spliced and unspliced Tx candidates
    seCombined <- SummarizedExperiment::rbind(seFilteredSpliced, seFilteredUnspliced)
    exonRangesCombined<- c(exonsByReadClass, exonsByReadClassUnspliced)
    names(exonRangesCombined) <- 1:length(exonRangesCombined)
  } else {
    seCombined <- seFilteredSpliced
    exonRangesCombined<- exonsByReadClass
    names(exonRangesCombined) <- 1:length(exonRangesCombined)
  }
  exonMatchGene <- findOverlaps(exonRangesCombined,
                                annotationGrangesList,
                                select='arbitrary',
                                minoverlap=min.exonOverlap)
  geneIdByExon <- rep(NA,length(exonRangesCombined))
  geneIdByExon[!is.na(exonMatchGene)] <- mcols(annotationGrangesList)$GENEID[exonMatchGene[!is.na(exonMatchGene)]]
  geneIdByExon[!is.na(mcols(seCombined)$GENEID)] <-  mcols(seCombined)$GENEID[!is.na(mcols(seCombined)$GENEID)]
  exonMatchGene <- findOverlaps(exonRangesCombined[is.na(geneIdByExon)],
                                exonRangesCombined[!is.na(geneIdByExon)],
                                select = 'arbitrary',
                                minoverlap = min.exonOverlap)
  while(any(!is.na(exonMatchGene))) {
    geneIdByExon[is.na(geneIdByExon)][!is.na(exonMatchGene)] <- geneIdByExon[!is.na(geneIdByExon)][exonMatchGene[!is.na(exonMatchGene)]]
    exonMatchGene <- findOverlaps(exonRangesCombined[is.na(geneIdByExon)],
                                  exonRangesCombined[!is.na(geneIdByExon)],
                                  select = "arbitrary", minoverlap = min.exonOverlap)
  }
  mcols(seCombined)$GENEID[is.na(mcols(seCombined)$GENEID)] <- geneIdByExon[is.na(mcols(seCombined)$GENEID)]
  if(any(is.na(mcols(seCombined)$GENEID))){
    
    newGeneIds <- assignNewGeneIds(exonRangesCombined[is.na(mcols(seCombined)$GENEID)], minoverlap=5, ignore.strand=F)
    
    mcols(seCombined)$GENEID[as.integer(newGeneIds$readClassId)] <- newGeneIds$geneId
  }
  
  return(seCombined)
}

assignGeneIDs2 = function(se, annotations){
  mcols(se)$GENEID <- NA
  mcols(se)$GENEID[queryHits(ovExon[mcols(ovExon)$compatible][!duplicated(queryHits(ovExon[mcols(ovExon)$compatible]))])] <- mcols(annotationGrangesList[subjectHits(ovExon[mcols(ovExon)$compatible])[!duplicated(queryHits(ovExon[mcols(ovExon)$compatible]))]])$GENEID # annotate with compatible gene id,
  mcols(se)$GENEID[queryHits(ovExon[mcols(ovExon)$equal][!duplicated(queryHits(ovExon[mcols(ovExon)$equal]))])] <- mcols(annotationGrangesList[subjectHits(ovExon[mcols(ovExon)$equal])[!duplicated(queryHits(ovExon[mcols(ovExon)$equal]))]])$GENEID # annotate as identical,
  exonsByReadClass = getExonsByReadClass(se)
  se[(rowData(se)$confidenceType == 'highConfidenceJunctionReads' | rowData(se)$confidenceType == 'lowConfidenceJunctionReads'),]$GENEID = geneAssignmentByIntrons(se, annotations, exonsByReadClass$spliced)
  se = geneAssignmentArbitraryOverlap(se, annotations, c(exonsByReadClass$spliced, exonsByReadClass$unspliced))
  
  getExonsByReadClass = function(se){
    seFilteredSpliced <- se[(rowData(se)$confidenceType == 'highConfidenceJunctionReads' | rowData(se)$confidenceType == 'lowConfidenceJunctionReads'),]
    exonEndsShifted <-paste(rowData(seFilteredSpliced)$intronStarts,
                            rowData(seFilteredSpliced)$end.rc + 1,
                            sep=',')
    exonStartsShifted <- paste(rowData(seFilteredSpliced)$start.rc - 1,
                               rowData(seFilteredSpliced)$intronEnds,
                               sep=',')
    exonsByReadClass <- makeGRangesListFromFeatureFragments(seqnames=rowData(seFilteredSpliced)$chr.rc,
                                                            fragmentStarts=exonStartsShifted,
                                                            fragmentEnds=exonEndsShifted,
                                                            strand=rowData(seFilteredSpliced)$strand.rc)
    exonsByReadClass <- narrow(exonsByReadClass, start = 2, end = -2)  # correct junction to exon differences in coordinates
    names(exonsByReadClass) <- 1:length(exonsByReadClass)
    
    seFilteredUnspliced <- se[(rowData(se)$confidenceType == 'unsplicedNew' | rowData(se)$confidenceType == 'unsplicedWithin'), ]
    exonsByReadClassUnspliced= GRanges(seqnames=rowData(seFilteredUnspliced)$chr.rc,
                                       ranges=IRanges(start=rowData(seFilteredUnspliced)$start.rc,
                                                      end=rowData(seFilteredUnspliced)$end.rc),
                                       strand=rowData(seFilteredUnspliced)$strand.rc)
    #exonRangesCombined<- c(exonsByReadClass, exonsByReadClassUnspliced)
    #names(exonRangesCombined) <- 1:length(exonRangesCombined)
    return(list(spliced = exonsByReadClass, unspliced = exonsByReadClassUnspliced))
  }
  
  geneAssignmentByIntrons = function(se, annotationGrangesList, exonsByReadClass){
    ## using intron matches
    se <- se[(rowData(se)$confidenceType == 'highConfidenceJunctionReads' | rowData(se)$confidenceType == 'lowConfidenceJunctionReads'),]
    
    intronsByReadClass= makeGRangesListFromFeatureFragments(seqnames=rowData(se)$chr.rc,
                                                            fragmentStarts=rowData(se)$intronStarts,
                                                            fragmentEnds=rowData(se)$intronEnds,
                                                            strand=rowData(se)$strand.rc)
    unlistedIntrons <- unlist(intronsByReadClass, use.names=TRUE)
    
    unlistedIntronsAnnotations <- unlist(myGaps(annotationGrangesList))
    mcols(unlistedIntronsAnnotations)$GENEID <- mcols(annotationGrangesList)$GENEID[match(names(unlistedIntronsAnnotations), mcols(annotationGrangesList)$TXNAME)]
    
    ## assign gene ids based on the maximum number of matching introns/splice junctions
    overlapsNewIntronsAnnotatedIntrons <- findOverlaps(unlistedIntrons,
                                                       unlistedIntronsAnnotations,
                                                       type='equal',
                                                       select='all',
                                                       ignore.strand=FALSE)
    
    maxGeneCountPerNewTx <- tbl_df(data.frame(txId=names(unlistedIntrons)[queryHits(overlapsNewIntronsAnnotatedIntrons)],
                                              geneId=mcols(unlistedIntronsAnnotations)$GENEID[subjectHits(overlapsNewIntronsAnnotatedIntrons)],
                                              stringsAsFactors=FALSE)) %>%
      group_by(txId, geneId) %>%
      summarise(geneCount=n()) %>%
      group_by(txId) %>%
      filter(geneCount == max(geneCount)) %>%
      filter(!duplicated(txId)) %>%
      ungroup()


    geneIdByIntron <- rep(NA,length(exonsByReadClass))
    geneIdByIntron <- maxGeneCountPerNewTx$geneId[match(names(exonsByReadClass), maxGeneCountPerNewTx$txId)]
    mcols(se)$GENEID[is.na(mcols(se)$GENEID)] <- geneIdByIntron[is.na(mcols(se)$GENEID)]
    return(mcols(se)$GENEID)
  }
  
  geneAssignmentArbitraryOverlap = function(exonRangesCombined, annotationGrangesList){
    exonMatchGene <- findOverlaps(exonRangesCombined,
                                  annotationGrangesList,
                                  select='arbitrary',
                                  minoverlap=min.exonOverlap)
    geneIdByExon <- rep(NA,length(exonRangesCombined))
    geneIdByExon[!is.na(exonMatchGene)] <- mcols(annotationGrangesList)$GENEID[exonMatchGene[!is.na(exonMatchGene)]]
    geneIdByExon[!is.na(mcols(se)$GENEID)] <-  mcols(se)$GENEID[!is.na(mcols(se)$GENEID)]
    
    exonMatchGene <- findOverlaps(exonRangesCombined[is.na(geneIdByExon)],
                                  exonRangesCombined[!is.na(geneIdByExon)],
                                  select = 'arbitrary',
                                  minoverlap = min.exonOverlap)
    while(any(!is.na(exonMatchGene))) {
      geneIdByExon[is.na(geneIdByExon)][!is.na(exonMatchGene)] <- geneIdByExon[!is.na(geneIdByExon)][exonMatchGene[!is.na(exonMatchGene)]]
    }
    mcols(se)$GENEID[is.na(mcols(se)$GENEID)] <- geneIdByExon[is.na(mcols(se)$GENEID)]
    
    if(any(is.na(mcols(se)$GENEID))){
      
      newGeneIds <- assignNewGeneIds(exonRangesCombined[is.na(mcols(se)$GENEID)], minoverlap=5, ignore.strand=F)
      
      mcols(se)$GENEID[as.integer(newGeneIds$readClassId)] <- newGeneIds$geneId
    }
    return(se)
  }
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
  rowData(resultOutput)$totalGeneReadProp = rowSums(assays(resultOutput)$counts / rowSums(countsTBL))
  return(resultOutput)
}

convertSEtoGRangesList = function(combinedOutputs){
  
  exonEndsShifted <- paste(rowData(combinedOutputs)$intronStarts, rowData(combinedOutputs)$end.rc + 1, sep = ',')
  exonStartsShifted <- paste(rowData(combinedOutputs)$start.rc - 1, rowData(combinedOutputs)$intronEnds, sep = ',')
  #deal with single exon RCs
  singleExonIndex = which(is.na(rowData(combinedOutputs)$intronStarts))
  exonEndsShifted[singleExonIndex] = rowData(combinedOutputs)$end.rc[singleExonIndex] + 1
  exonStartsShifted[singleExonIndex] = rowData(combinedOutputs)$start.rc[singleExonIndex] - 1
  
  readClassList= makeGRangesListFromFeatureFragments(seqnames=rowData(combinedOutputs)$chr.rc,
                                                     fragmentStarts=exonStartsShifted,
                                                     fragmentEnds=exonEndsShifted,
                                                     strand=rowData(combinedOutputs)$strand.rc)
  
  readClassList <- narrow(readClassList, start = 2, end = -2)  # correct junction to exon differences in coordinates
  
  unlistData <- unlist(readClassList, use.names = FALSE)
  partitioning <- PartitioningByEnd(cumsum(elementNROWS(readClassList)), names = NULL)
  
  exon_rank <- sapply(width((partitioning)), seq, from = 1)
  exon_rank[which(rowData(combinedOutputs)$strand.rc == '-')] <- lapply(exon_rank[which(rowData(combinedOutputs)$strand.rc == '-')], rev)  # * assumes positive for exon ranking
  exon_endRank <- lapply(exon_rank, rev)
  unlistData$exon_rank <- unlist(exon_rank)
  unlistData$exon_endRank <- unlist(exon_endRank)
  
  readClassList <- relist(unlistData, partitioning)
  names(readClassList) = 1:length(readClassList)
  
  return(readClassList)
}

getReadClassClassifications = function(query, subject, maxDist = 5){

    subjectExtend <- extendGrangesListElements(subject, by = maxDist)
    queryForOverlap <- dropGrangesListElementsByWidth(query, minWidth = maxDist, cutStartEnd = T)
    olap = findOverlaps(queryForOverlap, subjectExtend, ignore.strand = F, type = 'within')
    compatible = rep(F, length(query))
    compatible[queryHits(olap)] = T
    
    olapEqual = findOverlaps(queryForOverlap, cutStartEndFromGrangesList(subject), ignore.strand = F, type = 'equal')
    equal = rep(F, length(query))
    equal[queryHits(olapEqual)] = T
    
    # query2 <- query[queryHits(olap)]
    # subject <- subject[subjectHits(olap)]
    # splice <- myGaps(query2)
    # compatible <- rangesDist(query2, subject, splice, maxDist)
    # equal2 <- (!is.na(S4Vectors::match(olap, olapEqual)))
    
    return(list(equal = equal, compatible = compatible))
}

getEqualReadClasses = function(distNewTx, RClength){
  equal = rep("NOVEL", RClength)
  matchingTranscripts=by(cbind(distNewTx$annotationTxId, distNewTx$equal), distNewTx$queryHits, function(x){
    if("TRUE" %in% x[,2]){
      x = x[which(x[,2]==TRUE),1]
      return(paste(as.character(x),collapse=";"))
    } else {
      return("NOVEL")
    }
  })
  equal[as.numeric(names(matchingTranscripts))]=as.character(matchingTranscripts)

  return(equal)
}

getCompatibleReadClasses = function(distNewTx, RClength){
  #get an arbitrary compatible transcript
  compatible = rep("NOVEL", RClength)
  non_compatible_match = rep(F, RClength)
  matchingTranscripts=by(cbind(distNewTx$annotationTxId, distNewTx$compatible), distNewTx$queryHits, function(x){
    if("TRUE" %in% x[,2]){
      x=x[which(x[,2]==TRUE),1]
      return(paste(as.character(x),collapse=";"))
    } else {
      ids = paste(as.character(x[,1]),collapse=";")
      return(paste0("non_compatible_match_",ids))
    }
  })
  compatible[as.numeric(names(matchingTranscripts))]=as.character(matchingTranscripts)
  
  #assign a gene to the non_compaitible_matches 
  
  selection = which(grepl("non_compatible_match",compatible))
  non_compatible_match[selection] = gsub("non_compatible_match_","",compatible[selection])
  compatible[selection] = "non_compatible_match"
  
  return(list(compatible=compatible, non_compatible_match=non_compatible_match))
}

getTranscriptProp = function(se, readClassList){
  #takes a se and a list of semi colon delinated compatible transcripts/read classes and calculates the max transcript proportion it has
  #align compatible read classes with their counts, their read class name and the transcript they subset
  
  #or just an SE which causes it to find overlaps within the read classes themselves

  start_time = Sys.time()
  
  allOverlaps = NULL
  i= 0
  for(chr in seqlevels(readClassList)){
    selection = which(unlist(runValue(seqnames(readClassList))==chr))
    readClassListTemp = readClassList[selection]
    query <- cutStartEndFromGrangesList(readClassListTemp)
    
    overlaps = countOverlaps(query, readClassListTemp, type = 'within')
    megaRCs = readClassList[which(overlaps == 1),]
    overlaps = findOverlaps(query, megaRCs, type = 'within')
    
    #account for the index change due to seperating the readclasses by chromosome
    overlaps = as.matrix(overlaps)+i
    i = i + length(readClassListTemp)
    allOverlaps = rbind(allOverlaps, overlaps)
  }
  
  print(1)
  print(Sys.time()-start_time)
  
  allCompats = allOverlaps[,'subjectHits']
  compatIndex = allOverlaps[,'queryHits']
  compatCounts = rowSums(assays(se)$counts)[allOverlaps[,'queryHits']]
  
  print(2)
  print(Sys.time()-start_time)
  
  #inverse the data so read classes are grouped now by transcript
  annotationSubsets = by(compatIndex,allCompats, FUN = c)
  #calculate the total expression of a transcript (all its subset read class counts) and the proportion from each RC
  transcriptSums = by(cbind(compatCounts,compatIndex),allCompats, FUN = function(x){
    counts = as.numeric(x[,1])
    transcriptProp = counts/sum(counts)
    names(transcriptProp) =  x[,2]
    return(transcriptProp)
  })
  
  
  print(3)
  print(Sys.time()-start_time)
  
  #align the proportions with the read class names
  transcriptProp = unlist(transcriptSums)
  indexes = unlist(sapply(transcriptSums, FUN = names))
  x = as.data.frame(cbind(transcriptProp, indexes))
  
  print(4)
  print(Sys.time()-start_time)
  
  #get top transcript prop for all the transcripts an RC is part of then break ties by just picking any
  #a RC can have multiple transcript props and the reason I pick the top is so that read classes of subset transcripts are fairly represented
  topTranscriptProp_RC = x %>% group_by(indexes) %>% top_n(1, transcriptProp) %>% distinct(indexes, transcriptProp, .keep_all = TRUE) 
  bottomTranscriptProp_RC = x %>% group_by(indexes) %>% top_n(-1, transcriptProp) %>% distinct(indexes, transcriptProp, .keep_all = TRUE) 

  rowData(se)$transcriptProp = rep(0, nrow(se))
  rowData(se)$transcriptPropMin = rep(0, nrow(se))
  rowData(se)[as.numeric(topTranscriptProp_RC$indexes),]$transcriptProp = topTranscriptProp_RC$transcriptProp
  rowData(se)[as.numeric(bottomTranscriptProp_RC$indexes),]$transcriptPropMin = bottomTranscriptProp_RC$transcriptProp
  
  print(5)
  print(Sys.time()-start_time)
  
  return(se)
}


# generates features counting the presence of A's and T's at RC aligned starts and ends (different from polyA tail which is not aligned)
countPolyATerminals = function(se, genomeSequence){
  
  start_time = Sys.time()
  ##################
  #calculate the presence of polyA and polyT seqences at the start and end of the alignment (not polyA tail)
  #1. get the "start" and "end" of each read class (use .9 percentile like bambu)
  #2. get the sequence of the read class
  #3. in the first/last 10 bp how many A's are there
  #4. in the first/last 10 bp how many T's are there
  #5. Does aligning the A's and T's with strand help
  #6. does extending the window into the non-align sequence help
  
  #get all first/last exons
  RCranges = makeGRangesFromDataFrame(rowData(se),start.field = "start.rc", end.field = "end.rc",
                                      seqnames.field = "chr.rc", strand.field = "strand.rc")
  strand(RCranges)[which(as.character(strand(RCranges))=='*')]='+'
  genomeSequence = loadGenomeAnnotation(genomeSequence)
  
  print(1)
  print(Sys.time()-start_time) 
  
  startSeqs = as.character(BSgenome::getSeq(genomeSequence,resize(RCranges,10, fix="start")[which(as.character(seqnames(RCranges)) %in% names(genomeSequence))]))
  endSeqs = as.character(BSgenome::getSeq(genomeSequence,resize(RCranges,10, fix="end")[which(as.character(seqnames(RCranges)) %in% names(genomeSequence))]))
  
  print(2)
  print(Sys.time()-start_time) 
  
  #startSeqs = BSgenome::getSeq(genomeSequence,resize(firstExons,10, fix="start")[which(seqnames(firstExons) %in% seqnames(genomeSequence))], as.character = T)
  #endSeqs = BSgenome::getSeq(genomeSequence,resize(lastExons,10, fix="end")[which(seqnames(lastExons) %in% seqnames(genomeSequence))], as.character = T)
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
  
  print(3)
  print(Sys.time()-start_time) 
  
  #count highest conseq A's in first/last 10 bp
  alignedPolyAstart = rep(NA,nrow(se))  
  alignedPolyAstart[index] = sapply(startSeqs, FUN = function(x){
    a_rle = rle(strsplit(x,"")[[1]])
    max(a_rle$lengths[which(a_rle$values=="A")])
  })
  alignedPolyAend = rep(NA,nrow(se)) 
  alignedPolyAend[index] = sapply(endSeqs, FUN = function(x){
    a_rle = rle(strsplit(x,"")[[1]])
    max(a_rle$lengths[which(a_rle$values=="A")])
  })
  alignedPolyTstart = rep(NA,nrow(se))  
  alignedPolyTstart[index] = sapply(startSeqs, FUN = function(x){
    a_rle = rle(strsplit(x,"")[[1]])
    max(a_rle$lengths[which(a_rle$values=="T")])
  })
  alignedPolyTend = rep(NA,nrow(se)) 
  alignedPolyTend[index] = sapply(endSeqs, FUN = function(x){
    a_rle = rle(strsplit(x,"")[[1]])
    max(a_rle$lengths[which(a_rle$values=="T")])
  })
  
  print(4)
  print(Sys.time()-start_time) 
  
  rowData(se)$numAstart = numAstart
  rowData(se)$numAend = numAend
  rowData(se)$alignedPolyAstart = alignedPolyAstart
  rowData(se)$alignedPolyAend = alignedPolyAend
  rowData(se)$numTstart = numTstart
  rowData(se)$numTend = numTend
  rowData(se)$alignedPolyTstart = alignedPolyTstart
  rowData(se)$alignedPolyTend = alignedPolyTend
  return(se)
}

prepareTranscriptModelFeatures = function(input, withAdapters = F){

  #labels = rowData(input)$equal != "NOVEL"
  labels = rowData(input)$equal
  
  numReads = rowSums(assays(input)$counts)
  logNumReads = log(numReads, 2)
  
  geneReadProp=rowData(input)$totalGeneReadProp
  geneReadProp[is.na(geneReadProp)]=0
  
  tx_strand_bias=(1-abs(0.5-(rowSums(assays(input)$strand_bias)/numReads)))
  
  twoRead = ifelse(numReads>=2,0,1)
  threeRead = ifelse(numReads>=3,0,1)
  fourRead = ifelse(numReads>=4,0,1)
  fiveRead = ifelse(numReads>=5,0,1)
  
  #numExons = rowData(input)$numExons
  
  SD = apply(cbind(assays(input)$startSD, assays(input)$counts),MARGIN = 1, FUN = function(x){
    weighted.mean(x[1:(length(x)/2)],x[((length(x)/2)+1):length(x)], na.rm = T)
  })*-1
  SD[which(is.na(SD))] = 1
  SDend = apply(cbind(assays(input)$endSD, assays(input)$counts),1,FUN = function(x){
    weighted.mean(x[1:(length(x)/2)],x[((length(x)/2)+1):length(x)], na.rm = T)
  })*-1
  SDend[which(is.na(SDend))] = 1
  
  # uniqueReads = rowSums(assays(input)$uniqueReads)/numReads
  
  numAstart = rowData(input)$numAstart
  numAstart[is.infinite(numAstart)]=0
  numAend = rowData(input)$numAend
  numAend[is.infinite(numAend)]=0
  alignedPolyAstart = rowData(input)$alignedPolyAstart
  alignedPolyAstart[is.infinite(alignedPolyAstart)]=0
  alignedPolyAend = rowData(input)$alignedPolyAend
  alignedPolyAend[is.infinite(alignedPolyAend)]=0
  
  numTstart = rowData(input)$numTstart
  numTstart[is.infinite(numTstart)]=0
  numTend = rowData(input)$numTend
  numTend[is.infinite(numTend)]=0
  alignedPolyTstart = rowData(input)$alignedPolyTstart
  alignedPolyTstart[is.infinite(alignedPolyTstart)]=0
  alignedPolyTend = rowData(input)$alignedPolyTend
  alignedPolyTend[is.infinite(alignedPolyTend)]=0
  
  # TODO better solution to mean of TESscore?
  # TSSscore = rowMeans(assays(input)$TSSscore, na.rm = T)
  # TSSscore[is.na(TSSscore)] = 0
  # TESscore = rowMeans(assays(input)$TESscore, na.rm = T)
  # TESscore[is.na(TESscore)] = 0
  # bothTSSandTES = ifelse((TSSscore > 0.5 & TESscore > 0.5), 1, 0)
  
  transcriptProp = as.numeric(rowData(input)$transcriptProp)
  transcriptProp[is.na(transcriptProp)] = 0
  transcriptPropMin = as.numeric(rowData(input)$transcriptPropMin)
  transcriptPropMin[is.na(transcriptPropMin)] = 0
  
  features = cbind(numReads, SD, SDend, geneReadProp, tx_strand_bias, numAstart, numAend, alignedPolyAstart, alignedPolyAend, numTstart, numTend, alignedPolyTstart, alignedPolyTend, transcriptProp, transcriptPropMin)
  #features = cbind(numReads, SD, SDend, geneReadProp, tx_strand_bias, uniqueReads, numAstart, numAend, alignedPolyAstart, alignedPolyAend, numTstart, numTend, alignedPolyTstart, alignedPolyTend, TSSscore, TESscore, transcriptProp, transcriptPropMin, bothTSSandTES)
  #features = cbind(numReads, SD, SDend, geneReadProp, tx_strand_bias)
  
  
  if("adapNumReads" %in% names(assays(input))){
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
    
    adapter_distance = apply(cbind(assays(input)$adaptDist, assays(input)$counts),1,FUN = function(x){
      weighted.mean(x[1:(length(x)/2)],x[((length(x)/2)+1):length(x)])
    })
    #adapter_distance = assays(input)$adaptDist
    adapter_distance[is.na(adapter_distance)]=100
    adapter_distance[adapter_distance<=40]=1
    adapter_distance[adapter_distance>40]=0
    
    adapterend_distance = apply(cbind(assays(input)$adaptEndDist, assays(input)$counts),1,FUN = function(x){
      weighted.mean(x[1:(length(x)/2)],x[((length(x)/2)+1):length(x)])
    })
    #adapterend_distance = assays(input)$adapEndDist
    adapterend_distance[is.na(adapterend_distance)]=100
    adapterend_distance[adapterend_distance<=30]=1
    adapterend_distance[adapterend_distance>30]=0
    
    bothAdapters = rowSums(assays(input)$bothAdapters)
    
    polyAEnd = rowSums(assays(input)$polyAEnd)
    
    bothAdaptersProp = rowSums(assays(input)$bothAdapters)/numReads
    
    polyAEndProp = rowSums(assays(input)$polyAEnd)/numReads
    
    features = cbind(features, adapter_prop, adapterend_prop, adapter_distance, adapterend_distance, bothAdapters, polyAEnd, bothAdaptersProp, polyAEndProp)
    
    #features = cbind(features, adapter_prop, adapterend_prop, adapter_distance, adapterend_distance)
  }
  
  return(list(features = features, labels = labels))
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

  strand_bias = 1-abs(0.5-(as.numeric(by(assays(se)$strand_bias, rowData(se)$GENEID, sum))/numReads))
  strand_bias[is.na(strand_bias)]=0
  #as.numeric(by(assays(se3)$startSD, rowData(se3)$GENEID, mean))
  #as.numeric(by(assays(se3)$endSD, rowData(se3)$GENEID, mean))
  #how many read classes does a gene have
  numRCs = table(rowData(se)$GENEID)
  
  #number of non-subset read classes
  subset = (sapply(str_split(rowData(se)$resultOutput.compatible, ';'),length)>=2 | (sapply(str_split(rowData(se)$resultOutput.compatible, ';'),length)==1 & rowData(se)$resultOutput.equal == "NOVEL" & rowData(se)$resultOutput.compatible != "non_compatible_match"))
  subsetRCs = table(rowData(se)$GENEID[subset])
  temp = numRCs
  temp[names(subsetRCs)] = temp[names(subsetRCs)]-subsetRCs
  numNonSubsetRCs = temp
  #max number of exons of the read classes for each gene
  numExons = as.numeric(by(rowData(se)$numExons, rowData(se)$GENEID, max))
  
  #is spliced?
  isSpliced = ifelse(numExons>1,1,0)
  
  #highCOnfidence
  highConfidence = by(rowData(se)$confidenceType, rowData(se)$GENEID, function(x){
    ifelse("highConfidenceJunctionReads" %in% x, 1,0)
  })
  highConfidence = as.numeric(highConfidence)

  #percentage of reads that are non-multialigned
  # uniqueReadsProp = by(assays(se)$uniqueReads, rowData(se)$GENEID, sum)
  # uniqueReadsProp = as.numeric(highConfidence)/numReads
  
  features = cbind(numReadsLog,strand_bias, numRCs, numNonSubsetRCs, numExons, isSpliced, highConfidence)
  
  #features = cbind(numReadsLog,strand_bias, numRCs, numNonSubsetRCs, numExons, isSpliced, highConfidence)
  
  return(list(features=features, labels = labels, names = geneIDs))
}

trainGeneModel = function(features, labels, names, plot = NULL, saveFig = T){
  
  if(sum(labels)==length(labels) | sum(labels)==0){return(NULL)}
  
  geneScore = getGeneScore(features, labels, names)
  
  if(!is.null(plot)){
    
    if(saveFig){
      svg(paste0(savePath,'/',plot,"_ROC.svg"), width = 7, height = 7)
      plotGeneModel(features, labels)
      dev.off()
    } else{
      plotGeneModel(features, labels)
    }
    
    if(saveFig){
      svg(paste0(savePath,plot, "_Pres_Sens.svg"), width = 7, height = 7)
    }
    measure = measureit(score = geneScore$model$score, class = geneScore$model$testLabels, measure = c("SENS", "PREC"))
    measure$PREC[1]=1
    base::plot(measure$PREC~measure$SENS, type = "l")
    colours = c("red","orange","yellow","green", "turquoise", "blue","purple", "plum", "gray","black")
    thresholds = c(.9, .8, .7, .6, .5, .4, .3, .2, .1, 0)
    for(i in 1:10){
      abline(v = measure$SENS[sum(unique(geneScore$model$score)>thresholds[i])+1], lwd=1, lty = 2, col = colours[i])
    }
    if(saveFig){
      dev.off()
    }
  }
  
  return(geneScore$model$cvfit)
}

getGeneScore = function(features, labels, names, model = NULL){

  if(!is.null(model)){
    
    score = as.numeric(predict(model, newx = as.matrix(features), s = "lambda.min",type="response"))
    
  } else{
    
    result = calculateScore(cbind(features,labels))
    model = result
    score = as.numeric(predict(result$cvfit, newx = as.matrix(features), s = "lambda.min",type="response"))

  }
  
  names(score) = names
  #scoreByRC = score[rowData(se)$GENEID]
  
  return(list(score=score, model = model))
}

trainTranscriptModel = function(features, labels, plot = NULL, saveFig = T){
  
  result=calculateScore(cbind(features, labels))

  if(!is.null(plot)){
    #print(paste0(savePath,'/',plot,"_ROC.svg"))
    if(saveFig){
      svg(paste0('/',plot,"_ROC.svg"), width = 7, height = 7)
      plotTranscriptModel(features, labels, plot)
      dev.off()
    } else{
      plotTranscriptModel(features, labels, plot)
    }

    if(saveFig){
      svg(paste0(savePath,'/',plot, "_Pres_Sens.svg"), width = 7, height = 7)
    }
    measure = measureit(score = result$score, class = result$testLabels, measure = c("SENS", "PREC"))
    base::plot(measure$PREC~measure$SENS, type = "l")
    #draw lines to where the thresholds are at the moment
    colours = c("red","orange","yellow","green", "turquoise", "blue","purple", "plum", "gray","black")
    thresholds = c(.9, .8, .7, .6, .5, .4, .3, .2, .1, 0)
    for(i in 1:10){
      abline(v = measure$SENS[sum(unique(result$score)>thresholds[i])+1], lwd=1, lty = 2, col = colours[i])
    }
    if(saveFig){
      dev.off()
    }
    
  }
  
  return(result$cvfit)
}

plotGeneModel = function(features, labels){
  colours = c('black','gray','red','blue', 'green', 'purple', 'orange', 'yellow', 'brown', 'pink', 'teal', 'darksalmon', 'darkslategray4', 'deeppink2', 'goldenrod', 'burlywood1')
  legend = NULL
  i=1
  
  result = calculateScore(cbind(features,labels))
  x=rocit(score=result$score, class=result$testLabels)
  #x=rocit(score=score, class=labels)
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
  
  plotLine(features[,c("numReadsLog","numNonSubsetRCs")],labels, colours[i])
  i = i+1
  legend = c(legend,"numNonSubsetRCs")
  plotLine(features[,c("numReadsLog","numExons")],labels, colours[i])
  i = i+1
  legend = c(legend,"numExons")
  plotLine(features[,c("numReadsLog","isSpliced")],labels, colours[i])
  i = i+1
  legend = c(legend,"isSpliced")
  plotLine(features[,c("numReadsLog","strand_bias")],labels, colours[i])
  i = i+1
  legend = c(legend,"strand_bias")
  # plotLine(features[,c("numReadsLog","uniqueReadsProp")],labels, colours[i])
  # i = i+1
  # legend = c(legend,"uniqueReadsProp")
  plotLine(features[,c("numReadsLog","numRCs")],labels, colours[i])
  i = i+1
  legend = c(legend,"numRCs")
  
  legend("bottomright", legend=legend, col = colours[1:i],lty=1, cex=1)
}

plotTranscriptModel = function(features, labels, plot){
  legend = NULL
  i= 1
  colours = c('black','gray','red','blue', 'green', 'purple', 'orange', 'yellow', 'brown', 'pink', 'darksalmon', 'darkslategray4', 'deeppink2', 'goldenrod', 'burlywood1')
  
  result = calculateScore(cbind(features, labels))
  x=rocit(score=result$score, class=result$testLabels)
  base::plot(x$FPR, x$TPR, col=colours[i],  main="Transcript Score ROC", type ='l')
  text(0.4,0.2,c(paste0("Trained on:",as.character(result$trainCount))), cex=1.5)
  text(0.4,0.1,c(paste0("Tested on:",as.character(length(result$testLabels)))), cex=1.5)
  text(0.8,0.6,c(paste0("All-AUC:",as.character(signif(x$AUC,3)))), cex=1.5)
  legend = c(legend, c("All"))
  i = i+1
  
  x=rocit(score=result$testData[,1], class=result$testLabels)
  lines(x$FPR, x$TPR, col=colours[i],  main="reads")
  text(0.8,0.5,c(paste0("Read-AUC:",as.character(signif(x$AUC,3)))), cex=1.5)
  legend = c(legend, c("Number of Reads"))
  i = i+1
  
  lines(c(0,1),c(0,1))
  
  if(T){
    plotLine(features[,c("numReads","geneReadProp")],labels, colours[i])
    legend = c(legend,"geneReadProp")
    i = i+1
    plotLine(features[,c("numReads","tx_strand_bias")],labels, colours[i])
    legend = c(legend,"tx_strand_bias")
    i = i+1
    plotLine(features[,c("numReads","SD")],labels, colours[i])
    legend = c(legend,"SD")
    i = i+1
    plotLine(features[,c("numReads","SDend")],labels, colours[i])
    legend = c(legend,"SDend")
    i = i+1
    # plotLine(features[,c("numReads","uniqueReads")],labels, colours[i])
    # legend = c(legend,"uniqueReads")
    # i = i+1
    plotLine(features[,c("numReads","numAstart")],labels, colours[i])
    legend = c(legend,"numAstart")
    i = i+1
    plotLine(features[,c("numReads","numAend")],labels, colours[i])
    legend = c(legend,"numAend")
    i = i+1
    plotLine(features[,c("numReads","alignedPolyAstart")],labels, colours[i])
    legend = c(legend,"alignedPolyAstart")
    i = i+1
    plotLine(features[,c("numReads","alignedPolyAend")],labels, colours[i])
    legend = c(legend,"alignedPolyAend")
    i = i+1
    plotLine(features[,c("numReads","numTstart")],labels, colours[i])
    legend = c(legend,"numTstart")
    i = i+1
    plotLine(features[,c("numReads","numTend")],labels, colours[i])
    legend = c(legend,"numTend")
    i = i+1
    # plotLine(features[,c("numReads","alignedPolyTstart")],labels, colours[i])
    # legend = c(legend,"alignedPolyTstart")
    # i = i+1
    # plotLine(features[,c("numReads","alignedPolyTend")],labels, colours[i])
    # legend = c(legend,"alignedPolyTend")
    # i = i+1
    # plotLine(features[,c("numReads","bothAdapters")],labels, colours[i])
    # legend = c(legend,"bothAdapters")
    # i = i+1
    # plotLine(features[,c("numReads","polyAEnd")],labels, colours[i])
    # legend = c(legend,"polyAEnd")
    # i = i+1
    # plotLine(features[,c("numReads","bothAdaptersProp")],labels, colours[i])
    # legend = c(legend,"bothAdaptersProp")
    # i = i+1
    # plotLine(features[,c("numReads","polyAEndProp")],labels, colours[i])
    # legend = c(legend,"polyAEndProp")
    # i = i+1
    # plotLine(features[,c("numReads","TSSscore")],labels, colours[i])
    # legend = c(legend,"TSSscore")
    # i = i+1
    # plotLine(features[,c("numReads","TESscore")],labels, colours[i])
    # legend = c(legend,"TESscore")
    # i = i+1
    # plotLine(features[,c("numReads","bothTSSandTES")],labels, colours[i])
    # legend = c(legend,"bothTSSandTES")
    # i = i+1
    plotLine(features[,c("numReads","transcriptProp")],labels, colours[i])
    legend = c(legend,"transcriptProp")
    i = i+1
    plotLine(features[,c("numReads","transcriptPropMin")],labels, colours[i])
    legend = c(legend,"transcriptPropMin")
    i = i+1
  }
  if(F){
    plotLine(features[,"logNumReads"],labels, colours[i])
    legend = c(legend,"logNumReads")
    i = i+1
    plotLine(features[,"geneReadProp"],labels, colours[i])
    legend = c(legend,"geneReadProp")
    i = i+1
    plotLine(features[,"tx_strand_bias"],labels, colours[i])
    legend = c(legend,"tx_strand_bias")
    i = i+1
    plotLine(features[,"SD"],labels, colours[i])
    legend = c(legend,"SD")
    i = i+1
    plotLine(features[,"SDend"],labels, colours[i])
    legend = c(legend,"SDend")
    i = i+1
    # plotLine(features[,"uniqueReads"],labels, colours[i])
    # legend = c(legend,"uniqueReads")
    # i = i+1
    plotLine(features[,"numAstart"],labels, colours[i])
    legend = c(legend,"numAstart")
    i = i+1
    plotLine(features[,"numAend"],labels, colours[i])
    legend = c(legend,"numAend")
    i = i+1
    plotLine(features[,"alignedPolyAstart"],labels, colours[i])
    legend = c(legend,"alignedPolyAstart")
    i = i+1
    plotLine(features[,"alignedPolyAend"],labels, colours[i])
    legend = c(legend,"alignedPolyAend")
    i = i+1
    plotLine(features[,"numTstart"],labels, colours[i])
    legend = c(legend,"numTstart")
    i = i+1
    plotLine(features[,"numTend"],labels, colours[i])
    legend = c(legend,"numTend")
    i = i+1
    # plotLine(features[,"alignedPolyTstart"],labels, colours[i])
    # legend = c(legend,"alignedPolyTstart")
    # i = i+1
    # plotLine(features[,"alignedPolyTend"],labels, colours[i])
    # legend = c(legend,"alignedPolyTend")
    # i = i+1
    # plotLine(features[,"bothAdapters"],labels, colours[i])
    # legend = c(legend,"bothAdapters")
    # i = i+1
    # plotLine(features[,"polyAEnd"],labels, colours[i])
    # legend = c(legend,"polyAEnd")
    # i = i+1
    # plotLine(features[,"bothAdaptersProp"],labels, colours[i])
    # legend = c(legend,"bothAdaptersProp")
    # i = i+1
    # plotLine(features[,"polyAEndProp"],labels, colours[i])
    # legend = c(legend,"polyAEndProp")
    # i = i+1
    # plotLine(features[,"TSSscore"],labels, colours[i])
    # legend = c(legend,"TSSscore")
    # i = i+1
    # plotLine(features[,"TESscore"],labels, colours[i])
    # legend = c(legend,"TESscore")
    # i = i+1
    # plotLine(features[,"bothTSSandTES"],labels, colours[i])
    # legend = c(legend,"bothTSSandTES")
    # i = i+1
    plotLine(features[,"transcriptProp"],labels, colours[i])
    legend = c(legend,"transcriptProp")
    i = i+1
    plotLine(features[,"transcriptPropMin"],labels, colours[i])
    legend = c(legend,"transcriptPropMin")
    i = i+1
  }
  
  
  legend("bottomright", legend=legend, col = colours[1:i],lty=1, cex=1)
  
}

plotLine = function(features, labels, colour, log = F){
  if(!is.matrix(features)){
    score = features
  } else {
    result = calculateScore(cbind(features,labels), log = log)
    score = result$score
    labels = result$testLabels
  }
  #score = as.numeric(predict(result$cvfit, newx = result$testData, s = "lambda.min",type="response"))
  x=rocit(score=score, class=labels)
  lines(x$FPR, x$TPR, col=colour,  main="reads")
}

calculateScore = function(input, finalFeatureIndex, model = NULL, log = F){
  finalFeatureIndex = ncol(input)-1
  sampleSizeTraining=0
  sampleSizeValidation=0
  if(is.null(model)){
    #split data for training
    fractionTraining   <- 0.60
    fractionValidation <- 0.20
    fractionTest       <- 0.20
    
    # Compute sample sizes.
    sampleSizeTraining   <- floor(fractionTraining   * nrow(input))
    sampleSizeValidation <- floor(fractionValidation * nrow(input))
    sampleSizeTest       <- floor(fractionTest       * nrow(input))
    
    # Create the randomly-sampled indices for the dataframe. Use setdiff() to
    # avoid overlapping subsets of indices.
    indicesTraining    <- sort(sample(seq_len(nrow(input)), size=sampleSizeTraining))
    indicesNotTraining <- setdiff(seq_len(nrow(input)), indicesTraining)
    indicesValidation  <- sort(sample(indicesNotTraining, size=sampleSizeValidation))
    indicesTest        <- setdiff(indicesNotTraining, indicesValidation)
    
    # Finally, output the three dataframes for training, validation and test.
    dfTraining   <- input[indicesTraining, ]
    dfValidation <- input[indicesValidation, ]
    dfTest       <- input[indicesTest, ]
    
    #fit=glmnet(dfTraining[,1:finalFeatureIndex],dfTraining[,finalFeatureIndex+1])
    if(!log){
      cvfit = cv.glmnet(as.matrix(rbind(dfTraining[,1:finalFeatureIndex],dfValidation[,1:finalFeatureIndex])), c(dfTraining[,finalFeatureIndex+1],dfValidation[,finalFeatureIndex+1]),family = "binomial", type.measure = "class")
    } else {
      cvfit = glmnet(as.matrix(rbind(dfTraining[,1:finalFeatureIndex],dfValidation[,1:finalFeatureIndex])), c(dfTraining[,finalFeatureIndex+1],dfValidation[,finalFeatureIndex+1]),family = "binomial", type.measure = "class")
    }
    #cvfit = cv.glmnet(as.matrix(rbind(dfTraining[,1:finalFeatureIndex],dfValidation[,1:finalFeatureIndex])), c(dfTraining[,finalFeatureIndex+1],dfValidation[,finalFeatureIndex+1]))
  } else{
    dfTest = input
    cvfit = model
  }
  
  if(!log){
    score = as.numeric(predict(cvfit, newx = as.matrix(dfTest[,1:finalFeatureIndex]), s="lambda.min",type="response"))
  } else {
    score = as.numeric(predict(cvfit, newx = as.matrix(dfTest[,1:finalFeatureIndex]),s = 0.01, type="response"))
  }
  #score = as.numeric(predict(cvfit, newx = as.matrix(dfTest[,1:finalFeatureIndex]), s="lambda.min",type="response"))
  testLabels = dfTest[,finalFeatureIndex+1]
  testData = dfTest[,1:finalFeatureIndex]
  
  return(list(score=score, cvfit=cvfit, testData = testData, testLabels=testLabels, trainCount = sampleSizeTraining + sampleSizeValidation, indicesTest = indicesTest))
}
