#' Assigns each read class geneScore and txScore
#' @param se summerized experiment object with read classes/ranges
#' @param genomeSequence genomeSequence
#' @param annotations GRangesList of annotations
scoreReadClasses = function(se, genomeSequence, annotations, defaultModels, 
                            fit = TRUE, min.readCount = 2, verbose = FALSE){
    start.ptm <- proc.time()
    options(scipen = 999) #maintain numeric basepair locations not sci.notfi.
    rowData(se)$GENEID = assignGeneIds(rowRanges(se), annotations)
    rowData(se)$novel = grepl("gene.", rowData(se)$GENEID)
    rowData(se)$numExons = unname(elementNROWS(rowRanges(se)))
    countsTBL = calculateGeneProportion(counts=mcols(se)$readCount,
                                        geneIds=mcols(se)$GENEID)
    rowData(se)$geneReadProp = countsTBL$geneReadProp
    rowData(se)$geneReadCount = countsTBL$geneReadCount

    thresholdIndex = which(rowData(se)$readCount
                        >=min.readCount)
    compTable <- isReadClassCompatible(rowRanges(se[thresholdIndex,]), 
                                    annotations)
    polyATerminals = countPolyATerminals(rowRanges(se[thresholdIndex,]), 
                                        genomeSequence)
    newRowData = data.frame(equal = compTable$equal,
                        compatible = compTable$compatible,
                        numAstart = polyATerminals$numAstart,
                        numAend = polyATerminals$numAend,
                        numTstart = polyATerminals$numTstart,
                        numTend = polyATerminals$numTend)
    rowData(se)[names(newRowData)] = NA
    rowData(se)[thresholdIndex,names(newRowData)] = newRowData

    geneScore = getGeneScore(rowData(se)[thresholdIndex,], defaultModels, fit)
    rowData(se)$geneScore = rep(NA,nrow(se))
    rowData(se)$geneFDR = rep(NA,nrow(se))
    rowData(se)$geneScore[thresholdIndex] = geneScore$geneScore
    rowData(se)$geneFDR[thresholdIndex] = geneScore$geneFDR
    txScore = getTranscriptScore(rowData(se)[thresholdIndex,], 
                                 defaultModels, fit)
    rowData(se)$txScore = rep(NA,nrow(se))
    rowData(se)$txFDR = rep(NA,nrow(se))
    rowData(se)$txScore[thresholdIndex] = txScore$txScore
    rowData(se)$txFDR[thresholdIndex] = txScore$txFDR
    end.ptm <- proc.time()
    if (verbose) 
        message("Finished generating scores for read classes in ", 
        round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
    return(se)
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
    outData <- data.frame(compatible=rep(0, length(query)), 
        equal = rep(FALSE, length(query)))
    query <- cutStartEndFromGrangesList(query)
    subject <- cutStartEndFromGrangesList(subject)

    # reduce memory and speed footprint by reducing number of queries
    # based on all intron match prefilter
    unlistIntronsQuery <- unlistIntrons(query, use.names = FALSE, 
        use.ids = FALSE)
    intronMatchesQuery <- unlistIntronsQuery %in% unlistIntrons(subject,
                        use.names = FALSE, use.ids = FALSE)

    partitioningQuery <- 
        PartitioningByEnd(cumsum(elementNROWS(gaps(ranges(query)))),
        names = NULL)
    allIntronMatchQuery <- all(relist(intronMatchesQuery, partitioningQuery))

    olap = findOverlaps(query[allIntronMatchQuery],subject, 
    ignore.strand = FALSE, type = 'within')
    query <- query[allIntronMatchQuery][queryHits(olap)]

    subject <- subject[subjectHits(olap)]
    splice <- myGaps(query)

    comp <- myCompatibleTranscription(query = query, subject = subject,
    splice = splice)
    equal <- elementNROWS(query)==elementNROWS(subject) & comp

    outData$compatible[allIntronMatchQuery] <- countQueryHits(olap[comp])
    outData$equal[allIntronMatchQuery] <- countQueryHits(olap[equal])>0

    return(outData)
}

#' returns number of A/T's each read class aligned 5' and 3' end
countPolyATerminals = function(grl, genomeSequence){
    start <- resize(granges(unlist(selectStartExonsFromGrangesList(grl, 
        exonNumber = 1),use.names = FALSE)), 
        width = 10, fix = 'start', ignore.strand=FALSE)
    end <- resize(granges(unlist(selectEndExonsFromGrangesList(grl, 
        exonNumber = 1), use.names = FALSE)), 
        width = 10, fix = 'end', ignore.strand=FALSE)
    startSeqs = BSgenome::getSeq(genomeSequence,start)
    endSeqs = BSgenome::getSeq(genomeSequence,end)
    numATstart = letterFrequency(startSeqs, c("A","T"))
    numATend= letterFrequency(endSeqs, c("A","T"))
    return(data.frame(numAstart=numATstart[,"A"], numAend= numATend[,"A"],  
                    numTstart=numATstart[,"T"], numTend=numATend[,"T"]))
}

#' calculates a score based on how likely the read class is associated with a 
#' real gene
getGeneScore = function(rowData, defaultModels, fit = TRUE){
    geneFeatures = prepareGeneModelFeatures(rowData)
    features = dplyr::select(geneFeatures,!c(labels, GENEID))
    if(checkFeatures(geneFeatures) & fit){
        geneModel = fitXGBoostModel(labels.train=geneFeatures$labels,
        data.train=as.matrix(features), show.cv=FALSE)
        geneScore = as.numeric(predict(geneModel, as.matrix(features)))
        geneFDR = calculateFDR(geneScore, geneFeatures$labels)
    } else {
        message("Gene Model not trained. Using prior models")
        geneScore = as.numeric(predict(defaultModels$geneModel, 
                                       as.matrix(features)))
        geneFDR = 1-geneScore
    }
    geneRCMap = match(rowData$GENEID, geneFeatures$GENEID)
    geneScore = geneScore[geneRCMap]
    geneFDR = geneFDR[geneRCMap]
    return(data.frame(geneScore = geneScore, geneFDR = geneFDR))   
}

#' calculates the minimum FDR for each score 
calculateFDR = function(score, labels){
    scoreOrder = order(score, decreasing = TRUE)
    labels = labels[scoreOrder]
    score = score[scoreOrder]
    FDR = cumsum(!labels)/(seq_len(length(score)))
    FDR = rev(cummin(rev(FDR)))
    return(FDR[order(scoreOrder)])
}

#' calculate and format features by gene for model
prepareGeneModelFeatures = function(rowData){
    geneReadCount = NA
    scalingFactor = sum(rowData$readCount)/1000000
    outData <- as_tibble(rowData) %>% group_by(GENEID) %>% 
    summarise(numReads = geneReadCount[1],
        labels = !novel[1], 
        strand_bias=1-abs(0.5-(sum(readCount.posStrand, na.rm=TRUE)/numReads)),
        numRCs=n(), 
        numExons = max(numExons, na.rm=TRUE), 
        isSpliced = numExons>1, 
        # subsets are compatible with 2 or more RCs
        # or compat with only 1 but are not equal
        numNonSubsetRCs = numRCs-sum(compatible>=2 | (compatible==1 & !equal)),
        highConfidence=any(confidenceType=='highConfidenceJunctionReads')) %>%
    mutate(numReads = log2(pmax(1,1+(numReads/scalingFactor))))
    return(outData)
}

#' ensures that the data is trainable after filtering
checkFeatures = function(features){
    labels = features$labels
    if(sum(labels)==length(labels) | sum(labels)==0){
    message("Missing presence of both TRUE and FALSE labels.")
    return(FALSE)
    }
    if(length(labels)<1000){
        message("Not enough data points")
        return(FALSE)
    }
    return(TRUE)
}

#' calculates a score based on how likely a read class is full length
getTranscriptScore = function(rowData, defaultModels, fit = TRUE){
    txFeatures = prepareTranscriptModelFeatures(rowData)
    features = dplyr::select(txFeatures,!c(labels))
    if(checkFeatures(txFeatures) & fit){
        ## Multi-Exon
        indexME = which(!rowData$novel & rowData$numExons>1)
        transcriptModelME = fitXGBoostModel(
                    data.train=as.matrix(features[indexME,]),
                    labels.train=txFeatures$labels[indexME], show.cv=FALSE)
        txScoreME = predict(transcriptModelME, as.matrix(features))

        ## Single-Exon
        indexSE = which(!rowData$novel & rowData$numExons==1)
        transcriptModelSE = fitXGBoostModel(
                    data.train=as.matrix(features[indexSE,]),
                    labels.train=txFeatures$labels[indexSE], show.cv=FALSE)
        txScoreSE = predict(transcriptModelSE, as.matrix(features))
        txScore = txScoreME
        txScore[which(rowData$numExons==1)] =
            txScoreSE[which(rowData$numExons==1)]
        txFDR = calculateFDR(txScore, txFeatures$labels)
    } else {
        message("Transcript model not trained. Using prior models")
        txScoreME = predict(defaultModels$txModel.dcDNA.ME, as.matrix(features))
        txScoreSE = predict(defaultModels$txModel.dcDNA.SE, as.matrix(features))
        txScore = txScoreME
        txScore[which(rowData$numExons==1)] =
            txScoreSE[which(rowData$numExons==1)]
        txFDR = 1-txScore
    }
    return(data.frame(txScore = txScore, txFDR = txFDR))
}

#' calculate and format read class features for model training
prepareTranscriptModelFeatures = function(rowData){
    scalingFactor = sum(rowData$readCount)/1000000
    outData <- as_tibble(rowData) %>%  
    dplyr::select(numReads = readCount, geneReadProp, startSD, endSD,
        numAstart, numAend, numTstart,numTend, 
        tx_strand_bias = readCount.posStrand, labels = equal) %>%
    mutate(
        tx_strand_bias=(1-abs(0.5-(tx_strand_bias/numReads))),
        numReads = log2(pmax(1,1+(numReads/scalingFactor)))
        )
    return(outData)
}