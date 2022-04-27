#' Assigns each read class geneScore and txScore
#' @param se summerized experiment object with read classes/ranges
#' @param genomeSequence genomeSequence
#' @param annotations GRangesList of annotations
scoreReadClasses = function(se, genomeSequence, annotations, defaultModels, 
                            fit = TRUE, returnModel = FALSE, 
                            min.readCount = 2, verbose = FALSE){
    start.ptm <- proc.time()
    options(scipen = 999) #maintain numeric basepair locations not sci.notfi.
    geneIds = assignGeneIds(rowRanges(se), annotations)
    rowData(se)$GENEID = geneIds$GENEID
    rowData(se)$novelGene = geneIds$novelGene
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
    
    model = NULL
    if (fit) model = bambu.train(se)
    if(returnModel) metadata(se)$model = model
    txScore = getTranscriptScore(rowData(se)[thresholdIndex,], model,
                                 defaultModels, fit=fit)
    rowData(se)$txScore = rep(NA,nrow(se))

    if(!is.null(txScore))  rowData(se)$txScore[thresholdIndex] = txScore

    #calculate using the default model for NDR recommendation
    txScore.noFit = getTranscriptScore(rowData(se)[thresholdIndex,], 
                                model, defaultModels, fit=FALSE)
    rowData(se)$txScore.noFit = rep(NA,nrow(se))
    rowData(se)$txScore.noFit[thresholdIndex] = txScore.noFit

    defaultSuitabilityScore = mean(rowData(se)$txScore.noFit[rowData(se)$equal & rowData(se)$readCount
                           >=min.readCount] - 
        rowData(se)$txScore[rowData(se)$equal & rowData(se)$readCount
                           >=min.readCount])
    if(defaultSuitabilityScore > 0.1) {
        message("Model agreement score: ", defaultSuitabilityScore)
        message("Pre-trained model is performing better than the fitted model")
        message("We recommend running it with fit=FALSE, or train a new model")
    }

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
    strand(start)[which(strand(start)=='*')] = "+"
    strand(end)[which(strand(end)=='*')] = "+"
    startSeqs = BSgenome::getSeq(genomeSequence,start)
    endSeqs = BSgenome::getSeq(genomeSequence,end)
    numATstart = BSgenome::letterFrequency(startSeqs, c("A","T"))
    numATend= BSgenome::letterFrequency(endSeqs, c("A","T"))
    return(data.frame(numAstart=numATstart[,"A"], numAend= numATend[,"A"],  
                      numTstart=numATstart[,"T"], numTend=numATend[,"T"]))
}


#' calculates a score based on how likely a read class is full length
getTranscriptScore = function(rowData, model = NULL, defaultModels, nrounds = 50, fit = TRUE){
    txFeatures = prepareTranscriptModelFeatures(rowData)
    features = dplyr::select(txFeatures,!c(labels))
    if(checkFeatures(txFeatures) & fit & !is.null(model)){
        ## Multi-Exon
        indexME = which(!rowData$novelGene & rowData$numExons>1)
        if(length(indexME)>0){
            txScore = predict(model$transcriptModelME, as.matrix(features))
        } else txScore = NULL
        
        ## Single-Exon
        indexSE = which(!rowData$novelGene & rowData$numExons==1)
        if(length(indexSE)>0){
            txScoreSE = predict(model$transcriptModelSE, as.matrix(features))
        } else txScoreSE = NULL
        
    } else {
        if (!is.null(defaultModels)){
            txScore = predict(defaultModels$transcriptModelME, 
                as.matrix(features))
            txScoreSE = predict(defaultModels$transcriptModelSE, 
                as.matrix(features))
        } else {
            warning("Transcript model not trained. ",
                "No pre-trained models provided. ",
                "Scores will not be calculated and ",
                "transcript discovery will not happen")
            txScore = rep(0, nrow(features))
            txScoreSE = rep(0, nrow(features))
        }
    }
    txScore[which(rowData$numExons==1)] =
        txScoreSE[which(rowData$numExons==1)]
    return(txScore)
}

#' Function to train a model for use on other data
#' @noRd
bambu.train <- function(rcFile = NULL, min.readCount = 2, nrounds = 50, NDR.threshold = 0.1) {
    rowData = rowData(rcFile)[which(rowData(rcFile)$readCount>=min.readCount),]
    txFeatures = prepareTranscriptModelFeatures(rowData)
    features = dplyr::select(txFeatures,!c(labels))
    if(!checkFeatures(txFeatures)){
        message("Transcript model not trained. Using pre-trained models")
        return(NULL)
    }
    ## Multi-Exon
    indexME = which(!rowData$novelGene & rowData$numExons>1)
    if(length(indexME)>0){
        transcriptModelME = fitXGBoostModel(
            data.train=as.matrix(features[indexME,]),
            labels.train=txFeatures$labels[indexME], 
            nrounds = nrounds, show.cv=FALSE)
        txScore = predict(transcriptModelME, as.matrix(features))[indexME]
    }
    ## Single-Exon
    indexSE = which(!rowData$novelGene & rowData$numExons==1)
    if(length(indexSE)>0){
        transcriptModelSE = fitXGBoostModel(
            data.train=as.matrix(features[indexSE,]),
            labels.train=txFeatures$labels[indexSE], 
            nrounds = nrounds, show.cv=FALSE)
        txScoreSE = predict(transcriptModelSE, as.matrix(features))[indexSE]
    }
    ##Calculate the txScore baseline
    NDR = calculateNDR(txScore, txFeatures$labels[indexME])
    NDR.SE = calculateNDR(txScoreSE, txFeatures$labels[indexSE])
    #lm of NDR vs txScore
    lmNDR = lm(txScore~poly(NDR,3,raw=TRUE))
    txScoreBaseline = predict(lmNDR, newdata=data.frame(NDR=NDR.threshold))
    lmNDR.SE = lm(txScoreSE~NDR.SE)
    txScoreBaselineSE = predict(lmNDR.SE, newdata=data.frame(NDR.SE=NDR.threshold))

    return(list(transcriptModelME = transcriptModelME, 
                transcriptModelSE = transcriptModelSE,
                txScoreBaseline = txScoreBaseline,
                txScoreBaselineSE = txScoreBaselineSE,
                lmNDR = lmNDR,
                lmNDR.SE = lmNDR.SE))
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

#' ensures that the data is trainable after filtering
checkFeatures = function(features){
    labels = features$labels
    trainable = TRUE
    if(sum(labels)==length(labels) | sum(labels)==0){
        message("Missing presence of both TRUE and FALSE labels.")
        trainable = FALSE
    }
    if(length(labels)<1000){
        message("Not enough data points")
        trainable = FALSE
    }
    if(sum(labels)<50 | sum(!labels)<50){
        message("Not enough TRUE/FALSE labels")
        trainable = FALSE
    }
    return(trainable)
}
