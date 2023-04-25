#' Assigns each read class geneScore and txScore
#' @param se summerized experiment object with read classes/ranges
#' @param genomeSequence genomeSequence
#' @param annotations GRangesList of annotations
#' @noRd
scoreReadClasses = function(se, genomeSequence, annotations, defaultModels, 
                            fit = TRUE, returnModel = FALSE, 
                            min.readCount = 2, min.exonOverlap = 10,
                            fusionMode = FALSE, verbose = FALSE){
    start.ptm <- proc.time()
    options(scipen = 999) #maintain numeric basepair locations not sci.notfi.
    geneIds = assignGeneIds(rowRanges(se), annotations, min.exonOverlap, fusionMode)
    rowData(se)$GENEID = geneIds$GENEID
    rowData(se)$novelGene = geneIds$novelGene
    rowData(se)$numExons = unname(elementNROWS(rowRanges(se)))
    countsTBL = calculateGeneProportion(counts=mcols(se)$readCount,
                                        geneIds=mcols(se)$GENEID)
    rowData(se)$geneReadProp = countsTBL$geneReadProp
    rowData(se)$geneReadCount = countsTBL$geneReadCount
    
    thresholdIndex = which(rowData(se)$readCount
                           >=min.readCount)
    if(length(thresholdIndex)==0){
        warningText = "No read classes with more than 1 read. Unable to train or score any."
        metadata(se)$warnings = c(metadata(se)$warnings, warningText)
        if(verbose) warning(warningText)
    }

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
    
    #calculate using the pretrained model for NDR recommendation
    rowData(se)$txScore.noFit = rep(NA,nrow(se))
    if(length(thresholdIndex)>0){
        txScore.noFit = getTranscriptScore(rowData(se)[thresholdIndex,], 
                                    model = NULL, defaultModels)
        
        rowData(se)$txScore.noFit[thresholdIndex] = txScore.noFit
    }

    model = NULL
    rowData(se)$txScore = rowData(se)$txScore.noFit
    if (fit & length(thresholdIndex)>0){ 
        model = trainBambu(se, verbose = verbose, min.readCount = min.readCount)
        if(returnModel) metadata(se)$model = model
        txScore = getTranscriptScore(rowData(se)[thresholdIndex,], model,
                                 defaultModels)
        rowData(se)$txScore = rep(NA,nrow(se))
        if(!is.null(txScore))  rowData(se)$txScore[thresholdIndex] = txScore
    }
    if(is.null(model) & fit) {
        warningText = "Bambu was unable to train a model on this sample, and is using a pretrained model"
        metadata(se)$warnings = c(metadata(se)$warnings, warningText)
        if(verbose) warning(warningText)
    }
    end.ptm <- proc.time()
    if (verbose) 
        message("Finished generating scores for read classes in ", 
                round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
    return(se)
}

#' % of a genes read counts assigned to each read class
#' @noRd
calculateGeneProportion = function(counts, geneIds){
    countsTBL <- tibble(counts, geneIds) %>%
        group_by(geneIds) %>% mutate(geneReadCount = sum(counts),
                                     geneReadProp = counts/geneReadCount)
    return(countsTBL)
    
}

#' returns number of ref anno each read class is a subset of
#' @noRd
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
#' @noRd
countPolyATerminals = function(grl, genomeSequence){
    start <- resize(granges(unlist(selectStartExonsFromGrangesList(grl, 
            exonNumber = 1),use.names = FALSE)), 
                    width = 10, fix = 'start', ignore.strand=FALSE)
    end <- resize(granges(unlist(selectEndExonsFromGrangesList(grl, 
            exonNumber = 1), use.names = FALSE)), 
                  width = 10, fix = 'end', ignore.strand=FALSE)
    strand(start)[which(strand(start)=='*')] = "+"
    strand(end)[which(strand(end)=='*')] = "+"
    seqlevels(start) = seqlevels(genomeSequence) #needed for windows DNAStringSet
    seqlevels(end) = seqlevels(genomeSequence)
    startSeqs = BSgenome::getSeq(genomeSequence,start)
    endSeqs = BSgenome::getSeq(genomeSequence,end)
    numATstart = BSgenome::letterFrequency(startSeqs, c("A","T"))
    numATend= BSgenome::letterFrequency(endSeqs, c("A","T"))
    return(data.frame(numAstart=numATstart[,"A"], numAend= numATend[,"A"],  
                      numTstart=numATstart[,"T"], numTend=numATend[,"T"]))
}


#' calculates a score based on how likely a read class is full length
#' @noRd
getTranscriptScore = function(rowData, model = NULL, defaultModels){
    txFeatures = prepareTranscriptModelFeatures(rowData)
    features = dplyr::select(txFeatures,!c(labels))
    if(!is.null(model)){
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
#' @title Function to train a model for use on other data
#' @description This function train a model for use on other data
#' @param rcFile summerized experiment object with read classes/ranges produced from bambu(quant = FALSE, discovery = FALSE) or rcOutdir
#' @param min.readCount the minimum number of reads a read class is required to be have to be used for training
#' @param nrounds xgboost hyperparameter - the number of decision trees in the final mode
#' @param NDR.threshold the effective NDR threshold that bambu will try and match on other samples when using this model
#' @param verbose if additional messages should be output
#' Output - A list containing 6 objects which is passed directly into bambu(opt.discovery=list(defaultModels=trainBambu()))
#'      transcriptModelME - the model for multi-exon transcripts 
#'      transcriptModelSE - the model for single-exon transcripts 
#'      txScoreBaseline - the txScore used for NDR calibration for multi-exon transcripts
#'      txScoreBaselineSE - [DEPRECATED] the txScore used for NDR calibration for single-exon transcripts
#'      lmNDR = lmNDR - the linear model of the reletionship between txScore and NDR used to calculate the baseline for multi-exon transcripts
#'      lmNDR.SE = lmNDR.SE - the linear model of the reletionship between txScore and NDR used to calculate the baseline for single-exon transcripts
#'      NDR.threshold - the NDR threshold usd to calculate the txScoreBaseline on the lmNDR (baselineFDR)
#' @details 
#' @return It returns a model object to use in \link{bambu}
#' @export
trainBambu <- function(rcFile = NULL, min.readCount = 2, nrounds = 50, NDR.threshold = 0.1, verbose = TRUE) {
    rowData = rowData(rcFile)[which(rowData(rcFile)$readCount>=min.readCount),]
    txFeatures = prepareTranscriptModelFeatures(rowData)
    features = dplyr::select(txFeatures,!c(labels))
    if(!checkFeatures(txFeatures, verbose)){
        if(verbose) message("Transcript model not trained. Using pre-trained models")
        return(NULL)
    }
    transcriptModelME = NULL
    transcriptModelSE = NULL
    txScoreBaseline = NA
    txScoreBaselineSE = NA
    lmNDR = NULL
    lmNDR.SE = NULL
    ## Multi-Exon
    indexME = which(!rowData$novelGene & rowData$numExons>1)
    if(length(indexME)>0){
        transcriptModelME = fitXGBoostModel(
            data.train=as.matrix(features[indexME,]),
            labels.train=txFeatures$labels[indexME], 
            nrounds = nrounds, show.cv=FALSE)
        txScore = predict(transcriptModelME, as.matrix(features))[indexME]

        ##Calculate the txScore baseline
        NDR = calculateNDR(txScore, txFeatures$labels[indexME])
        #lm of NDR vs txScore
        lmNDR = lm(txScore~poly(NDR,3,raw=TRUE))
        txScoreBaseline = predict(lmNDR, newdata=data.frame(NDR=NDR.threshold))

        ## Compare the trained model AUC to the default model AUC
        txScore.default = predict(defaultModels$transcriptModelME, as.matrix(features))[indexME] 
        newPerformance = evaluatePerformance(txFeatures$labels[indexME],txScore)
        currentPerformance = evaluatePerformance(txFeatures$labels[indexME],txScore.default)
        if(verbose){
        message("On the dataset the new trained model achieves a ROC AUC of ",
            signif(newPerformance$AUC,3),  " and a Precision-Recall AUC of ", signif(newPerformance$PR.AUC,3), ".", 
            "This is compared to the Bambu pretrained model trained on human ONT RNA-Seq data model which achiveves a ROC AUC of ",
            signif(currentPerformance$AUC,3), " and a Precision-Recall AUC of ", signif(currentPerformance$PR.AUC,3))
        }
        #shrink size of lm
        lmNDR = trim_lm(lmNDR)
    }
    ## Single-Exon
    indexSE = which(!rowData$novelGene & rowData$numExons==1)
    if(length(indexSE)>0){
        transcriptModelSE = fitXGBoostModel(
            data.train=as.matrix(features[indexSE,]),
            labels.train=txFeatures$labels[indexSE], 
            nrounds = nrounds, show.cv=FALSE)
        txScoreSE = predict(transcriptModelSE, as.matrix(features))[indexSE]

        NDR.SE = calculateNDR(txScoreSE, txFeatures$labels[indexSE])
        lmNDR.SE = glm(txScoreSE~NDR.SE)
        txScoreBaselineSE = predict(lmNDR.SE, newdata=data.frame(NDR.SE=NDR.threshold))
        lmNDR.SE = trim_lm(lmNDR.SE)
    }


    

    return(list(transcriptModelME = transcriptModelME, 
                transcriptModelSE = transcriptModelSE,
                txScoreBaseline = txScoreBaseline,
                txScoreBaselineSE = txScoreBaselineSE,
                lmNDR = lmNDR,
                lmNDR.SE = lmNDR.SE,
                NDR.threshold = NDR.threshold))
}

#' reduces the size of a lm so it can be saved with a lower footprint for prediction
#' @noRd
trim_lm = function(lm){
    lm$residuals = c()
    lm$effects = c()
    lm$fitted.values = c()
    lm$model = c()
    lm$qr$qr=c()
    attr(lm$terms, ".Environment") = c()
    lm$linear.predictors = NULL
    lm$weights = NULL
    lm$prior.weights = NULL
    lm$y = NULL
    lm$data = NULL
    lm$formula = NULL
    return(lm)
}

#' calculate and format read class features for model training
#' @noRd
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
#' @noRd
checkFeatures = function(features, verbose = FALSE){
    labels = features$labels
    trainable = TRUE
    print(features)
    if(sum(labels)==length(labels) | sum(labels)==0){
        if (verbose) message("Sample is missing presence of both TRUE and FALSE labels.")
        trainable = FALSE
    }
    if(length(labels)<1000){
        if (verbose) message("Sample has less than 1000 labeled read classes")
        trainable = FALSE
    }
    if(sum(labels)<50 | sum(!labels)<50){
        if (verbose) message("Sample does not have more than 50 of both read class labels")
        trainable = FALSE
    }
    return(trainable)
}
