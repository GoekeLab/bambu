#' calculate stranded read counts
#' @param uniqueJunctions uniqueJunctions
#' @param junctionMatchList junctionMatchList
#' @param genomeSequence genomeSequence
#' @param unstranded_unlisted_junctions unstranded_unlisted_junctions
#' @param unlisted_junction_granges unlisted_junction_granges
#' @noRd
calculateStrandedReadCounts <- function(uniqueJunctions,
    genomeSequence,unstranded_unlisted_junctions,
    unlisted_junction_granges) {
    junctionMatchList <- as(findMatches(uniqueJunctions,
        unstranded_unlisted_junctions),"List")
    uniqueJunctions_score <- elementNROWS(junctionMatchList)
    junctionStrandList <- extractList(strand(unlisted_junction_granges),
        junctionMatchList)
    junctionSeqStart <- BSgenome::getSeq(genomeSequence,
        IRanges::shift(flank(uniqueJunctions,width = 2), 2))#shift from IRanges
    junctionSeqEnd <- BSgenome::getSeq(genomeSequence,
        IRanges::shift(flank(uniqueJunctions,width = 2, start = FALSE), -2))
    junctionMotif <- paste(junctionSeqStart, junctionSeqEnd, sep = "-")
    junctionStartName <- paste(seqnames(uniqueJunctions),start(uniqueJunctions),
        sep = ":")
    junctionEndName <- paste(seqnames(uniqueJunctions), end(uniqueJunctions),
        sep = ":")
    startScore <- as.integer(tapply(uniqueJunctions_score,
        junctionStartName, sum)[junctionStartName])
    endScore <- as.integer(tapply(uniqueJunctions_score,
        junctionEndName, sum)[junctionEndName])
    mcols(uniqueJunctions) <- DataFrame(
        score = uniqueJunctions_score,
        plus_score = sum(junctionStrandList == "+"),
        minus_score = sum(junctionStrandList == "-"),
        spliceMotif = junctionMotif,
        spliceStrand = spliceStrand(junctionMotif),
        junctionStartName = junctionStartName,
        junctionEndName = junctionEndName,
        startScore = startScore,
        endScore = endScore,
        id = seq_along(uniqueJunctions))
    strand(uniqueJunctions) <- uniqueJunctions$spliceStrand
    return(uniqueJunctions)
}


#' Create Junction tables from unlisted junction granges
#' @importFrom BiocParallel bppram bpvec
#' @noRd
createJunctionTable <- function(unlisted_junction_granges,
                                genomeSequence = NULL, ncore = 1) {
    # License note: This function is adopted from the GenomicAlignments package 
    if (is.null(genomeSequence)) stop("Reference genome sequence is missing,
        please provide fasta file or BSgenome name, see available.genomes()")
    if (is(genomeSequence, "character")) {
        if (grepl(".fa", genomeSequence)) {
            if (.Platform$OS.type == "windows") {
                genomeSequence <- Biostrings::readDNAStringSet(genomeSequence)
                newlevels <- unlist(lapply(strsplit(names(genomeSequence)," "),
                    "[[", 1))
                names(genomeSequence) <- newlevels
            } else {
                genomeSequence <- Rsamtools::FaFile(genomeSequence)
            }
        } else {
            genomeSequence <- BSgenome::getBSgenome(genomeSequence)
        }
    }
    
    if (!all(seqlevels(unlisted_junction_granges) %in%
            seqlevels(genomeSequence))) {
        message("not all chromosomes present in reference genome sequence,
            ranges are dropped")
        unlisted_junction_granges <- keepSeqlevels(unlisted_junction_granges,
            value = seqlevels(unlisted_junction_granges)[seqlevels(
            unlisted_junction_granges) %in% seqlevels(genomeSequence)],
            pruning.mode = "coarse")
    }
    unstranded_unlisted_junctions <- unstrand(unlisted_junction_granges)
    uniqueJunctions <- sort(unique(unstranded_unlisted_junctions))
    names(uniqueJunctions) <- paste("junc", seq_along(uniqueJunctions),
                                    sep = ".")
    uniqueJunctions <- calculateStrandedReadCounts(uniqueJunctions,
        genomeSequence, unstranded_unlisted_junctions,
        unlisted_junction_granges)
    return(uniqueJunctions)
}

#' update junctions object if strand prediction improves overlap 
#' with annotations
#' @param annotatedIntronNumber annotatedIntronNumber
#' @param uniqueJunctions uniqueJunctions
#' @param uniqueJunctionsUpdate uniqueJunctionsUpdate
#' @param uniqueAnnotatedIntrons uniqueAnnotatedIntrons
#' @param strandStep strandStep
#' @param verbose verbose
#' @noRd
updateJunctionwimprove <- function(annotatedIntronNumber, uniqueJunctions,
    uniqueJunctionsUpdate, uniqueAnnotatedIntrons, strandStep, verbose) {
    annotatedIntronNumberNew <- evalAnnotationOverlap(uniqueJunctionsUpdate,
        uniqueAnnotatedIntrons, ignore.strand = FALSE)["TRUE"]
    if (annotatedIntronNumberNew > annotatedIntronNumber & !is.na(
        annotatedIntronNumber)) {
        # update junctions object if strand prediction improves overlap
        # with annotations
        if (verbose) {
            message("after strand correction, annotated introns:")
            message(annotatedIntronNumberNew)
            message(annotatedIntronNumberNew / length(uniqueJunctionsUpdate))
        }
        annotatedIntronNumber <- annotatedIntronNumberNew
        uniqueJunctions <- uniqueJunctionsUpdate
    } else {
        strandStep <- FALSE
    }
    outputList <- list(
        "strandStep" = strandStep,
        "annotatedIntronNumber" = annotatedIntronNumber,
        "uniqueJunctions" = uniqueJunctions
    )
    return(outputList)    
}

#' JUNCTIONSTRANDCORRECTION
#' @noRd
junctionStrandCorrection <- function(uniqueJunctions, unlisted_junction_granges,
    uniqueAnnotatedIntrons, stranded, verbose = FALSE) {
    # note: strand sometimes incorrectly infered based on motifs, might 
    # introduce systematic errors due to alignment (biased to splice motifs)
    allJunctionToUniqueJunctionOverlap <- 
        findOverlaps(unlisted_junction_granges, uniqueJunctions,
            type = "equal", ignore.strand = TRUE)
    uniqueJunctionsUpdate <- uniqueJunctions
    # make a copy to revert to if strand correction does not improve results
    annotatedIntronNumber <- evalAnnotationOverlap(uniqueJunctions,
        uniqueAnnotatedIntrons,ignore.strand = FALSE)["TRUE"]
    if (verbose) {
        message("before strand correction, annotated introns:")
        message(annotatedIntronNumber)
        message(annotatedIntronNumber / length(uniqueJunctions))
    }
    # infer strand for each read based on strand of junctions
    strandStep <- TRUE
    while (strandStep) { # iterate twice to improve strand prediction w.t.
        # mean junction counts, annotate junction strand with read strand
        if (stranded) { # update junction strand score
            uniqueJunctionsUpdate <- 
                updateStrandScoreByRead(unlisted_junction_granges,
                    uniqueJunctionsUpdate, uniqueJunctions, stranded, 
                    allJunctionToUniqueJunctionOverlap)
            strandScoreByRead <-
                uniqueJunctionsUpdate$minus_score_inferedByRead -
                uniqueJunctionsUpdate$plus_score_inferedByRead
        } else{ 
            strandScoreByRead <- uniqueJunctionsUpdate$minus_score -
                uniqueJunctionsUpdate$plus_score
        }
        # overwrite info from motif which increases overlap with known junc
        strand(uniqueJunctionsUpdate[strandScoreByRead < 0]) <- "+"
        strand(uniqueJunctionsUpdate[strandScoreByRead > 0]) <- "-"
        updatedList <- updateJunctionwimprove(annotatedIntronNumber,
            uniqueJunctions, uniqueJunctionsUpdate, uniqueAnnotatedIntrons,
            strandStep, verbose)
        annotatedIntronNumber <- updatedList$annotatedIntronNumber
        uniqueJunctions <- updatedList$uniqueJunctions
        strandStep <- updatedList$strandStep
    }
    return(list(uniqueJunctions = uniqueJunctions,
                unlisted_junctions = unlisted_junction_granges))
}

#' updateStrandScoreByRead
#' @noRd
updateStrandScoreByRead <- function(unlisted_junction_granges,
    uniqueJunctionsUpdate, uniqueJunctions, stranded, 
    allJunctionToUniqueJunctionOverlap){
    unlisted_junction_granges_strand <-
        as.character(strand(uniqueJunctionsUpdate)[subjectHits(
            allJunctionToUniqueJunctionOverlap)])
    unlisted_junction_granges_strandList <-
        splitAsList(unlisted_junction_granges_strand,
            mcols(unlisted_junction_granges)$id)
    strandJunctionSum <-
        sum(unlisted_junction_granges_strandList == "-") -
        sum(unlisted_junction_granges_strandList == "+")
    readStrand <- rep("*", length(unlisted_junction_granges_strandList))
    readStrand[strandJunctionSum < 0] <- "+"
    readStrand[strandJunctionSum > 0] <- "-"
    strand(unlisted_junction_granges) <-
        readStrand[match(mcols(unlisted_junction_granges)$id,
            as.integer(names(unlisted_junction_granges_strandList)))]
    unstranded_unlisted_junction_granges <-
        unstrand(unlisted_junction_granges)
    junctionMatchList <- as(findMatches(uniqueJunctions,
        unstranded_unlisted_junction_granges), "List")
    tmp <- extractList(strand(unlisted_junction_granges),
        junctionMatchList)
    uniqueJunctionsUpdate$plus_score_inferedByRead <- sum(tmp == "+")
    uniqueJunctionsUpdate$minus_score_inferedByRead <- sum(tmp == "-")
    return(uniqueJunctionsUpdate)
}



#' Evaluate annoation overlap
#' @noRd
evalAnnotationOverlap <- function(intronRanges, uniqueAnnotatedIntrons,
    ignore.strand = FALSE) {
    return(table(!is.na(GenomicRanges::match(intronRanges,
    uniqueAnnotatedIntrons, ignore.strand = ignore.strand ))))
}

#' Test splice sites
#' @noRd
testSpliceSites <- function(data, splice = "Start", prime = "start", 
    junctionModel = NULL, verbose = FALSE){ 
    distSplice.prime <- data[, paste0("dist",splice,".",prime)]
    spliceStrand <- data[, "spliceStrand"]
    spliceStrand.prime <- data[, paste0("spliceStrand.",prime)]
    spliceScore <- data[, paste0(tolower(splice),"Score")]
    spliceScore.prime <-  data[, paste0(tolower(splice),"Score.",prime)]
    annotatedSplice.prime <- data[, paste0("annotated",splice,".",prime)]
    annotatedSplice <- data[, paste0("annotated",splice)]
    predSplice.primeName <- paste0('spliceSitePrediction',splice,'.',prime)
    if(prime == "start"){
        mySet.all <- which((distSplice.prime != 0) & (spliceStrand != "*") &
            (spliceScore > 0) & (distSplice.prime < 15))
    }else{
        mySet.all <- which((distSplice.prime != 0) & (spliceStrand != "*") &
            (spliceScore > 0) & (abs(distSplice.prime) < 15))
    }
    predictionSplice.prime <- rep(NA, nrow(data))
    if (any(mySet.all)) {
        mySet.training = 
            which((annotatedSplice.prime | annotatedSplice)[mySet.all])
        myData = data.frame(spliceScore / (spliceScore + spliceScore.prime),
            spliceScore, distSplice.prime, (spliceStrand.prime == "+"),
            (spliceStrand.prime == "-"), (spliceStrand == "+"))[mySet.all,]
        colnames(myData) <- paste('A',seq_len(ncol(myData)),sep = '.')
        modelmatrix =
            model.matrix(~A.1+A.2+A.3+A.4+A.5, data = data.frame(myData))
        predSplice.prime <- NULL
        if (is.null(junctionModel)) { 
            myResults = fitBinomialModel(labels.train = 
                as.integer(annotatedSplice)[mySet.all][mySet.training], 
                data.train = modelmatrix[mySet.training,],
                data.test = modelmatrix, show.cv = verbose, maxSize.cv = 10000)
            predSplice.prime <- myResults[[2]]
            predictions = myResults[[1]]
        } else{
            predictions = glmnet:::predict.cv.glmnet(
                junctionModel[[predSplice.primeName]],
                newx = modelmatrix,s = 'lambda.min')
        }
        predictionsSplice.prime <- rep(NA, nrow(data))
        names(predictionSplice.prime) <-
            data[, paste0("junction",splice,"Name")]
        predictionSplice.prime[mySet.all] = predictions
        return(list(preds = predictionSplice.prime, predsJ = predSplice.prime))
    } else {
        return(list(preds = rep(NA, nrow(data)), predsJ = NULL))
    }
}

#' Create metadata for splice information
#' @noRd
createSpliceMetadata <- function(annotatedJunctions, splice){
    metadata <- mcols(annotatedJunctions)[,c(paste0(tolower(splice),'Score'),
        paste0('junction',splice,'Name'), paste0('annotated',splice),
        'spliceStrand','spliceMotif')]
    len <- length(annotatedJunctions)
    distdata <- data.frame(
        dist.start = c(0,(diff(start(annotatedJunctions)) *
            as.integer((seqnames(annotatedJunctions[-1]) == 
            seqnames(annotatedJunctions[-len]))))),
        annotated.start = c(FALSE,metadata[,3][-len]),
        Score.start = c(FALSE,metadata[,1][-len]),
        spliceStrand.start = c(FALSE,metadata[,4][-len]),
        spliceMotif.start = c(FALSE,metadata[,5][-len]),
        dist.end = c(-diff(end(annotatedJunctions)) *
            as.integer((seqnames(annotatedJunctions[-len]) == 
            seqnames(annotatedJunctions[-1]))),0),
        annotated.end = c(metadata[,3][-1],FALSE),
        Score.end = c(metadata[,1][-1],FALSE),
        spliceStrand.end = c(metadata[,4][-1],FALSE),
        spliceMotif.end = c(metadata[,5][-1],FALSE))
    colnames(distdata) <- 
        c(paste0('dist',splice,'.start'),
          paste0('annotated',splice,'.start'),
          paste0(tolower(splice),'Score.start'),
          'spliceStrand.start','spliceMotif.start',
          paste0('dist',splice,'.end'), paste0('annotated',splice,'.end'),
          paste0(tolower(splice),'Score.end'),
          'spliceStrand.end', 'spliceMotif.end')
    metadata <- cbind(metadata, distdata)
    return(metadata)
}

#' Predict splicing junctions
#' @noRd
predictSpliceJunctions <- function(annotatedJunctions, junctionModel=NULL,
    verbose = FALSE){
    spliceVec <- c("Start","End")
    ## if needed this can be a single function
    metadataList <- lapply(spliceVec, function(splice){
        annotatedJunctionsTmp <- GRanges(seqnames = seqnames(annotatedJunctions),
            ranges = IRanges(start = get(tolower(splice))(annotatedJunctions),
            end = get(tolower(splice))(annotatedJunctions)), strand = '*')
        mcols(annotatedJunctionsTmp) <- mcols(annotatedJunctions)
        if(splice == "Start"){
            annotatedJunctionsTmp <- unique(annotatedJunctionsTmp)
        }else{
            annotatedJunctionsTmp <- sort(unique(annotatedJunctionsTmp))
        }
        return(createSpliceMetadata(annotatedJunctionsTmp, splice))})
    names(metadataList) <- spliceVec
    if ( is.null(junctionModel)) junctionModelList <- list()
    preds <- lapply(spliceVec, function(splice){
        preds <- lapply(tolower(spliceVec), function(prime){
            return(testSpliceSites(metadataList[[splice]], splice = splice,
                                   prime = prime, junctionModel, verbose))})
        names(preds) <- tolower(spliceVec) 
        return(preds)})
    names(preds) <- spliceVec
    annotatedJunctions$spliceSitePredictionStart.start <- 
        preds[[1]][[1]]$preds[annotatedJunctions$junctionStartName]
    annotatedJunctions$spliceSitePredictionStart.end <- 
        preds[[1]][[2]]$preds[annotatedJunctions$junctionStartName]
    annotatedJunctions$spliceSitePredictionEnd.start <- 
        preds[[2]][[1]]$preds[annotatedJunctions$junctionEndName]
    annotatedJunctions$spliceSitePredictionEnd.end <- 
        preds[[2]][[2]]$preds[annotatedJunctions$junctionEndName]
    if (is.null(junctionModel)) {
        for (splice in spliceVec) {
            for (prime in tolower(spliceVec)) {
                junctionModelList[[paste0('spliceSitePrediction',
                splice,'.',prime)]] <- preds[[splice]][[prime]]$predsJ
            }}}
    if (is.null(junctionModel)) junctionModel <- junctionModelList
    return(list(annotatedJunctions, junctionModel))
}

#' Fit binomial model
#' @noRd
fitBinomialModel <- function(labels.train, data.train, data.test, 
    show.cv=TRUE, maxSize.cv=10000){
    if (show.cv) {
        mySample <- sample(seq_along(labels.train),
            min(floor(length(labels.train)/2),maxSize.cv))
        data.train.cv <- data.train[mySample,]
        labels.train.cv <- labels.train[mySample]
        data.train.cv.test <- data.train[-mySample,]
        labels.train.cv.test <- labels.train[-mySample]
        cv.fit <- glmnet::cv.glmnet(x = data.train.cv,
            y = labels.train.cv, family = 'binomial')
        predictions <-
            glmnet:::predict.cv.glmnet(cv.fit, newx = data.train.cv.test,
            s = 'lambda.min')
        message('prediction accuracy (CV) (higher for splice 
                donor than splice acceptor)')
        testResults <- fisher.test(table(predictions > 0,labels.train.cv.test))
        show(testResults$estimate)
        show(testResults$p.value)
        show(evalutePerformance(labels.train.cv.test == 1,predictions)$AUC)
    }
    cv.fit <- 
        glmnet::cv.glmnet(x = data.train, y = labels.train,family = 'binomial')
    predictions <-
        glmnet:::predict.cv.glmnet(cv.fit, newx = data.test, s = 'lambda.min')
    return(list(predictions,cv.fit))
}

#' if conflict (very rare) use reference junction with higher read count/score
#' @noRd
useRefJunctionForConflict <- function(junctions, candidateJunctionsMinus, 
    candidateJunctionsPlus){
    conflictJunctions <- junctions[names(candidateJunctionsMinus[which(!is.na(
        candidateJunctionsMinus$mergedHighConfJunctionId))][which(
            names(candidateJunctionsMinus)[!is.na(
                candidateJunctionsMinus$mergedHighConfJunctionId)] %in%
                names(candidateJunctionsPlus)[!is.na(
                    candidateJunctionsPlus$mergedHighConfJunctionId)])])]
    conflictNames <- names(conflictJunctions)
    scoreDiff <- junctions[candidateJunctionsPlus[
        conflictNames]$mergedHighConfJunctionId]$score - 
        junctions[candidateJunctionsMinus[
            conflictNames]$mergedHighConfJunctionId]$score
    if (any(scoreDiff < 0)) 
        candidateJunctionsPlus[conflictNames]$mergedHighConfJunctionId[
            scoreDiff < 0 ] <- candidateJunctionsMinus[
                conflictNames]$mergedHighConfJunctionId[ scoreDiff < 0 ]
    if (any(scoreDiff > 0)) 
        candidateJunctionsMinus[conflictNames]$mergedHighConfJunctionId[
            scoreDiff > 0 ] <- candidateJunctionsPlus[
                conflictNames]$mergedHighConfJunctionId[scoreDiff > 0 ]
    if (0 %in% scoreDiff) 
        candidateJunctionsPlus[conflictNames]$mergedHighConfJunctionId[
            scoreDiff == 0] <- candidateJunctionsMinus[
                conflictNames]$mergedHighConfJunctionId[scoreDiff == 0] <- NA
    mergedHighConfJunctionId <- rep(NA, length(junctions))
    names(mergedHighConfJunctionId) <- names(junctions)
    mergedHighConfJunctionId[names(candidateJunctionsPlus)] <- 
        candidateJunctionsPlus$mergedHighConfJunctionId
    mergedHighConfJunctionId[names(candidateJunctionsMinus)] <- 
        candidateJunctionsMinus$mergedHighConfJunctionId
    junctions$mergedHighConfJunctionId <- as.character(mergedHighConfJunctionId)
    return(junctions)
}

#' find junctions by plus and minus strand
#' @noRd
findJunctionsByStrand <- function(candidateJunctions,highConfidentJunctionSet,
    junctionModel, verbose){
    highConfJunctions <- predictSpliceJunctions(candidateJunctions[
        which(highConfidentJunctionSet)], junctionModel = junctionModel)[[1]]
    candidateJunctions$highConfJunctionPrediction = rep(FALSE,
        length(candidateJunctions))
    ## replace all NA's with 0
    spliceSitePredictionList <- 
        cbind(highConfJunctions$spliceSitePredictionStart.start,
              highConfJunctions$spliceSitePredictionStart.end,
              highConfJunctions$spliceSitePredictionEnd.start,
              highConfJunctions$spliceSitePredictionEnd.end)
    spliceSitePredictionList[is.na(spliceSitePredictionList)] <- 2 # NA
    setReferenceJunctions <- (apply(spliceSitePredictionList>0,1,sum)==4)| 
        highConfJunctions$annotatedJunction
    candidateJunctions$highConfJunctionPrediction[highConfidentJunctionSet] <- 
        setReferenceJunctions
    mergedHighConfJunctionId <- rep(NA,length(candidateJunctions))
    ##max Distance can be a parameter that can be set by users
    #here: assign reference junction to all junctions based on prediciton score
    for (maxDist in 0:10) {
        if (verbose) message(maxDist)
        overlapByDist = findOverlaps(candidateJunctions[which(is.na(
            mergedHighConfJunctionId))], candidateJunctions[which(
            candidateJunctions$highConfJunctionPrediction)],
            maxgap = maxDist,type = 'equal')
        mergedHighConfJunctionId[which(is.na(mergedHighConfJunctionId))][
            queryHits(overlapByDist)] <- names(candidateJunctions)[
                which(candidateJunctions$highConfJunctionPrediction)][
                    subjectHits(overlapByDist)]
    }
    candidateJunctions$mergedHighConfJunctionId <- mergedHighConfJunctionId
    return(candidateJunctions)
}

#' this function adds "mergedHighConfJunctionId" to the junciton list
#' which contains the ID of the most likely high confident junction that
#' each junction originates from
#' @noRd
findHighConfidenceJunctions <- function(junctions, junctionModel,
    verbose = FALSE) {
    if (verbose) {
        message('reads count for all annotated junctions')
        message(sum(junctions$score[junctions$annotatedJunction]))
        message(sum(junctions$score[junctions$annotatedJunction]) /
                    sum(junctions$score))
    }
    ##note: the output can be visualised (bed/bigbed track) 
    #calculation is based on distance and properties of next junctions
    candidateJunctionsPlus <-
        junctions[which(strand(junctions) == '+' | strand(junctions) == '*')]
    candidateJunctionsMinus <-
        junctions[which(strand(junctions) == '-' | strand(junctions) == '*')]
    highConfidentJunctionSetPlus <- candidateJunctionsPlus$score > 1 &
        candidateJunctionsPlus$spliceStrand == '+'
    highConfidentJunctionSetMinus <- candidateJunctionsMinus$score > 1 &
        candidateJunctionsMinus$spliceStrand == '-'
    if (sum(highConfidentJunctionSetPlus) > 0 &
        sum(highConfidentJunctionSetMinus) > 0) {
        candidateJunctionsPlus <-
            findJunctionsByStrand(candidateJunctionsPlus, 
                highConfidentJunctionSetPlus, junctionModel, verbose)
        candidateJunctionsMinus <-
            findJunctionsByStrand(candidateJunctionsMinus, 
                highConfidentJunctionSetMinus,junctionModel, verbose)
        
        #if conflict (very rare) use ref junction with higher read count/score
        junctions <- useRefJunctionForConflict(junctions,
                candidateJunctionsMinus, candidateJunctionsPlus)
    } else {
        warning('no junction correction as no high confidence
                reference junctions found')
        junctions$mergedHighConfJunctionId <-  names(junctions)
    }
    rm(candidateJunctionsPlus, candidateJunctionsMinus)
    if (verbose) {
        message('reads count for all annotated junctions after 
                correction to reference junction')
        sumByJuncId <- tapply(junctions$score,
                junctions$mergedHighConfJunctionId, sum)
        message(
            sum(sumByJuncId[junctions[names(sumByJuncId)]$annotatedJunction]))
        message(sum(sumByJuncId[junctions[names(sumByJuncId)]$
            annotatedJunction]) / sum(junctions$score))
    }
    return(junctions[,'mergedHighConfJunctionId'])
}

#' Evaluate performance
#' @noRd
evalutePerformance <- function(labels, scores, decreasing = TRUE){
    labels <- labels[order(scores, decreasing = decreasing)]
    results <- list()
    # TP/(TP+FP); True Positive Rate;Sensitivity; recall
    results[['TPR']] <- cumsum(labels)/sum(labels)
    # FP/(FP+TN); False Positive Rate;1-Specificity
    results[['FPR']] <- cumsum(!labels)/sum(!labels)
    # TP/(TP+FP); positive predictive value;precision
    results[['precision']] <- cumsum(labels)/(seq_along(labels))
    results[['AUC']] <- sum(results[['TPR']][!duplicated( results[['FPR']],
        fromLast = TRUE)] / sum(!duplicated( results[['FPR']],
        fromLast = TRUE)))
    return(results)
}