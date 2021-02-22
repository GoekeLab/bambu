#' correct junction from prediction
#' @param uniqueJunctions uniqueJunctions
#' @param verbose verbose
#' @noRd
junctionErrorCorrection <- function(uniqueJunctions, verbose) {
    start.ptm <- proc.time()
    if (sum(uniqueJunctions$annotatedJunction) > 5000 &
        sum(!uniqueJunctions$annotatedJunction) > 4000) {
        uniqJunctionsNmodels <-
            findUniqueJunctions(uniqueJunctions, NULL, verbose)
        uniqueJunctions <- uniqJunctionsNmodels$uniqueJunctions
        junctionModel <- uniqJunctionsNmodels$predictSpliceSites[[2]]
    } else {
        junctionModel <- standardJunctionModels_temp
        uniqJunctionsNmodels <-
            findUniqueJunctions(uniqueJunctions, junctionModel, verbose)
        uniqueJunctions <- uniqJunctionsNmodels$uniqueJunctions
        message("Junction correction with not enough data,
            precalculated model is used")
    }
    end.ptm <- proc.time()
    if (verbose) 
        message("Model to predict true splice sites built in ",
                round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    start.ptm <- proc.time()
    uniqueJunctions <- findHighConfidenceJunctions(junctions = uniqueJunctions,
        junctionModel = junctionModel, verbose = verbose)
    uniqueJunctions$mergedHighConfJunctionIdAll_noNA <- 
        uniqueJunctions$mergedHighConfJunctionId
    uniqueJunctions$mergedHighConfJunctionIdAll_noNA[
        is.na(uniqueJunctions$mergedHighConfJunctionId)] <- 
        names(uniqueJunctions[is.na(uniqueJunctions$mergedHighConfJunctionId)])
    uniqueJunctions$strand.mergedHighConfJunction <- 
        as.character(strand(
            uniqueJunctions[uniqueJunctions$mergedHighConfJunctionIdAll_noNA]))
    end.ptm <- proc.time()
    if (verbose) 
        message("Finished correcting junction based on set of high confidence
            junctions in ", round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    return(uniqueJunctions)
}


#' find unique junctions
#' @noRd
findUniqueJunctions <- function(uniqueJunctions, junctionModel, verbose){
    predictSpliceSites <- predictSpliceJunctions(
        annotatedJunctions = uniqueJunctions,
        junctionModel = junctionModel,
        verbose = verbose)
    uniqueJunctions <- predictSpliceSites[[1]][, c(
        "score", "spliceMotif",
        "spliceStrand", "junctionStartName", "junctionEndName",
        "startScore", "endScore", "annotatedJunction",
        "annotatedStart", "annotatedEnd")]
    return(list("uniqueJunctions" = uniqueJunctions, 
                "predictSpliceSites" = predictSpliceSites))
}



#' Test splice sites
#' @importFrom stats model.matrix
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
    if (prime == "start") {
        mySet.all <- which((distSplice.prime != 0) & (spliceStrand != "*") &
            (spliceScore > 0) & (distSplice.prime < 15))
    }else{
        mySet.all <- which((distSplice.prime != 0) & (spliceStrand != "*") &
            (spliceScore > 0) & (abs(distSplice.prime) < 15))
    }
    predictionSplice.prime <- rep(NA, nrow(data))
    if (any(mySet.all)) {
        mySet.training <- 
            which((annotatedSplice.prime | annotatedSplice)[mySet.all])
        myData <- data.frame(spliceScore / (spliceScore + spliceScore.prime),
            spliceScore, distSplice.prime, (spliceStrand.prime == "+"),
            (spliceStrand.prime == "-"), (spliceStrand == "+"))[mySet.all,]
        colnames(myData) <- paste('A',seq_len(ncol(myData)),sep = '.')
        modelmatrix <- 
            model.matrix(~A.1+A.2+A.3+A.4+A.5, data = data.frame(myData))
        predSplice.prime <- NULL
        if (is.null(junctionModel)) { 
            myResults = fitXGBoostModel(labels.train = 
            as.integer(annotatedSplice)[mySet.all][mySet.training], 
            data.train = modelmatrix[mySet.training,],
            data.test = modelmatrix, show.cv = verbose, maxSize.cv = 10000)
            predSplice.prime <- myResults[[2]]
            predictions <- myResults[[1]]
        } else {
            predictions = xgboost:::predict.xgb.Booster(
                junctionModel[[predSplice.primeName]], modelmatrix)
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
#' @importFrom GenomicRanges GRanges
#' @noRd
predictSpliceJunctions <- function(annotatedJunctions, junctionModel=NULL,
    verbose = FALSE){
    spliceVec <- c("Start","End")
    ## if needed this can be a single function
    metadataList <- lapply(spliceVec, function(splice){
        annotatedJunctionsTmp <- 
            GRanges(seqnames = seqnames(annotatedJunctions),
            ranges = IRanges(start = get(tolower(splice))(annotatedJunctions),
            end = get(tolower(splice))(annotatedJunctions)), strand = '*')
        mcols(annotatedJunctionsTmp) <- mcols(annotatedJunctions)
        if (splice == "Start") {
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

#' Fit xgboost model
#' @importFrom xgboost xgboost
#' @importFrom stats fisher.test
#' @noRd
fitXGBoostModel <- function(labels.train, data.train, data.test, 
                            show.cv=TRUE, maxSize.cv=10000){
    if (show.cv) {
        mySample <- sample(seq_along(labels.train),
                           min(floor(length(labels.train)/2),maxSize.cv))
        data.train.cv <- data.train[mySample,]
        labels.train.cv <- labels.train[mySample]
        data.train.cv.test <- data.train[-mySample,]
        labels.train.cv.test <- labels.train[-mySample]

        negative_labels = sum(labels.train == 0)
        positive_labels = sum(labels.train == 1)
        cv.fit <- xgboost(data = data.train.cv, 
            label = labels.train, nthread=2, eta=1, max.depth=5, 
            min_child_weight=5,lambda=0, alpha=10, gamma=0, subsample=0.7, 
            colsample_bytree=0.7, nround= 300, objective = "binary:logistic", 
            eval_metric='error',
            scale_pos_weight=negative_labels/positive_labels, verbose = 0)
        predictions <- predict(cv.fit, data.train.cv.test)
        message('prediction accuracy (CV) (higher for splice 
                donor than splice acceptor)')
        testResults <- fisher.test(table(predictions > 0,labels.train.cv.test))
        show(testResults$estimate)
        show(testResults$p.value)
        show(evalutePerformance(labels.train.cv.test == 1,predictions)$AUC)
    }
    negative_labels = sum(labels.train == 0)
    positive_labels = sum(labels.train == 1)
    cv.fit <- xgboost(data = data.train, 
            label = labels.train, nthread=2, eta=1, max.depth=5, 
            min_child_weight=5,lambda=0, alpha=10, gamma=0, subsample=0.7,
            colsample_bytree=0.7, nround= 300, objective = "binary:logistic", 
            eval_metric='error',
            scale_pos_weight=negative_labels/positive_labels, verbose = 0)
    predictions <- predict(cv.fit, data.test)

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
    setReferenceJunctions <- (apply(spliceSitePredictionList > 0,1,sum) == 4) | 
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
    junctionStrand = as.character(strand(junctions))
    ##note: the output can be visualised (bed/bigbed track) 
    #calculation is based on distance and properties of next junctions
    candidateJunctionsPlus <-
        junctions[which(junctionStrand == '+' | junctionStrand == '*')]
    candidateJunctionsMinus <-
        junctions[which(junctionStrand == '-' | junctionStrand == '*')]
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
    if (verbose) {
        message('reads count for all annotated junctions after 
                correction to reference junction')
        sumByJuncId <- tapply(junctions$score,
                junctions$mergedHighConfJunctionId, sum)
        message(
            sum(sumByJuncId[junctions[names(sumByJuncId)]$annotatedJunction]))
        message(
            sum(sumByJuncId[junctions[names(sumByJuncId)]$annotatedJunction]) / 
                sum(junctions$score))
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
