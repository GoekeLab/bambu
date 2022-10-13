#' Perform quantification
#' @inheritParams bambu
#' @import data.table
#' @noRd
bambu.quantify <- function(readClass, annotations, emParameters, 
                           trackReads = FALSE, returnDistTable = FALSE, ncore = 1,
                           verbose = FALSE, isoreParameters = setIsoreParameters(NULL)) {
    min.exonDistance = isoreParameters[["min.exonDistance"]]
    min.primarySecondaryDist =
        isoreParameters[['min.primarySecondaryDist']] 
    min.primarySecondaryDistStartEnd =
        isoreParameters[['min.primarySecondaryDistStartEnd2']]
    if (is.character(readClass)) readClass <- readRDS(file = readClass)
    readClassDist <- isore.estimateDistanceToAnnotations(readClass, annotations,
                                                         min.exonDistance = min.exonDistance,
                                                         min.primarySecondaryDist = min.primarySecondaryDist,
                                                         min.primarySecondaryDistStartEnd = min.primarySecondaryDistStartEnd,
                                                         verbose = verbose)
    metadata(readClassDist)$distTable <- modifyIncompatibleAssignment(metadata(readClassDist)$distTable)
    incompatibleCounts <- processIncompatibleCounts(readClassDist)
    readClassDt <- genEquiRCs(readClassDist, annotations, verbose) 
    tx_len <- rbind(data.table(txid = mcols(annotations)$txid,
                               txlen = sum(width(annotations))))
    readClassDt <- tx_len[readClassDt, on = "txid"] %>% distinct()
    compatibleCounts <- bambu.quantDT(readClassDt, emParameters = emParameters,
                               ncore = ncore, verbose = verbose)
    compatibleCounts[, tx_name := names(annotations)[as.numeric(txid)]]
    incompatibleCounts <- incompatibleCounts[data.table(GENEID = unique(mcols(annotations)$GENEID)), on = "GENEID"]
    incompatibleCounts[is.na(counts), counts := 0]
    setnames(incompatibleCounts, "counts", colnames(readClass))
    counts <- compatibleCounts[match(names(annotations),tx_name)]
    colNameRC <- colnames(readClass)
    colDataRC <- colData(readClass)
    seOutput <- SummarizedExperiment(
        assays = SimpleList(counts = matrix(counts$counts, ncol = 1,
        dimnames = list(NULL, colNameRC)), CPM = matrix(counts$CPM,
        ncol =  1, dimnames = list(NULL, colNameRC)),
        fullLengthCounts = matrix(counts$FullLengthCounts, ncol = 1,
            dimnames = list(NULL, colNameRC)),
        uniqueCounts = matrix(counts$UniqueCounts, 
            ncol = 1, dimnames = list(NULL, colNameRC))), colData = colDataRC)
    metadata(seOutput)$incompatibleCounts = incompatibleCounts
    if (returnDistTable) metadata(seOutput)$distTable = metadata(readClassDist)$distTable
    if (trackReads) metadata(seOutput)$readToTranscriptMap = 
        generateReadToTranscriptMap(readClass, metadata(readClassDist)$distTable, 
                                    annotations)
    return(seOutput)
}





#' Process data.table object
#' @param readClassDt A data.table object
#' @inheritParams bambu
#' @noRd
bambu.quantDT <- function(readClassDt = readClassDt, 
                          emParameters = list(degradationBias = TRUE, maxiter = 10000, conv = 10^(-2),
                                              minvalue = 10^(-8)), ncore = 1, verbose = FALSE) {
    if (is.null(readClassDt)) {
        stop("Input object is missing.")
    } else if (any(!(c("GENEID", "txid", "eqClassId","nobs") %in% 
                     colnames(readClassDt)))) {
        stop("Columns GENEID, txid, eqClassId, nobs,
            are missing from object.")
    }
    ## ----step2: match to simple numbers to increase claculation efficiency
    geneVec <- unique(readClassDt$GENEID)
    txVec <- unique(readClassDt$txid)
    eqClassMap <- readClassDt[,.(eqClassById, eqClassId)] %>% distinct()
    readClassDt <- 
        simplifyNames(readClassDt,geneVec)
    d_mode <- emParameters[["degradationBias"]]
    start.ptm <- proc.time()
    if (d_mode) {
        d_rateOut <- calculateDegradationRate(readClassDt)
    }else{
        d_rateOut <- rep(NA,2)
    }
    end.ptm <- proc.time()
    if (verbose) message("Finished estimate degradation bias in ",
                         round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    readClassDt <- modifyAvaluewithDegradation_rate(readClassDt, 
                                                    d_rateOut[1], d_mode = d_mode)
    removeList <- removeUnObservedGenes(readClassDt)
    readClassDt <- removeList[[1]] # keep only observed genes for estimation
    outList <- removeList[[2]] #for unobserved genes, set estimates to 0 
    start.ptm <- proc.time()
    outListEst <- abundance_quantification(readClassDt,
                                           ncore = ncore,
                                           maxiter = emParameters[["maxiter"]],
                                           conv = emParameters[["conv"]], minvalue = emParameters[["minvalue"]])
    end.ptm <- proc.time()
    if (verbose) message("Finished EM estimation in ",
                         round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    theta_est <- formatOutput(rbind(outList,outListEst),txVec,geneVec)
    theta_est <- removeDuplicates(theta_est)
    return(theta_est)
}


#' Generate read to transcript mapping
#' @noRd
generateReadToTranscriptMap <- function(readClass, distTable, annotations){
    if(!is.null(metadata(readClass)$readNames)) { 
        read_id = metadata(readClass)$readNames}
    else { read_id = metadata(readClass)$readId}
    #unpack and reverse the read class to read id relationship
    readOrder = order(unlist(rowData(readClass)$readIds))
    lens = lengths(rowData(readClass)$readIds)
    rcIndex = seq_along(readClass)
    readToRC = rep(rcIndex, lens)[readOrder]
    read_id = read_id[(match(unlist(rowData(readClass)$readIds)[readOrder], metadata(readClass)$readId))]
    #get annotation indexs
    distTable$annotationTxId = match(distTable$annotationTxId, names(annotations))
    #match read classes with transcripts
    readClass_id = rownames(readClass)[readToRC]
    equalMatches = distTable %>% 
        filter(equal) %>% 
        group_by(readClassId) %>% summarise(annotationTxIds = list(annotationTxId))
    equalMatches = equalMatches$annotationTxIds[match(readClass_id, equalMatches$readClassId)]
    compatibleMatches = distTable %>% 
        filter(!equal & compatible) %>% 
        group_by(readClassId) %>% summarise(annotationTxIds = list(annotationTxId))
    compatibleMatches = compatibleMatches$annotationTxIds[match(readClass_id, compatibleMatches$readClassId)]
    readToTranscriptMap = tibble(readId=read_id, equalMatches = equalMatches, compatibleMatches = compatibleMatches)
    
    return(readToTranscriptMap)
}


#' Process incompatible counts
#' @noRd
processIncompatibleCounts <- function(readClassDist){
    distTable <- data.table(as.data.frame(metadata(readClassDist)$distTable))[, 
        .(readClassId, annotationTxId, readCount, GENEID, dist,equal)]
    distTableIncompatible <- distTable[grep("unidentified", annotationTxId)]
    # filter out multiple geneIDs mapped to the same readClass using rowData(se)
    geneRCMap <- as.data.table(as.data.frame(rowData(readClassDist)),
                                    keep.rownames = TRUE)
    setnames(geneRCMap, old = c("rn", "geneId"),
             new = c("readClassId", "GENEID"))
    distTable <- distTable[geneRCMap[ readClassId %in% 
                                               unique(distTableIncompatible$readClassId), .(readClassId, GENEID)],
                           on = c("readClassId", "GENEID")]
    counts <- distTable[,.(GENEID, readCount)]
    setnames(counts, "readCount", "counts")
    return(counts)
    
}