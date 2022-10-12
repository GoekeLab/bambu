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
    readClassMod <- modifyIncompatibleAssignment(readClassDist)
    annotationsUpdated <- updateAnnotations(readClassMod, annotations, verbose)
    readClassDt <- genEquiRCs(readClassMod, annotationsUpdated, verbose)
    tx_len <- rbind(data.table(tx_id = names(annotationsUpdated),
                               tx_len = sum(width(annotationsUpdated))),
                    data.table(tx_id = paste0(names(annotationsUpdated),"Start"),
                               tx_len = sum(width(annotationsUpdated))))
    readClassDt <- unique(tx_len[readClassDt, on = "tx_id"])
    countsOut <- bambu.quantDT(readClassDt, emParameters = emParameters,
                               ncore = ncore, verbose = verbose)
    counts <- countsOut[[1]]
    annoDt <- data.table(tx_name = names(annotationsUpdated))
    if(any(!(counts$tx_name %in% annoDt$tx_name))) 
        warning("For developer: transcript names in counts table not
    found in updated annotations, please check!")
    incompatibleCounts <- data.table(GENEID = unique(mcols(annotations)$GENEID))
    counts_incompatible <- counts[grepl("_unidentified",tx_name)]
    counts_incompatible[, GENEID := gsub("_unidentified","",tx_name)]
    
    incompatibleCounts <- counts_incompatible[,.(GENEID, counts)][incompatibleCounts, on = "GENEID"]
    incompatibleCounts[is.na(counts), counts := 0]
    setnames(incompatibleCounts, "counts", colnames(readClass))
    counts <- counts[match(names(annotations),tx_name)]
    colNameRC <- colnames(readClass)
    colDataRC <- cbind(colData(readClass), d_rate = countsOut[[2]],
                       nGeneFordRate = countsOut[[3]])
    seOutput <- SummarizedExperiment(
        assays = SimpleList(counts = matrix(counts$counts, ncol = 1,
        dimnames = list(NULL, colNameRC)), CPM = matrix(counts$CPM,
        ncol =  1, dimnames = list(NULL, colNameRC)),
        fullLengthCounts = matrix(counts$FullLengthCounts, ncol = 1,
            dimnames = list(NULL, colNameRC)),
        uniqueCounts = matrix(counts$UniqueCounts, 
            ncol = 1, dimnames = list(NULL, colNameRC))), colData = colDataRC)
    metadata(seOutput)$incompatibleCounts = incompatibleCounts
    if (returnDistTable) metadata(seOutput)$distTable = metadata(readClassMod)$distTable
    if (trackReads) metadata(seOutput)$readToTranscriptMap = 
        generateReadToTranscriptMap(readClass, metadata(readClassMod)$distTable, 
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
    } else if (any(!(c("gene_id", "tx_id", "read_class_id","nobs") %in% 
                     colnames(readClassDt)))) {
        stop("Columns gene_id, tx_id, read_class_id, nobs,
            are missing from object.")
    }
    ## ----step2: match to simple numbers to increase claculation efficiency
    geneVec <- unique(readClassDt$gene_id)
    ori_txvec <- unique(gsub("Start","",readClassDt$tx_id))
    txVec <- unique(readClassDt$tx_id)
    readclassVec <- unique(readClassDt$read_class_id)
    readClassDt <- 
        simplifyNames(readClassDt,txVec, geneVec,ori_txvec, readclassVec)
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
    theta_est <- formatOutput(rbind(outList,outListEst),ori_txvec,geneVec)
    theta_est <- removeDuplicates(theta_est)
    return(list(theta_est, d_rateOut[1], d_rateOut[2]))
}

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

