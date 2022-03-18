#' Perform quantification
#' @inheritParams bambu
#' @import data.table
#' @noRd
bambu.quantify <- function(readClass, annotations, emParameters,ncore = 1,
    verbose = FALSE, isoreParameters = setIsoreParameters(NULL)) {
    min.exonDistance = isoreParameters[["min.exonDistance"]]
    min.primarySecondaryDist =
        isoreParameters[['min.primarySecondaryDist']] 
    min.primarySecondaryDistStartEnd =
        isoreParameters[['min.primarySecondaryDistStartEnd2']]
    if (is.character(readClass)) readClass <- readRDS(file = readClass)
    readClass <- isore.estimateDistanceToAnnotations(readClass, annotations,
        min.exonDistance = min.exonDistance,
        min.primarySecondaryDist = min.primarySecondaryDist,
        min.primarySecondaryDistStartEnd = min.primarySecondaryDistStartEnd,
        verbose = verbose)
    readClassDt <- genEquiRCs(readClass, annotations)
    tx_len <- rbind(data.table(tx_id = names(annotations),
        tx_len = sum(width(annotations))),
                    data.table(tx_id = paste0(names(annotations),"Start"),
                               tx_len = sum(width(annotations))))
    readClassDt <- unique(tx_len[readClassDt, on = "tx_id"])
    countsOut <- bambu.quantDT(readClassDt, emParameters = emParameters,
        ncore = ncore, verbose = verbose)
    counts <- countsOut[[1]]
    counts <- counts[match(names(annotations), tx_name)]
    colNameRC <- colnames(readClass)
    colDataRC <- cbind(colData(readClass), d_rate = countsOut[[2]],
        nGeneFordRate = countsOut[[3]])
    seOutput <- SummarizedExperiment(
        assays = SimpleList(counts = matrix(counts$counts, ncol = 1,
        dimnames = list(NULL, colNameRC)), CPM = matrix(counts$CPM,
        ncol =  1, dimnames = list(NULL, colNameRC)),
        fullLengthCounts = matrix(counts$FullLengthCounts, ncol = 1,
            dimnames = list(NULL, colNameRC)),
        partialLengthCounts = matrix(counts$PartialLengthCounts, 
            ncol = 1, dimnames = list(NULL, colNameRC)),
        uniqueCounts = matrix(counts$UniqueCounts, 
            ncol = 1, dimnames = list(NULL, colNameRC)),
        theta = matrix(counts$theta, 
            ncol = 1, dimnames = list(NULL, colNameRC))), colData = colDataRC)
    metadata(seOutput)$distTable = metadata(readClass)$distTable
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

generateReadToTranscriptMap <- function(readClass, annotations){
    distTable = metadata(readClass[[2]])$distTable
    readClass = readClass[[1]]
    if (is.character(readClass)) readClass <- readRDS(file = readClass)
    if(!is.null(metadata(readClass)$readNames)) { 
        read_id = metadata(readClass)$readNames}
    else { read_id = metadata(readClass)$readId}

    print(distTable)

    #unpack and reverse the read class to read id reletionship
    readOrder = order(unlist(rowData(readClass)$readIds))
    lens = lengths(rowData(readClass)$readIds)
    rcIndex = seq_along(readClass)
    readToRC = rep(rcIndex, lens)[readOrder]

    #get annotation indexs
    distTable$annotationTxId = match(distTable$annotationTxId, names(annotations))

    #match read classes with transcripts
    readClass_id = rownames(readClass)[readToRC]
    equal_matches = distTable %>% 
        filter(equal) %>% 
        group_by(readClassId) %>% summarise(annotationTxIds = list(annotationTxId))
    equal_matches = equal_matches$annotationTxIds[match(readClass_id, equal_matches$readClassId)]
    compatible_matches = distTable %>% 
        filter(!equal & compatible) %>% 
        group_by(readClassId) %>% summarise(annotationTxIds = list(annotationTxId))
    compatible_matches = compatible_matches$annotationTxIds[match(readClass_id, compatible_matches$readClassId)]
    readToTranscriptMap = data.table(read_id=read_id, equal_matches = equal_matches, compatible_matches = compatible_matches)
    #readToTranscriptMap = readToTranscriptMap[!is.na(transcript_id),]

    return(readToTranscriptMap)
}
