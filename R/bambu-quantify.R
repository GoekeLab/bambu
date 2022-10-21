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
    compatibleCounts <- bambu.quantDT(readClassDt, emParameters = emParameters,
                               ncore = ncore, verbose = verbose)
    incompatibleCounts <- incompatibleCounts[data.table(GENEID = unique(mcols(annotations)$GENEID)), on = "GENEID"]
    incompatibleCounts[is.na(counts), counts := 0]
    setnames(incompatibleCounts, "counts", colnames(readClass))
    counts <- compatibleCounts[match(mcols(annotations)$txid), txid)]
    colNameRC <- colnames(readClass)
    colDataRC <- colData(readClass)
    seOutput <- SummarizedExperiment(
        assays = SimpleList(counts = matrix(counts$counts, ncol = 1,
        dimnames = list(NULL, colNameRC)), CPM = matrix(counts$CPM,
        ncol =  1, dimnames = list(NULL, colNameRC)),
        fullLengthCounts = matrix(counts$fullLengthCounts, ncol = 1,
            dimnames = list(NULL, colNameRC)),
        uniqueCounts = matrix(counts$uniqueCounts, 
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
    rcPreOut <- addAval(readClassDt, emParameters, verbose)
    readClassDt <- rcPreOut[[1]]
    outIni <- initialiseOutput(readClassDt)
    readClassDt <- filterTxRc(readClassDt) 
    readClassDt <- assignGroups(readClassDt)
    inputRcDt <- getInputList(readClassDt)
    readClassDt <- split(readClassDt, by = "gene_grp_id")
    start.ptm <- proc.time()
    outEst <- abundance_quantification(inputRcDt, readClassDt,
                                       ncore = ncore,
                                       maxiter = emParameters[["maxiter"]],
                                       conv = emParameters[["conv"]], minvalue = emParameters[["minvalue"]])
    end.ptm <- proc.time()
    if (verbose) message("Finished EM estimation in ",
                         round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    outEst <- modifyQuantOut(outEst,outIni)
    theta_est <- formatOutput(rbind(rcPreOut[[2]],outEst))
    theta_est <- removeDuplicates(theta_est)
    return(theta_est)
}


