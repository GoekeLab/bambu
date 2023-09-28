#' Perform quantification
#' @inheritParams bambu
#' @import data.table
#' @noRd
bambu.quantify <- function(i, readClassDist, readClassDt, countMatrix, annotations, emParameters, 
                           trackReads = FALSE, returnDistTable = FALSE,
                           verbose = FALSE, isoreParameters = setIsoreParameters(NULL)) {
    metadata(readClassDist)$distTable$readCount <- countMatrix[,i] 
    metadata(readClassDist)$distTable = metadata(readClassDist)$distTable[metadata(readClassDist)$distTable$readCount != 0,]
    #distTable$readCount <- assays(seReadClass)$counts
    

    #calculate equivilent class counts
    start.ptm <- proc.time()
    y = metadata(readClassDist)$distTable %>% group_by(eqClassById) %>%
    mutate(anyEqual = any(equal)) %>%
    select(eqClassById, firstExonWidth,totalWidth, readCount,GENEID,anyEqual) %>% #eqClassByIdTemp,
    distinct() %>%
    mutate(nobs = sum(readCount),
            rcWidth = ifelse(anyEqual, max(totalWidth), 
                            max(firstExonWidth))) %>%
    select(eqClassById,GENEID,nobs,rcWidth) %>% 
    ungroup()  %>%
    distinct()
    readClassDt$nobs = y$nobs[match(readClassDt$eqClassById,y$eqClassById)]
    readClassDt$nobs[is.na(readClassDt$nobs)] = 0
    end.ptm <- proc.time()
    if (verbose) message("Finished eqv.class counts in ",
                         round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    ####

    incompatibleCounts <- processIncompatibleCounts(readClassDist)
    compatibleCounts <- bambu.quantDT(readClassDt, emParameters = emParameters,verbose = verbose)
    incompatibleCounts <- incompatibleCounts[data.table(GENEID = unique(mcols(annotations)$GENEID)), on = "GENEID"]
    incompatibleCounts[is.na(counts), counts := 0]
    compatibleCounts <- calculateCPM(compatibleCounts, incompatibleCounts)
    setnames(incompatibleCounts, "counts", colnames(countMatrix)[i])
    counts <- compatibleCounts[match(mcols(annotations)$txid, txid)]
    colNameRC <- colnames(countMatrix)[i]
    colDataRC <- NULL # TODO carry over coldata??!!
    sig.digit <- emParameters[["sig.digit"]]
    seOutput <- SummarizedExperiment(
        assays = SimpleList(counts = matrix(round(counts$counts,sig.digit), ncol = 1,
        dimnames = list(NULL, colNameRC)), CPM = matrix(round(counts$CPM,sig.digit),
        ncol =  1, dimnames = list(NULL, colNameRC)),
        fullLengthCounts = matrix(round(counts$fullLengthCounts,sig.digit), ncol = 1,
            dimnames = list(NULL, colNameRC)),
        uniqueCounts = matrix(counts$uniqueCounts, 
            ncol = 1, dimnames = list(NULL, colNameRC))))
    metadata(seOutput)$incompatibleCounts = incompatibleCounts
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
                                       maxiter = emParameters[["maxiter"]],
                                       conv = emParameters[["conv"]], minvalue = emParameters[["minvalue"]])
    end.ptm <- proc.time()
    if (verbose) message("Finished EM estimation in ",
                         round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    outEst <- modifyQuantOut(outEst,outIni)
    theta_est <- rbind(rcPreOut[[2]],outEst)
    theta_est <- removeDuplicates(theta_est)
    return(theta_est)
}


