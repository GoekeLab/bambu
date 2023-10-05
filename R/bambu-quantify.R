#' Perform quantification
#' @inheritParams bambu
#' @import data.table
#' @noRd
bambu.quantify <- function(distTable, readClassDt, countMatrix, txid.index, GENEIDs, emParameters, 
                           trackReads = FALSE, returnDistTable = FALSE,
                           verbose = FALSE, isoreParameters = setIsoreParameters(NULL)) {
    distTable$readCount <- countMatrix
    distTable = distTable[distTable$readCount != 0,]
    readClassDt$nobs = calculateEqClassCounts(distTable, readClassDt)
    incompatibleCounts <- processIncompatibleCounts(distTable)
    compatibleCounts <- bambu.quantDT(readClassDt, emParameters = emParameters,verbose = verbose)
    incompatibleCounts <- incompatibleCounts[data.table(GENEID.i = GENEIDs), on = "GENEID.i"]
    incompatibleCounts[is.na(counts), counts := 0]
    compatibleCounts <- calculateCPM(compatibleCounts, incompatibleCounts)
    counts <- compatibleCounts[match(txid.index, txid)]
    sig.digit <- emParameters[["sig.digit"]]
    seOutput <- list(incompatibleCounts = as(incompatibleCounts$counts, "sparseVector"),
                    counts = as(round(counts$counts,sig.digit), "sparseVector"),
                    CPM = as(round(counts$CPM,sig.digit), "sparseVector"),
                    fullLengthCounts = as(round(counts$fullLengthCounts,sig.digit), "sparseVector"),
                    uniqueCounts = as(round(counts$uniqueCounts,sig.digit), "sparseVector"))                   
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


