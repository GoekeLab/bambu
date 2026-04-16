#' Perform quantification
#' @inheritParams bambu
#' @import data.table
#' @noRd
bambu.quantify <- function(readClassDt, columnIdx, incompatibleCounts, nonuniqueCounts, txid.index, GENEIDs, emParameters,
                           trackReads = FALSE, returnDistTable = FALSE,
                           verbose = FALSE, isoreParameters = setIsoreParameters(NULL)) {
    start.ptm <- proc.time()

    # Calculate nobs for sample(s) columnIdx
    # Use data.table syntax for in-place creation of nobs column for the scope of this function
    readClassDt[, nobs := {
        ids <- columnIds[[1]]
        if (is.null(ids) || length(ids) == 0 || is.na(ids[1])) {
            0L
        } else if (length(columnIdx) == 1) {
            match_idx <- match(columnIdx, ids)
            if (is.na(match_idx)) 0L else columnCounts[[1]][match_idx]
        } else {
            sum(columnCounts[[1]][ids %in% columnIdx])
        }
    }, by = eqClassId]

    compatibleCounts <- bambu.quantDT(readClassDt, emParameters = emParameters,verbose = verbose)
    incompatibleCounts <- incompatibleCounts[data.table(GENEID.i = GENEIDs), on = "GENEID.i"]
    incompatibleCounts[is.na(counts), counts := 0]
    compatibleCounts <- calculateCPM(compatibleCounts, incompatibleCounts)
    counts <- compatibleCounts[match(txid.index, txid)]
    sig.digit <- emParameters[["sig.digit"]]
    seOutput <- list(incompatibleCounts = as(incompatibleCounts$counts, "sparseVector"),
                    nonuniqueCounts = as(nonuniqueCounts, "sparseVector"),
                    counts = as(round(counts$counts,sig.digit), "sparseVector"),
                    CPM = as(round(counts$CPM,sig.digit), "sparseVector"),
                    fullLengthCounts = as(round(counts$fullLengthCounts,sig.digit), "sparseVector"),
                    uniqueCounts = as(round(counts$uniqueCounts,sig.digit), "sparseVector"))              
    end.ptm <- proc.time()
    # if (verbose) message("bambu.quantify ", round((end.ptm - start.ptm)[3] / 60, 3), " mins.")
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
    # if (verbose) message("Finished EM estimation in ",
    #                     round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    outEst <- modifyQuantOut(outEst,outIni)
    theta_est <- rbind(rcPreOut[[2]],outEst)
    theta_est <- removeDuplicates(theta_est)
    end.ptm <- proc.time()
    return(theta_est)
}


