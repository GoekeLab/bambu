## Internal functions for bambu =========================
#' Extend annotations
#' @inheritParams bambu
#' @noRd
bambu.extendAnnotations <- function(readClassList, annotations,
    isoreParameters, stranded, bpParameters, verbose = FALSE) {
    start.ptm <- proc.time()
    combinedTxCandidates <- isore.combineTranscriptCandidates(readClassList,
        stranded, ## stranded used for unspliced reduce  
        min.readCount = isoreParameters[["min.readCount"]], 
        min.readFractionByGene = isoreParameters[["min.readFractionByGene"]],
        min.geneFDR = isoreParameters[["min.geneFDR"]],
        min.txFDR = isoreParameters[["min.txFDR"]],
        bpParameters,
        verbose)
    end.ptm <- proc.time()
    if (verbose) message("combining transcripts in ",
        round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
    start.ptm <- proc.time()
    annotations <- isore.extendAnnotations(
        combinedTranscripts = combinedTxCandidates,
        annotationGrangesList = annotations,
        remove.subsetTx = isoreParameters[["remove.subsetTx"]],
        min.sampleNumber = isoreParameters[["min.sampleNumber"]],
        min.exonDistance = isoreParameters[["min.exonDistance"]],
        min.exonOverlap = isoreParameters[["min.exonOverlap"]],
        min.primarySecondaryDist = 
        isoreParameters[['min.primarySecondaryDist']], 
        min.primarySecondaryDistStartEnd = 
        isoreParameters[['min.primarySecondaryDistStartEnd1']],
        prefix = isoreParameters[["prefix"]],
        verbose = verbose)
    end.ptm <- proc.time()
    if (verbose) message("extend annotations in ",
                         round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
    return(annotations)
}
