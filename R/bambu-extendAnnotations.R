## Internal functions for bambu =========================
#' Extend annotations
#' @inheritParams bambu
#' @noRd
bambu.extendAnnotations <- function(readClassList, annotations, NDR,
    isoreParameters, stranded, bpParameters, verbose = FALSE) {
    start.ptm_all <- proc.time()
    combinedTxCandidates <- isore.combineTranscriptCandidates(readClassList,
        stranded, ## stranded used for unspliced reduce  
        min.readCount = isoreParameters[["min.readCount"]], 
        min.readFractionByGene = isoreParameters[["min.readFractionByGene"]],
        min.txScore.multiExon = isoreParameters[["min.txScore.multiExon"]],
        min.txScore.singleExon = isoreParameters[["min.txScore.singleExon"]],
        bpParameters,
        verbose)
    end.ptm_all <- proc.time()
    if (verbose) message("combining transcripts in ",
        round((end.ptm_all - start.ptm_all)[3] / 60, 1)," mins.")
    start.ptm_all <- proc.time()
    annotations <- isore.extendAnnotations(
        combinedTranscripts = combinedTxCandidates,
        annotationGrangesList = annotations,
        remove.subsetTx = isoreParameters[["remove.subsetTx"]],
        min.sampleNumber = isoreParameters[["min.sampleNumber"]],
        NDR = NDR,
        min.exonDistance = isoreParameters[["min.exonDistance"]],
        min.exonOverlap = isoreParameters[["min.exonOverlap"]],
        min.primarySecondaryDist = 
        isoreParameters[['min.primarySecondaryDist']], 
        min.primarySecondaryDistStartEnd = 
        isoreParameters[['min.primarySecondaryDistStartEnd1']],
        prefix = isoreParameters[["prefix"]],
        txScoreBaseline = isoreParameters[["txScoreBaseline"]],
        verbose = verbose)
    end.ptm_all <- proc.time()
    if (verbose) message("extend annotations in ",
                         round((end.ptm_all - start.ptm_all)[3] / 60, 1)," mins.")
    return(annotations)
}
