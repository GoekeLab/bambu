## Internal functions for bambu =========================
#' Extend annotations
#' @inheritParams bambu
#' @noRd
bambu.extendAnnotations <- function(readClassList, annotations, NDR,
    isoreParameters, rcAssignmentParameters, stranded, bpParameters, fusionMode = FALSE, verbose = FALSE) {
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
        min.exonDistance = rcAssignmentParameters[["min.exonDistance"]],
        min.exonOverlap = isoreParameters[["min.exonOverlap"]],
        min.primarySecondaryDist =
        rcAssignmentParameters[['min.primarySecondaryDist']],
        min.primarySecondaryDistStartEnd =
        rcAssignmentParameters[['min.primarySecondaryDistStartEnd1']],
        min.readFractionByEqClass =  
        isoreParameters[['min.readFractionByEqClass']],
        fusionMode = fusionMode,
        prefix = isoreParameters[["prefix"]],
        baselineFDR = isoreParameters[["baselineFDR"]],
        defaultModels = isoreParameters[["defaultModels"]],
        verbose = verbose)
    end.ptm_all <- proc.time()
    if (verbose) message("extend annotations in ",
                         round((end.ptm_all - start.ptm_all)[3] / 60, 1)," mins.")
    return(annotations)
}
