## Internal functions for bambu =========================
#' Extend annotations
#' @inheritParams bambu
#' @noRd
bambu.extendAnnotations <- function(readClassList, annotations, NDR,
    discoveryParameters, stranded, bpParameters, fusionMode = FALSE, verbose = FALSE) {
    start.ptm_all <- proc.time()
    combinedTxCandidates <- isore.combineTranscriptCandidates(readClassList,
        stranded, ## stranded used for unspliced reduce  
        min.readCount = discoveryParameters[["min.readCount"]], 
        min.readFractionByGene = discoveryParameters[["min.readFractionByGene"]],
        min.txScore.multiExon = discoveryParameters[["min.txScore.multiExon"]],
        min.txScore.singleExon = discoveryParameters[["min.txScore.singleExon"]],
        bpParameters,
        verbose)
    end.ptm_all <- proc.time()
    if (verbose) message("combining transcripts in ",
        round((end.ptm_all - start.ptm_all)[3] / 60, 1)," mins.")
    start.ptm_all <- proc.time()
    annotations <- isore.extendAnnotations(
        combinedTranscripts = combinedTxCandidates,
        annotationGrangesList = annotations,
        remove.subsetTx = discoveryParameters[["remove.subsetTx"]],
        min.sampleNumber = discoveryParameters[["min.sampleNumber"]],
        NDR = NDR,
        min.exonDistance = discoveryParameters[["min.exonDistance"]],
        min.exonOverlap = discoveryParameters[["min.exonOverlap"]],
        min.primarySecondaryDist = 
        discoveryParameters[['min.primarySecondaryDist']], 
        min.primarySecondaryDistStartEnd = 
        discoveryParameters[['min.primarySecondaryDistStartEnd1']],
        min.readFractionByEqClass =  
        discoveryParameters[['min.readFractionByEqClass']],
        fusionMode = fusionMode,
        prefix = discoveryParameters[["prefix"]],
        baselineFDR = discoveryParameters[["baselineFDR"]],
        defaultModels = discoveryParameters[["defaultModels"]],
        verbose = verbose)
    end.ptm_all <- proc.time()
    if (verbose) message("extend annotations in ",
                         round((end.ptm_all - start.ptm_all)[3] / 60, 1)," mins.")
    return(annotations)
}
