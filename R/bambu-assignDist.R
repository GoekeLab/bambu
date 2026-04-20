#' Create equivilence classes and assign to transcripts
#' @inheritParams bambu
#' @import data.table
#' @noRd
assignReadClasstoTranscripts <- function(readClassList, annotations, isoreParameters, 
                                        verbose, sampleMetadata, demultiplexed,
                                        returnDistTable = FALSE, trackReads = TRUE) {
    if (is.character(readClassList)) readClassList <- readRDS(file = readClassList)
    metadata(readClassList)$readClassDist <- calculateDistTable(readClassList, annotations, isoreParameters, verbose, returnDistTable)
    readClassList <- splitReadClassFiles(readClassList)
    readClassDt <- genEquiRCs(metadata(readClassList)$readClassDist, annotations, verbose) 
    readClassDt$eqClass.match = match(readClassDt$eqClassById,metadata(readClassList)$eqClassById)

    # Add columnIds and columnCounts columns
    readClassDt[, `:=`(columnIds = vector("list", .N), columnCounts = vector("list", .N))]
    observedIdx <- which(!is.na(readClassDt$eqClass.match))
    readClassDt$columnIds[observedIdx] <- metadata(readClassList)$columnIds[readClassDt$eqClass.match[observedIdx]]
    readClassDt$columnCounts[observedIdx] <- metadata(readClassList)$columnCounts[readClassDt$eqClass.match[observedIdx]]

    readClassDt <- simplifyNames(readClassDt)
    readClassDt <- readClassDt %>% group_by(eqClassId, gene_sid) %>% 
        mutate(multi_align = length(unique(txid))>1) %>% 
        ungroup() %>% 
        mutate(aval = 1) %>%
        data.table()
    #return non-em counts
    ColData <- generateColData(readClassList, sampleMetadata, demultiplexed)
    
    incompatibleCountMatrix <- metadata(readClassList)$incompatibleCountMatrix
    incompatibleCounts <- generateIncompatibleCounts(incompatibleCountMatrix, annotations)

    distTable <- if(returnDistTable) {
        metadata(metadata(readClassList)$readClassDist)$distTableOld
    } else{
        NULL
    }

    readToTranscriptMap <- if(trackReads) {
        generateReadToTranscriptMap(readClassList, metadata(readClassList)$readClassDist,annotations)
    } else{
        NULL
    }

    quantData <- constructQuantData(
        sampleData          = data.frame(ColData),
        readClassDt         = readClassDt,
        incompatibleCounts  = incompatibleCounts,
        distTable           = distTable,
        readToTranscriptMap = readToTranscriptMap
    )

    return(quantData)     
}

#' Generate incompatible counts
#' @noRd
generateIncompatibleCounts <- function(incompatibleCountMatrix, annotations){
    genes <- levels(factor(unique(mcols(annotations)$GENEID)))
    rownames(incompatibleCountMatrix) <- genes[as.numeric(rownames(incompatibleCountMatrix))]
    geneMat <- sparseMatrix(length(genes), ncol(incompatibleCountMatrix), x = 0)
    rownames(geneMat) <- genes
    colnames(geneMat) <- colnames(incompatibleCountMatrix)
    geneMat[match(rownames(incompatibleCountMatrix), rownames(geneMat)), ] <- incompatibleCountMatrix
    return(geneMat[unique(mcols(annotations)$GENEID), , drop = FALSE])
}
