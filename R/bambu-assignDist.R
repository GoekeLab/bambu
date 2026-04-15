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
    readClassDt[, `:=`(columnIds = as.list(rep(NA, .N)), columnCounts = as.list(rep(NA, .N)))]
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
    incompatibleCounts <- if(sum(incompatibleCountMatrix)==0) NULL else generateIncompatibleCounts(incompatibleCountMatrix, annotations)
    
    distTable <- if(returnDistTable) metadata(metadata(readClassList)$readClassDist)$distTableOld else NULL
    
    readToTranscriptMap <- if(trackReads) generateReadToTranscriptMap(readClassList, 
                                        metadata(readClassList)$readClassDist, 
                                        annotations) else NULL

    quantData <- new("quantData",
        sampleData = data.frame(ColData),
        uniqueCounts = generateUniqueCounts(readClassDt, annotations, nrow(metadata(readClassList)$sampleData)),
        readClassDt = readClassDt,
        incompatibleCountMatrix = incompatibleCountMatrix,
        sampleNames = as.character(metadata(readClassList)$sampleData$sampleName),
        incompatibleCounts = incompatibleCounts,
        nonuniqueCounts = generateNonUniqueCounts(readClassDt, annotations, nrow(metadata(readClassList)$sampleData)),
        distTable = distTable,
        readToTranscriptMap = readToTranscriptMap
    )

    return(quantData)     
}

#' Generate unique counts
#' @noRd
generateUniqueCounts <- function(readClassDt, annotations, nSamples){
    x <- readClassDt %>% filter(!multi_align & !is.na(eqClass.match))
    
    uniqueCounts <- if(nrow(x) == 0) {
        sparseMatrix(i = 1, j = 1, x = 0, dims = c(length(annotations), nSamples))
    } else {
        txids <- mcols(annotations)$txid
        i <- rep(match(x$txid, txids), lengths(x$columnIds))
        j <- unlist(x$columnIds)
        x_vals <- unlist(x$columnCounts)
        sparseMatrix(i = i, j = j, x = x_vals, dims = c(length(annotations), nSamples))
    }
    rownames(uniqueCounts) <- names(annotations)
    return(uniqueCounts)
}


#' Generate incompatible counts
#' @noRd
generateIncompatibleCounts <- function(incompatibleCountMatrix, annotations){
    genes <- levels(factor(unique(mcols(annotations)$GENEID)))
    rownames(incompatibleCountMatrix) <- genes[as.numeric(rownames(incompatibleCountMatrix))]
    geneMat <- sparseMatrix(length(genes), ncol(incompatibleCountMatrix), x = 0)
    rownames(geneMat) <- genes
    geneMat[rownames(incompatibleCountMatrix),] <- incompatibleCountMatrix
    return(geneMat)
}


#' Generate non-unique counts
#' @noRd
generateNonUniqueCounts <- function(readClassDt, annotations, nSamples){
    #fuse multi align RCs by gene
    x <- readClassDt %>% filter(multi_align & !is.na(eqClass.match))
    x <- x %>% distinct(eqClassId, .keep_all = TRUE)
    
    if(nrow(x) == 0) {
        genes <- levels(factor(unique(mcols(annotations)$GENEID)))
        return(sparseMatrix(i = 1, j = 1, x = 0, dims = c(length(genes), nSamples), dimnames = list(genes, NULL)))
    }

    # x has columnIds and columnCounts
    i <- rep(seq_along(x$gene_sid), lengths(x$columnIds))
    j <- unlist(x$columnIds)
    x_vals <- unlist(x$columnCounts)
    nonuniqueCounts <- sparseMatrix(i = i, j = j, x = x_vals, dims = c(nrow(x), nSamples))

    if(nrow(x)>1 & length(unique(x$gene_sid))>1){
        nonuniqueCounts.gene <- sparse.model.matrix(~ factor(x$gene_sid) - 1)
        nonuniqueCounts <- t(nonuniqueCounts.gene) %*% nonuniqueCounts
    } else{
        warning("The factor variable 'gene_sid' has only one level. Adjusting output.")
        nonuniqueCounts.gene <- Matrix(1, nrow = nrow(x), ncol = 1, sparse = TRUE)
        nonuniqueCounts <- t(nonuniqueCounts.gene) %*% nonuniqueCounts
    }
    #covert ids into gene ids
    geneids <- as.numeric(levels(factor(x$gene_sid)))
    geneids <- x$txid[match(geneids, x$gene_sid)]
    geneids <- mcols(annotations)$GENEID[as.numeric(geneids)]
    rownames(nonuniqueCounts) <- geneids
    #create matrix for all annotated genes
    genes <- levels(factor(unique(mcols(annotations)$GENEID)))
    geneMat <- sparseMatrix(length(genes), nSamples, x = 0)
    rownames(geneMat) <- genes
    if(!is.null(rownames(nonuniqueCounts))){
      geneMat[rownames(nonuniqueCounts),] <- nonuniqueCounts
    }
    return(geneMat)
}
