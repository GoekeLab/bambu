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
    incompatibleCounts <- generateIncompatibleCounts(incompatibleCountMatrix, annotations)
    
    distTable <- if(returnDistTable) metadata(metadata(readClassList)$readClassDist)$distTableOld else NULL
    
    readToTranscriptMap <- if(trackReads) generateReadToTranscriptMap(readClassList, 
                                        metadata(readClassList)$readClassDist, 
                                        annotations) else NULL

    quantData <- new("quantData",
        sampleData = data.frame(ColData),
        readClassDt = readClassDt,
        incompatibleCounts = incompatibleCounts,
        distTable = distTable,
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
    return(geneMat)
}

#' Generate unique counts
#' @noRd
generateUniqueCountsFromQuantData <- function(quantData, annotations){
    uniqueCountsList <- lapply(quantData, function(x) {
        readClassDt <- x$readClassDt
        x_filtered <- readClassDt %>% filter(!multi_align & !is.na(eqClass.match))
        
        uniqueCounts <- if(nrow(x_filtered) == 0) {
            sparseMatrix(i = 1, j = 1, x = 0, dims = c(length(annotations), nrow(x$sampleData)))
        } else {
            txids <- mcols(annotations)$txid
            i <- rep(match(x_filtered$txid, txids), lengths(x_filtered$columnIds))
            j <- unlist(x_filtered$columnIds)
            x_vals <- unlist(x_filtered$columnCounts)
            sparseMatrix(i = i, j = j, x = x_vals, dims = c(length(annotations), nrow(x$sampleData)))
        }
        rownames(uniqueCounts) <- names(annotations)
        colnames(uniqueCounts) <- rownames(x$sampleData)
        return(uniqueCounts)
    })
    return(do.call(cbind, uniqueCountsList))
}

#' Generate non-unique counts
#' @noRd
generateNonUniqueCountsFromQuantData <- function(quantData, annotations){
    nonuniqueCountsList <- lapply(quantData, function(x_obj) {
        readClassDt <- x_obj$readClassDt
        #fuse multi align RCs by gene
        x <- readClassDt %>% filter(multi_align & !is.na(eqClass.match))
        x <- x %>% distinct(eqClassId, .keep_all = TRUE)
        
        if(nrow(x) == 0) {
            genes <- levels(factor(unique(mcols(annotations)$GENEID)))
            geneMat <- sparseMatrix(i = 1, j = 1, x = 0, dims = c(length(genes), nrow(x_obj$sampleData)), dimnames = list(genes, NULL))
            colnames(geneMat) <- rownames(x_obj$sampleData)
            return(geneMat)
        }

        # x has columnIds and columnCounts
        i <- rep(seq_along(x$gene_sid), lengths(x$columnIds))
        j <- unlist(x$columnIds)
        x_vals <- unlist(x$columnCounts)
        nonuniqueCounts <- sparseMatrix(i = i, j = j, x = x_vals, dims = c(nrow(x), nrow(x_obj$sampleData)))

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
        geneMat <- sparseMatrix(length(genes), nrow(x_obj$sampleData), x = 0)
        rownames(geneMat) <- genes
        if(!is.null(rownames(nonuniqueCounts))){
          geneMat[match(rownames(nonuniqueCounts), rownames(geneMat)), ] <- nonuniqueCounts
        }
        colnames(geneMat) <- rownames(x_obj$sampleData)
        return(geneMat)
    })
    return(do.call(cbind, nonuniqueCountsList))
}

#' Generate gene counts directly from quantData
#' @noRd
generateGeneCountsFromQuantData <- function(quantData, annotations){
    # 1. Get unique transcript counts and aggregate them to gene level
    uniqueCounts <- generateUniqueCountsFromQuantData(quantData, annotations)
    
    geneIDs <- factor(mcols(annotations)$GENEID, levels = unique(mcols(annotations)$GENEID))
    geneCounts <- Matrix::fac2sparse(geneIDs) %*% uniqueCounts
    
    # 2. Get non-unique counts and incompatible counts
    nonuniqueCounts <- generateNonUniqueCountsFromQuantData(quantData, annotations)
    incompatibleCounts <- do.call(cbind, lapply(quantData, function(x) x$incompatibleCounts))
    
    # 3. Align the rows/dimensions and sum them up
    if(!is.null(incompatibleCounts)){
        incompatibleCounts <- Matrix(incompatibleCounts[match(rownames(geneCounts), rownames(incompatibleCounts)), ], sparse = TRUE)
        geneCounts <- geneCounts + incompatibleCounts
    }
    
    if(!is.null(nonuniqueCounts)){
        nonuniqueCounts <- Matrix(nonuniqueCounts[match(rownames(geneCounts), rownames(nonuniqueCounts)), ], sparse = TRUE)
        geneCounts <- geneCounts + nonuniqueCounts
    }
    
    return(geneCounts)
}
