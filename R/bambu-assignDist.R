#' Create equivilence classes and assign to transcripts
#' @inheritParams bambu
#' @import data.table
#' @noRd
assignReadClasstoTranscripts <- function(readClassList, annotations, isoreParameters, 
                                        verbose, demultiplexed, spatial, 
                                        returnDistTable = FALSE, trackReads = TRUE) {
    metadata(readClassList)$readClassDist <- calculateDistTable(readClassList, annotations, isoreParameters, verbose)
    readClassList = splitReadClassFiles(readClassList)
    readClassDt <- genEquiRCs(metadata(readClassList)$readClassDist, annotations, verbose) 
    readClassDt$eqClass.match = match(readClassDt$eqClassById,metadata(readClassList)$eqClassById)
    readClassDt <- simplifyNames(readClassDt)
    readClassDt = readClassDt %>% group_by(eqClassId, gene_sid) %>% 
        mutate(multi_align = length(unique(txid))>1) %>% ungroup() %>% mutate(aval = 1) %>%
        data.table()

    #return non-em counts
    ColData = generateColData(colnames(metadata(readClassList)$countMatrix), clusters = NULL, demultiplexed, spatial)
    quantData <- SummarizedExperiment(assays = SimpleList(
        counts = generateUniqueCounts(readClassDt, metadata(readClassList)$countMatrix, annotations)),
        rowRanges = annotations,
        colData = ColData)
    colnames(quantData) = ColData$id
    if(sum(metadata(readClassList)$incompatibleCountMatrix)==0){
        metadata(quantData)$incompatibleCounts = NULL
    } else{
        metadata(quantData)$incompatibleCounts = generateIncompatibleCounts(metadata(readClassList)$incompatibleCountMatrix, annotations)       
    }
    metadata(quantData)$nonuniqueCounts = generateNonUniqueCounts(readClassDt, metadata(readClassList)$countMatrix, annotations)
    metadata(quantData)$readClassDt = readClassDt
    metadata(quantData)$countMatrix = metadata(readClassList)$countMatrix
    metadata(quantData)$incompatibleCountMatrix  = metadata(readClassList)$incompatibleCountMatrix 
    metadata(quantData)$sampleNames = metadata(readClassList)$sampleNames 
    if(returnDistTable){
        metadata(quantData)$distTable = metadata(readClassList)$readClassDist
    }
    if (trackReads){
        metadata(quantData)$readToTranscriptMap = 
            generateReadToTranscriptMap(readClassList, metadata(readClassList)$readClassDist, 
                                  annotations)
    } 
    return(quantData)     

}


generateUniqueCounts <- function(readClassDt, countMatrix, annotations){
    x = readClassDt %>% filter(!multi_align & !is.na(eqClass.match))
    uniqueCounts = countMatrix[x$eqClass.match,]
    uniqueCounts.tx = sparse.model.matrix(~ factor(x$txid) - 1)
    uniqueCounts = t(uniqueCounts.tx) %*% uniqueCounts
    rownames(uniqueCounts) = names(annotations)[match(as.numeric(levels(factor(x$txid))),mcols(annotations)$txid)]
    counts = sparseMatrix(length(annotations), ncol(uniqueCounts), x = 0)
    rownames(counts) = names(annotations)
    counts[rownames(uniqueCounts),] = uniqueCounts
    return(counts)

    counts.total = colSums(countMatrix) + colSums(incompatibleCountMatrix)
    counts.total[counts.total==0] = 1
    counts.CPM = counts/counts.total * 10^6

}

generateIncompatibleCounts <- function(incompatibleCountMatrix, annotations){
    genes = levels(factor(unique(mcols(annotations)$GENEID)))
    rownames(incompatibleCountMatrix) = genes[as.numeric(rownames(incompatibleCountMatrix))]
    geneMat = sparseMatrix(length(genes), ncol(incompatibleCountMatrix), x = 0)
    rownames(geneMat) = genes
    geneMat[rownames(incompatibleCountMatrix),] = incompatibleCountMatrix
    return(geneMat)
}

generateNonUniqueCounts <- function(readClassDt, countMatrix, annotations){
    #fuse multi align RCs by gene
    x = readClassDt %>% filter(multi_align & !is.na(eqClass.match))
    x = x %>% distinct(eqClassId, .keep_all = TRUE)
    nonuniqueCounts = countMatrix[x$eqClass.match,, drop = FALSE]
    if(nrow(x)>1){
    nonuniqueCounts.gene = sparse.model.matrix(~ factor(x$gene_sid) - 1)
    nonuniqueCounts = t(nonuniqueCounts.gene) %*% nonuniqueCounts
    }
    #covert ids into gene ids
    geneids = as.numeric(levels(factor(x$gene_sid)))
    geneids = x$txid[match(geneids, x$gene_sid)]
    geneids = mcols(annotations)$GENEID[as.numeric(geneids)]
    rownames(nonuniqueCounts) = geneids
    #create matrix for all annotated genes
    genes = levels(factor(unique(mcols(annotations)$GENEID)))
    geneMat = sparseMatrix(length(genes), ncol(nonuniqueCounts), x = 0)
    rownames(geneMat) = genes
    geneMat[rownames(nonuniqueCounts),] = nonuniqueCounts
    return(geneMat)
}
