
#' Generate nonunique count matrix (gene x sample) from multi-align reads in readClassDt
#' @noRd
generateNonUniqueCountMatrix <- function(readClassDt, annotations, sampleIds){
    genes <- unique(mcols(annotations)$GENEID)
    nSamples <- length(sampleIds)
    geneMat <- sparseMatrix(length(genes), nSamples, x = 0)
    rownames(geneMat) <- genes
    colnames(geneMat) <- sampleIds

    x <- readClassDt %>% filter(multi_align & !is.na(eqClass.match))
    x <- x %>% distinct(eqClassId, .keep_all = TRUE)
    if (nrow(x) == 0) return(geneMat)

    i <- rep(seq_along(x$gene_sid), lengths(x$columnIds))
    j <- unlist(x$columnIds)
    x_vals <- unlist(x$columnCounts)
    nonuniqueCounts <- sparseMatrix(i = i, j = j, x = x_vals, dims = c(nrow(x), nSamples))

    if (nrow(x) > 1 & length(unique(x$gene_sid)) > 1) {
        nonuniqueCounts.gene <- sparse.model.matrix(~ factor(x$gene_sid) - 1)
        nonuniqueCounts <- t(nonuniqueCounts.gene) %*% nonuniqueCounts
    } else {
        nonuniqueCounts.gene <- Matrix(1, nrow = nrow(x), ncol = 1, sparse = TRUE)
        nonuniqueCounts <- t(nonuniqueCounts.gene) %*% nonuniqueCounts
    }

    geneids <- as.numeric(levels(factor(x$gene_sid)))
    geneids <- x$txid[match(geneids, x$gene_sid)]
    geneids <- mcols(annotations)$GENEID[as.numeric(geneids)]
    rownames(nonuniqueCounts) <- geneids
    colnames(nonuniqueCounts) <- sampleIds

    geneMat[match(rownames(nonuniqueCounts), rownames(geneMat)), ] <- nonuniqueCounts
    return(geneMat)
}

#' rename runnames when there are duplicated names
#' @title rename_duplicatedNames
#' @param runnames sample names
#' @noRd
rename_duplicatedNames <- function(runnames){
    ## rename runnames when duplicated names are found
    if (length(which(duplicated(runnames)))) {
        iter <- 1
        while (length(which(duplicated(runnames)))) {
            if (iter == 1) {
                runnames[which(duplicated(runnames))] <-
                    paste0(runnames[which(duplicated(runnames))], "...", iter)
            } else {
                runnames[which(duplicated(runnames))] <-
                    gsub(paste0("...", iter - 1, "$"), paste0("...", iter),
                    runnames[which(duplicated(runnames))])
            }
            iter <- iter + 1
        }
    }
    return(runnames)
}


#' From tx ranges to gene ranges
#' @importFrom GenomicRanges reduce 
#' @noRd
reducedRangesByGenes <- function(annotations) {
    annotations <- annotations[order(mcols(annotations)$GENEID)]
    unlistData <- unlist(annotations)
    geneIds <- mcols(annotations)$GENEID[match(names(unlistData), 
        names(annotations))]
    partitioning <- PartitioningByEnd(cumsum(table(geneIds)),
                                    names = NULL)
    exonsByGene <- relist(unlistData, partitioning)
    exonsByGeneReduced <- reduce(exonsByGene)
    return(exonsByGeneReduced)
}

# From tx ranges to gene ranges
# txRangesToGeneRanges <- function(exByTx, TXNAMEGENEID_Map) {
#     # rename names to geneIDs
#     names(exByTx) <- as.data.table(TXNAMEGENEID_Map)[match(names(exByTx),
#         TXNAME)]$GENEID
# 
#     # combine gene exon ranges and reduce overlapping ones
#     unlistData <- unlist(exByTx, use.names = TRUE)
#     orderUnlistData <- unlistData[order(names(unlistData))]
#     orderUnlistData$exon_rank <- NULL
#     orderUnlistData$exon_endRank <- NULL
# 
#     exByGene <- splitAsList(orderUnlistData, names(orderUnlistData))
#     exByGene <- reduce(exByGene)
# 
#     # add exon_rank and endRank
#     unlistData <- unlist(exByGene, use.names = FALSE)
#     partitionDesign <- cumsum(elementNROWS(exByGene))
#     partitioning <- PartitioningByEnd(partitionDesign, names = NULL)
#     geneStrand <- as.character(strand(unlistData))[partitionDesign]
#     exon_rank <- lapply(width((partitioning)), seq, from = 1)
#     exon_rank[which(geneStrand == "-")] <-
#         lapply(exon_rank[which(geneStrand == "-")], rev)
#     # * assumes positive for exon ranking
#     exon_endRank <- lapply(exon_rank, rev)
#     unlistData$exon_rank <- unlist(exon_rank)
#     unlistData$exon_endRank <- unlist(exon_endRank)
#     exByGene <- relist(unlistData, partitioning)
# 
#     return(exByGene)
# }
