
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
        mcols(annotations)$TXNAME)]
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
