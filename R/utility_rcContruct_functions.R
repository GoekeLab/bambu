
#' calculate distance between first and last exon matches
#' @param uniqueJunctions uniqueJunctions
#' @param unlisted_junctions unlisted_junctions
#' @param readGrgList reads GRangesList
#' @param firstseg firstseg
#' @param intronStartTMP intronStartTMP
#' @param intronEndTMP intronEndTMP
#' @param readStrand readStrand
#' @param allJunctionToUniqueJunctionOverlap
#' allJunctionToUniqueJunctionOverlap
#' @noRd
createReadTable <- function(uniqueJunctions, unlisted_junctions, readGrgList,
                            firstseg, intronStartTMP, intronEndTMP, readStrand,
                            allJunctionToUniqueJunctionOverlap) {
    readTable <- as_tibble(data.frame(matrix(ncol = 7, 
        nrow = length(readGrgList))))
    colnames(readTable) <- c("chr", "start", "end", "strand", "intronEnds",
        "intronStarts", "confidenceType")
    # chr
    readTable[, "chr"] <- as.factor(seqnames(unlist(readGrgList)[firstseg]))
    # intron start and end
    readRanges <- unlist(range(ranges(readGrgList)), use.names = FALSE)
    readTable[, "intronStarts"] <- unstrsplit(splitAsList(as.character(
        intronStartTMP), mcols(unlisted_junctions)$id), sep = ",")
    readTable[, "intronEnds"] <- unstrsplit(splitAsList(as.character(
        intronEndTMP), mcols(unlisted_junctions)$id), sep = ",")
    # start and end
    intronStartCoordinatesInt <- as.integer(min(splitAsList(intronStartTMP,
            mcols(unlisted_junctions)$id)) - 2)
    intronEndCoordinatesInt <- as.integer(max(splitAsList(intronEndTMP,
            mcols(unlisted_junctions)$id)) + 2)
    readTable[, "start"] <- pmin(start(readRanges), intronStartCoordinatesInt)
    readTable[, "end"] <- pmax(end(readRanges), intronEndCoordinatesInt)
    # strand
    readTable[, "strand"] <- readStrand
    # confidence type (note: can be changed to integer encoding)
    readTable[, "confidenceType"] <- "highConfidenceJunctionReads"
    lowConfidenceReads <- which(sum(is.na(splitAsList(
            uniqueJunctions$mergedHighConfJunctionId[subjectHits(
            allJunctionToUniqueJunctionOverlap)],
            mcols(unlisted_junctions)$id))) > 0)
    ## currently the 80% and 20% quantile of reads is used to 
    ## identify start and end sites
    readTable[lowConfidenceReads, "confidenceType"] <- 
        "lowConfidenceJunctionReads"
    indices=group_indices(readTable %>% group_by(chr, strand, 
        intronEnds, intronStarts))
    names(indices) = mcols(readGrgList)$id
    readTable <- readTable %>% 
        group_by( chr, strand, intronEnds, intronStarts) %>% 
        summarise(readCount = n(), start = nth(
            x = start, n = ceiling(readCount / 5), order_by = start),
            end = nth(x = end, n = ceiling(readCount / 1.25), order_by = end),
            confidenceType = dplyr::first(confidenceType)) %>%
        ungroup()# %>% arrange(chr, start, end)
    readTable$readClassId <- paste("rc", seq_len(nrow(readTable)), sep = ".")
    return(list(readTable = readTable, indices = indices))
}

#' reconstruct spliced transripts
#' @importFrom unstrsplit getFromNamespace
#' @noRd
constructSplicedReadClassTables <- function(uniqueJunctions,
                        unlisted_junctions, readGrgList, stranded = FALSE) {
    options(scipen = 999)
    uniqueReadIds <- unique(mcols(unlisted_junctions)$id)
    if (any(order(uniqueReadIds) != seq_along(uniqueReadIds))) 
        warning("read Id not sorted, can result in wrong assignments.
            Please report error")
    readGrgList <- readGrgList[match(uniqueReadIds, mcols(readGrgList)$id)]
    firstseg <- start(PartitioningByWidth(readGrgList))
    allJunctionToUniqueJunctionOverlap <- findOverlaps(unlisted_junctions,
        uniqueJunctions, type = "equal", ignore.strand = TRUE)
    intronStartTMP <- createIntronTmp(uniqueJunctions,
    allJunctionToUniqueJunctionOverlap,unlisted_junctions)[[1]]
    intronEndTMP <- createIntronTmp(uniqueJunctions,
    allJunctionToUniqueJunctionOverlap,unlisted_junctions)[[2]]
    if (!stranded) {
        readStrand <- correctReadTableStrand(uniqueJunctions,
            unlisted_junctions, allJunctionToUniqueJunctionOverlap)
    }else{
        readStrand <- as.character(strand(unlist(readGrgList)[firstseg]))
    }
    createReadTableOutput <- createReadTable(
        uniqueJunctions, unlisted_junctions, readGrgList,
        firstseg, intronStartTMP, intronEndTMP, readStrand,
        allJunctionToUniqueJunctionOverlap)
    readTable = createReadTableOutput$readTable
    indices = createReadTableOutput$indices
    exonsByReadClass <- createExonsByReadClass(readTable)
    ## combine new transcripts with annotated transcripts
    ## based on identical intron pattern
    readTable <- readTable %>% dplyr::select(chr.rc = chr, strand.rc = strand,
        start.rc = start, end.rc = end, intronStarts, 
        intronEnds, confidenceType, readCount)
    readTable$readClassId <- paste("rc", seq_len(nrow(readTable)), sep = ".")
    mcols(exonsByReadClass) <- readTable
    options(scipen = 0)
    return(list(readClassList = exonsByReadClass, indices = indices))
}

#' @noRd
correctReadTableStrand <- function(uniqueJunctions,
    unlisted_junctions, allJunctionToUniqueJunctionOverlap){
    
        unlisted_junctions_strand <-
            uniqueJunctions$strand.mergedHighConfJunction[subjectHits(
                allJunctionToUniqueJunctionOverlap)]
        plusCount <- as.integer(sum(splitAsList( unlisted_junctions_strand,
            mcols(unlisted_junctions)$id) == "+"))
        minusCount <- as.integer(sum(splitAsList(unlisted_junctions_strand,
            mcols(unlisted_junctions)$id) == "-"))
        strandJunctionSum <- minusCount - plusCount
        readStrand <- rep("*", length(strandJunctionSum))
        readStrand[strandJunctionSum < 0] <- "+"
        readStrand[strandJunctionSum > 0] <- "-"
    return(readStrand)
}


#' @noRd
createIntronTmp <- function(uniqueJunctions,
    allJunctionToUniqueJunctionOverlap,unlisted_junctions){
    intronStartTMP <-
        start(uniqueJunctions[uniqueJunctions$mergedHighConfJunctionIdAll_noNA[
            subjectHits(allJunctionToUniqueJunctionOverlap)]])
    
    intronEndTMP <-
        end(uniqueJunctions[uniqueJunctions$mergedHighConfJunctionIdAll_noNA[
            subjectHits(allJunctionToUniqueJunctionOverlap)]])
    
    exon_0size <- 
        which(intronStartTMP[-1] <= intronEndTMP[-length(intronEndTMP)] &
        mcols(unlisted_junctions)$id[-1] == 
            mcols(unlisted_junctions)$id[-length(unlisted_junctions)])
    
    if (length(exon_0size))
        intronStartTMP[-1][exon_0size] <-
            intronEndTMP[-length(intronEndTMP)][exon_0size] + 1
    return(list(intronStartTMP, intronEndTMP))
}

#' @noRd
createExonsByReadClass <- function(readTable){
    exonsByReadClass <- GenomicRanges::makeGRangesListFromFeatureFragments(
        seqnames = readTable$chr, fragmentStarts = paste(readTable$start - 1,
            readTable$intronEnds,sep = ","), 
            fragmentEnds = paste(readTable$intronStarts, readTable$end + 1,
                sep = ","), strand = readTable$strand)
    exonsByReadClass <- narrow(exonsByReadClass, start = 2, end = -2)
    
    # correct junction to exon differences in coordinates
    names(exonsByReadClass) <- readTable$readClassId
    unlistData <- unlist(exonsByReadClass, use.names = FALSE)
    partitioning <- PartitioningByEnd(cumsum(elementNROWS(exonsByReadClass)),
        names = NULL)
    exon_rank <- lapply(width((partitioning)), seq, from = 1)
    # add exon rank and exon_endRank
    exon_rank[which(readTable$strand == "-")] <-
        lapply(exon_rank[which(readTable$strand == "-")], rev)
    # * assumes positive for exon ranking
    exon_endRank <- lapply(exon_rank, rev)
    unlistData$exon_rank <- unlist(exon_rank)
    unlistData$exon_endRank <- unlist(exon_endRank)
    exonsByReadClass <- relist(unlistData, partitioning)
    
    return(exonsByReadClass)
}

#' initiate the hits dataframe
#' @param hitsWithin hitsWithin
#' @param grangesReference grangesReference
#' @param stranded stranded
#' @noRd
initiateHitsDF <- function(hitsWithin, grangesReference, stranded) {
    hitsDF <- as_tibble(hitsWithin)
    hitsDF$chr <-
        as.factor(seqnames(grangesReference)[subjectHits(hitsWithin)])
    hitsDF$start <- start(grangesReference)[subjectHits(hitsWithin)]
    hitsDF$end <- end(grangesReference)[subjectHits(hitsWithin)]
    if (stranded == FALSE) {
        hitsDF$strand <- "*"
    } else {
        hitsDF$strand <-
            as.character(strand(grangesReference)[subjectHits(hitsWithin)])
    }
    return(hitsDF)
}

#' reconstruct read classes using unspliced reads that fall
#' within exons from annotations
#' @noRd
constructUnsplicedReadClasses <- function(granges, grangesReference,
    confidenceType = "unspliced", stranded = TRUE) {
    if (is.null(mcols(granges)$id))
        stop("ID column is missing from mcols(granges)")
    # seqLevelList <- unique(c(seqlevels(granges),
    # seqlevels(grangesReference)))
    hitsWithin <- findOverlaps(granges, grangesReference,
        ignore.strand = !stranded, type = "within", select = "all")
    # find reads that overlap with reference ranges
    hitsDF <- initiateHitsDF(hitsWithin, grangesReference, stranded)
    ## create single exon read class by using the minimum end
    ## and maximum start of all overlapping exons (identical to
    ## minimum equivalent class)
    indices2=group_indices(hitsDF %>% group_by(queryHits, chr, strand))
    names(indices2) = mcols(granges[hitsDF$queryHits])$id
    hitsDF <- hitsDF %>% dplyr::select(queryHits, chr, start, end, strand) %>%
        group_by(queryHits, chr, strand) %>%
        summarise(start = max(start), end = min(end)) %>%
        group_by(chr, start, end, strand) %>%
        mutate(readClassId = paste0("rc", confidenceType, ".",
        cur_group_id())) %>% ungroup()
    indices=group_indices(hitsDF %>% group_by(readClassId))
    names(indices) = mcols(granges[hitsDF$queryHits])$id
    indices = indices[indices2]
    names(indices) = names(indices2)
    
    hitsDF <- hitsDF %>% dplyr::select(chr, start, end, strand, 
        readClassId, queryHits) %>%
        group_by(readClassId) %>% mutate(readCount = n()) %>% distinct() %>%
        ungroup() %>% mutate(confidenceType = confidenceType, 
            intronStarts = NA, intronEnds = NA, queryHits = queryHits) %>% 
        dplyr::select(chr, start, end, strand, intronStarts,
            intronEnds, confidenceType, readClassId, readCount, queryHits)
    exByReadClassUnspliced <- GenomicRanges::GRanges(
        seqnames = hitsDF$chr,
        ranges = IRanges(
            start = hitsDF$start,
            end = hitsDF$end),
        strand = hitsDF$strand
    )
    exByReadClassUnspliced$exon_rank <- 1
    exByReadClassUnspliced$exon_endRank <- 1
    partitioning <- PartitioningByEnd(seq_along(exByReadClassUnspliced))
    exByReadClassUnspliced <- relist(exByReadClassUnspliced, partitioning)
    names(exByReadClassUnspliced) <- hitsDF$readClassId
    hitsDF <- dplyr::select(hitsDF, chr.rc = chr, strand.rc = strand,
        start.rc = start, end.rc = end, intronStarts, intronEnds, 
        confidenceType, readCount)
    mcols(exByReadClassUnspliced) <- hitsDF
    return(list(readClasses = exByReadClassUnspliced, indices = indices))
}