
#' Isoform reconstruction using genomic alignments
#' @param readGrgList readGrgList
#' @param unlisted_junctions unlisted_junctions
#' @param uniqueJunctions uniqueJunctions
#' @param runName runName
#' @param annotations annotations
#' @param stranded stranded
#' @param verbose verbose
#' @inheritParams bambu
#' @noRd
isore.constructReadClasses <- function(readGrgList, unlisted_junctions,
                                       uniqueJunctions,
                                       runName = "sample1", 
                                       annotations, stranded = FALSE, 
                                       verbose = FALSE) {
    #split reads into single exon and multi exon reads
    reads.singleExon <- unlist(readGrgList[elementNROWS(readGrgList) == 1],
                                  use.names = FALSE)
    mcols(reads.singleExon)$id <- mcols(readGrgList[
        elementNROWS(readGrgList) == 1])$id
    
  ###### HERE #####  readGrgList <- readGrgList[elementNROWS(readGrgList) > 1]
    uniqueReadIds <- unique(mcols(unlisted_junctions)$id)
    if (any(order(uniqueReadIds) != seq_along(uniqueReadIds))) 
        warning("read Id not sorted, can result in wrong assignments.
            Please report error")
    #only keep multi exons reads in readGrgList
    readGrgList <- readGrgList[match(uniqueReadIds, mcols(readGrgList)$id)]
    
    start.ptm <- proc.time()
    readClassListSpliced <- constructSplicedReadClassTables(
        uniqueJunctions = uniqueJunctions,
        unlisted_junctions = unlisted_junctions,
        readGrgList = readGrgList,
        stranded = stranded)
    end.ptm <- proc.time()
    if (verbose)
        message("Finished create transcript models (read classes) for reads with
    spliced junctions in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
    
    exonsByReadClass <- generateExonsByReadClass(reads.singleExon, 
                                                 annotations, 
                                                 readClassListSpliced, 
                                                 stranded, verbose)
    counts <- matrix(mcols(exonsByReadClass)$readCount,
                     dimnames = list(names(exonsByReadClass), runName))
    colDataDf <- DataFrame(name = runName, row.names = runName)
    mcols(exonsByReadClass) <- mcols(exonsByReadClass)[, c("chr.rc", 
                                                           "strand.rc", "intronStarts", "intronEnds", "confidenceType")]
    se <- SummarizedExperiment(assays = SimpleList(counts = counts),
                               rowRanges = exonsByReadClass, colData = colDataDf)
    return(se)
}


#' reconstruct spliced transripts
#' @importFrom unstrsplit getFromNamespace
#' @noRd
constructSplicedReadClassTables <- function(uniqueJunctions, unlisted_junctions, 
                                            readGrgList, stranded = FALSE) {
    options(scipen = 999)
    allToUniqueJunctionMatch <- match(unlisted_junctions,
                                                uniqueJunctions, ignore.strand = TRUE)
    correctedJunctionMatches <- match(uniqueJunctions$mergedHighConfJunctionIdAll_noNA[allToUniqueJunctionMatch], names(uniqueJunctions))
    
    unlisted_junctions <- correctIntronRanges(unlisted_junctions, 
                                              uniqueJunctions, 
                                              correctedJunctionMatches)
    if (!stranded) {
        readStrand <- correctReadStrandById(as.character(strand(unlisted_junctions)),
                                            id=mcols(unlisted_junctions)$id)
    } else {
        readStrand <- as.character(getStrandFromGrList(readGrgList))
    }
    
    # confidence type (note: can be changed to integer encoding)
    readConfidence <- rep("highConfidenceJunctionReads", length(readStrand))
    lowConfidenceReads <- which(sum(is.na(splitAsList(
        uniqueJunctions$mergedHighConfJunctionId[allToUniqueJunctionMatch],
        mcols(unlisted_junctions)$id))) > 0)
    readConfidence[lowConfidenceReads] <- "lowConfidenceJunctionReads"
    
    readTable <- createReadTable(unlisted_junctions, readGrgList, readStrand, readConfidence)
    exonsByReadClass <- createExonsByReadClass(readTable)
    ## combine new transcripts with annotated transcripts
    ## based on identical intron pattern
    readTable <- readTable %>% dplyr::select(chr.rc = chr, strand.rc = strand,
                                             intronStarts, intronEnds, confidenceType, readCount)
    mcols(exonsByReadClass) <- readTable
    options(scipen = 0)
    return(exonsByReadClass)
}

#' this functions uses the uniqueJunction table which has reference junctions
#' and replaces intron coordinates with coordinates from the reference junction
#' the strand of junctions is also changed to the reference junction strand
#' @noRd
correctIntronRanges <- function(unlisted_junctions, uniqueJunctions,
                                correctedJunctionMatches){
    intronStartTMP <-
        start(uniqueJunctions)[correctedJunctionMatches]
    
    intronEndTMP <-
        end(uniqueJunctions)[correctedJunctionMatches]
    
    exon_0size <- 
        which(intronStartTMP[-1] <= intronEndTMP[-length(intronEndTMP)] &
                  mcols(unlisted_junctions)$id[-1] == 
                  mcols(unlisted_junctions)$id[-length(unlisted_junctions)])
    
    if (length(exon_0size))
        intronStartTMP[-1][exon_0size] <-
        intronEndTMP[-length(intronEndTMP)][exon_0size] + 1
    
    start(unlisted_junctions) <- intronStartTMP
    end(unlisted_junctions) <- intronEndTMP
    strand(unlisted_junctions) <- uniqueJunctions$strand.mergedHighConfJunction[correctedJunctionMatches]
    return(unlisted_junctions)
}

#' This function returns the inferred strand based on the number of +(plus) and
#' -(minus) junctions in each read (majority vote)
#' @noRd
correctReadStrandById <- function(strand, id, stranded=FALSE){
    plusCount <- as.integer(sum(splitAsList(strand, id) == "+"))
    minusCount <- as.integer(sum(splitAsList(strand, id) == "-"))
    strandJunctionSum <- minusCount - plusCount
    readStrand <- rep("*", length(strandJunctionSum))
    readStrand[strandJunctionSum < 0] <- "+"
    readStrand[strandJunctionSum > 0] <- "-"
    return(readStrand)
}


#' This function generates a table that contains 1 row for each (spliced) read
#' This table will then be summarised into read classes
#' @param unlisted_junctions unlisted_junctions
#' @param readGrgList reads GRangesList
#' @param readStrand readStrand
#' @param readConfidence readConfidence
#' allToUniqueJunctionMatch
#' @noRd
createReadTable <- function(unlisted_junctions, readGrgList,
                            readStrand, readConfidence) {
    readTable <- as_tibble(data.frame(matrix(ncol = 7, 
        nrow = length(readGrgList))))
    colnames(readTable) <- c("chr", "start", "end", "strand", "intronEnds",
        "intronStarts", "confidenceType")
    readTable[, "chr"] <- as.factor(getChrFromGrList(readGrgList))
    # intron start and end
    readRanges <- unlist(range(ranges(readGrgList)), use.names = FALSE)
    readTable[, "intronStarts"] <- unstrsplit(splitAsList(as.character(
        start(unlisted_junctions)), mcols(unlisted_junctions)$id), sep = ",")
    readTable[, "intronEnds"] <- unstrsplit(splitAsList(as.character(
        end(unlisted_junctions)), mcols(unlisted_junctions)$id), sep = ",")
    # start and end
    intronStartCoordinatesInt <- as.integer(min(splitAsList(start(unlisted_junctions),
            mcols(unlisted_junctions)$id)) - 2)
    intronEndCoordinatesInt <- as.integer(max(splitAsList(end(unlisted_junctions),
            mcols(unlisted_junctions)$id)) + 2)
    readTable[, "start"] <- pmin(start(readRanges), intronStartCoordinatesInt)
    readTable[, "end"] <- pmax(end(readRanges), intronEndCoordinatesInt)
    readTable[, "strand"] <- readStrand
    # confidence type (note: can be changed to integer encoding)
    readTable[, "confidenceType"] <- readConfidence

    ## currently the 80% and 20% quantile of reads is used to 
    ## identify start and end sites
    readTable <- readTable %>% group_by(chr, strand, intronEnds,
        intronStarts) %>% summarise( readCount = n(), start = nth(
                x = start, n = ceiling(readCount / 5), order_by = start),
            end = nth(x = end, n = ceiling(readCount / 1.25), order_by = end),
            confidenceType = dplyr::first(confidenceType)) %>%
        ungroup() %>% arrange(chr, start, end)
    readTable$readClassId <- paste("rc", seq_len(nrow(readTable)), sep = ".")
    return(readTable)
}

#' @noRd
getChrFromGrList <- function(grl) { 
    return(unlist(seqnames(grl), use.names = FALSE)[cumsum(elementNROWS(grl))]) 
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



#' generate exonByReadClass
#' @noRd
generateExonsByReadClass <- function(reads.singleExon, annotationGrangesList, readClassListSpliced, stranded, verbose){
    
    start.ptm <- proc.time()
    # singleExonReads <- unlist(readGrgList[elementNROWS(readGrgList) == 1],
    #                           use.names = FALSE)
    # mcols(singleExonReads)$id <- mcols(readGrgList[
    #     elementNROWS(readGrgList) == 1])$id
    referenceExons <- unique(c(GenomicRanges::granges(unlist(
        readClassListSpliced[mcols(readClassListSpliced)$confidenceType ==
                                 "highConfidenceJunctionReads" &
                                 mcols(readClassListSpliced)$strand.rc != "*"], use.names = FALSE)), 
        GenomicRanges::granges(unlist(annotationGrangesList,
                                      use.names = FALSE))))
    readClassListUnsplicedWithAnnotation <- constructUnsplicedReadClasses(
        granges = reads.singleExon, grangesReference = referenceExons,
        confidenceType = "unsplicedWithin", stranded = stranded)
    reads.singleExon <- reads.singleExon[!mcols(reads.singleExon)$id %in%
                                           readClassListUnsplicedWithAnnotation$readIds]
    referenceExons <- reduce(reads.singleExon, ignore.strand = !stranded)
    readClassListUnsplicedReduced <- constructUnsplicedReadClasses(
        granges = reads.singleExon, grangesReference = referenceExons,
        confidenceType = "unsplicedNew", stranded = stranded)
    end.ptm <- proc.time()
    if (verbose) message("Finished create single exon transcript models
        (read classes) in ", round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    exonsByReadClass <- c(readClassListSpliced,
                          readClassListUnsplicedWithAnnotation$exonsByReadClass,
                          readClassListUnsplicedReduced$exonsByReadClass)
    return(exonsByReadClass)
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
    hitsDF <- hitsDF %>% dplyr::select(queryHits, chr, start, end, strand) %>%
        group_by(queryHits, chr, strand) %>%
        summarise(start = max(start), end = min(end)) %>%
        group_by(chr, start, end, strand) %>%
        mutate(readClassId = paste0("rc", confidenceType, ".",
        cur_group_id())) %>% ungroup()
    readIds <- mcols(granges[hitsDF$queryHits])$id
    hitsDF <- hitsDF %>% dplyr::select(chr, start, end, strand, readClassId) %>%
        group_by(readClassId) %>% mutate(readCount = n()) %>% distinct() %>%
        ungroup() %>% mutate(confidenceType = confidenceType, intronStarts = NA,
            intronEnds = NA) %>% dplyr::select(
            chr, start, end, strand, intronStarts,
            intronEnds, confidenceType, readClassId, readCount)
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
        intronStarts, intronEnds, confidenceType, readCount)
    mcols(exByReadClassUnspliced) <- hitsDF
    # seqlevels(exByReadClassUnspliced) <- seqLevelList
    return(list(exonsByReadClass = exByReadClassUnspliced, readIds = readIds))
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
