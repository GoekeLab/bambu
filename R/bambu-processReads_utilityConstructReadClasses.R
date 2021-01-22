
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
    uniqueJunctions, runName = "sample1",
    annotations, stranded = FALSE, verbose = FALSE) {
    #split reads into single exon and multi exon reads
    reads.singleExon <- unlist(readGrgList[elementNROWS(readGrgList) == 1],
        use.names = FALSE)
    mcols(reads.singleExon)$id <- mcols(readGrgList[
        elementNROWS(readGrgList) == 1])$id
    #only keep multi exons reads in readGrgList    
    readGrgList <- readGrgList[elementNROWS(readGrgList) > 1]
    if (!identical(mcols(readGrgList)$id,unique(mcols(unlisted_junctions)$id))) 
        warning("read Id not sorted, can result in wrong assignments.
            Please report error")

    start.ptm <- proc.time()
    exonsByRC.spliced <- constructSplicedReadClasses(
        uniqueJunctions = uniqueJunctions,
        unlisted_junctions = unlisted_junctions,
        readGrgList = readGrgList,
        stranded = stranded)
    end.ptm <- proc.time()
    if (verbose) 
    message("Finished create transcript models (read classes) for reads with
    spliced junctions in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")

    exonsByRC.unspliced <- constructUnsplicedReadClasses(reads.singleExon, 
        annotations, exonsByRC.spliced, stranded, verbose)
    exonsByRC <- c(exonsByRC.spliced, exonsByRC.unspliced)
    counts <- matrix(mcols(exonsByRC)$readCount,
        dimnames = list(names(exonsByRC), runName))
    colDataDf <- DataFrame(name = runName, row.names = runName)
    mcols(exonsByRC) <- mcols(exonsByRC)[, c("chr.rc", "strand.rc", 
        "intronStarts", "intronEnds", "confidenceType")]
    se <- SummarizedExperiment(assays = SimpleList(counts = counts),
        rowRanges = exonsByRC, colData = colDataDf)
    return(se)
}


#' reconstruct spliced transripts
#' @importFrom S4Vectors unstrsplit 
#' @importFrom dplyr select %>%
#' @importFrom GenomicRanges match
#' @noRd
constructSplicedReadClasses <- function(uniqueJunctions, unlisted_junctions, 
    readGrgList, stranded = FALSE) {
    options(scipen = 999)
    allToUniqueJunctionMatch <- GenomicRanges::match(unlisted_junctions,
        uniqueJunctions, ignore.strand = TRUE)
    correctedJunctionMatches <- 
        base::match(uniqueJunctions$mergedHighConfJunctionIdAll_noNA[
            allToUniqueJunctionMatch], names(uniqueJunctions))

    unlisted_junctions <- correctIntronRanges(unlisted_junctions, 
        uniqueJunctions, correctedJunctionMatches)
    if (isFALSE(stranded)) {
        readStrand <- correctReadStrandById(
            as.factor(strand(unlisted_junctions)),
            id = mcols(unlisted_junctions)$id)
    } else {
        readStrand <- as.factor(getStrandFromGrList(readGrgList))
    }

    # confidence type (note: can be changed to integer encoding)
    readConfidence <- factor(rep("highConfidenceJunctionReads",
        length(readStrand)), levels = c('highConfidenceJunctionReads',
        'lowConfidenceJunctionReads'))
    lowConfidenceReads <- which(sum(is.na(splitAsList(
        uniqueJunctions$mergedHighConfJunctionId[allToUniqueJunctionMatch],
        mcols(unlisted_junctions)$id))) > 0)
    readConfidence[lowConfidenceReads] <- "lowConfidenceJunctionReads"

    readTable <- createReadTable(unlisted_junctions, readGrgList,
        readStrand, readConfidence)
    exonsByReadClass <- createExonsByReadClass(readTable)
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
    strand(unlisted_junctions) <-
        uniqueJunctions$strand.mergedHighConfJunction[correctedJunctionMatches]
    return(unlisted_junctions)
}

#' This function returns the inferred strand based on the number of +(plus) and
#' -(minus) junctions in each read (majority vote)
#' @noRd
correctReadStrandById <- function(strand, id, stranded = FALSE){
    plusCount <- as.integer(sum(splitAsList(strand, id) == "+"))
    minusCount <- as.integer(sum(splitAsList(strand, id) == "-"))
    strandJunctionSum <- minusCount - plusCount
    readStrand <- 
        factor(rep("*", length(strandJunctionSum)), levels = c('+','-','*'))
    readStrand[strandJunctionSum < 0] <- "+"
    readStrand[strandJunctionSum > 0] <- "-"
    return(readStrand)
}


#' This function generates a table that contains 1 row for each (spliced) read
#' This table is then summarised into read classes (identical junction patterns)
#' The readClass table is returned
#' @param unlisted_junctions unlisted_junctions
#' @param readGrgList reads GRangesList
#' @param readStrand readStrand
#' @param readConfidence readConfidence
#' @importFrom dplyr tibble %>% group_by n summarise nth order_by arrange mutate
#'     row_number .groups
#' @noRd
createReadTable <- function(unlisted_junctions, readGrgList,
    readStrand, readConfidence) {
    readRanges <- unlist(range(ranges(readGrgList)), use.names = FALSE)
    intronStartCoordinatesInt <- 
        as.integer(min(splitAsList(start(unlisted_junctions),
        mcols(unlisted_junctions)$id)) - 2)
    intronEndCoordinatesInt <- 
        as.integer(max(splitAsList(end(unlisted_junctions),
        mcols(unlisted_junctions)$id)) + 2)
    
    readTable <- tibble(chr = as.factor(getChrFromGrList(readGrgList)), 
        intronStarts = 
            unstrsplit(splitAsList(as.character(start(unlisted_junctions)),
                mcols(unlisted_junctions)$id), sep = ","),
        intronEnds = 
            unstrsplit(splitAsList(as.character(end(unlisted_junctions)),
                mcols(unlisted_junctions)$id), sep = ","),
        start = pmin(start(readRanges), intronStartCoordinatesInt),
        end = pmax(end(readRanges), intronEndCoordinatesInt),
        strand = readStrand,
        confidenceType = readConfidence)
    ## currently 80%/20% quantile of reads is used to identify start/end sites
    readTable <- readTable %>% 
        group_by(chr, strand, intronEnds, intronStarts, confidenceType) %>% 
        summarise(readCount = n(),
        start = nth(x = start, n = ceiling(readCount / 5), order_by = start),
        end = nth(x = end, n = ceiling(readCount / 1.25), order_by = end),
        .groups = 'drop') %>% arrange(chr, start, end) %>%
        mutate(readClassId = paste("rc", row_number(), sep = "."))
    return(readTable)
}

#' @noRd
getChrFromGrList <- function(grl) { 
    return(unlist(seqnames(grl), use.names = FALSE)[cumsum(elementNROWS(grl))]) 
}

#' Create Exons By Read Class
#' @importFrom GenomicRanges makeGRangesListFromFeatureFragments narrow 
#' @noRd
createExonsByReadClass <- function(readTable){
    exonsByReadClass <- makeGRangesListFromFeatureFragments(
        seqnames = readTable$chr,
        fragmentStarts = 
            paste(readTable$start - 1, readTable$intronEnds, sep = ","),
        fragmentEnds = 
            paste(readTable$intronStarts, readTable$end + 1, sep = ","),
        strand = readTable$strand)
    # correct junction to exon differences in coordinates
    exonsByReadClass <- narrow(exonsByReadClass, start = 2, end = -2)
    names(exonsByReadClass) <- readTable$readClassId
    
    # add exon rank and exon_endRank
    unlistData <- unlist(exonsByReadClass, use.names = FALSE)
    partitioning <- PartitioningByEnd(cumsum(elementNROWS(exonsByReadClass)),
        names = NULL)
    exon_rank <- lapply(width(partitioning), seq, from = 1)
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
#' @importFrom GenomicRanges granges unlist reduce 
#' @noRd
constructUnsplicedReadClasses <- function(reads.singleExon, annotations, 
    readClassListSpliced, stranded, verbose){

    start.ptm <- proc.time()
    referenceExons <- unique(c(granges(unlist(
        readClassListSpliced[mcols(readClassListSpliced)$confidenceType ==
            "highConfidenceJunctionReads" &
            mcols(readClassListSpliced)$strand.rc != "*"], 
        use.names = FALSE)),
        granges(unlist(annotations, use.names = FALSE))))
    #(1) reads fall into annotations or spliced read classes are summarised
    # by their minimum read class coordinates
    rcUnsplicedAnnotation <- getUnsplicedReadClassByReference(
        granges = reads.singleExon, grangesReference = referenceExons,
        confidenceType = "unsplicedWithin", stranded = stranded)
    reads.singleExon <- reads.singleExon[!mcols(reads.singleExon)$id %in%
        rcUnsplicedAnnotation$readIds]
    referenceExons <- reduce(reads.singleExon, ignore.strand = !stranded)
    #(2) reads do not fall within a annotated exon/high confidence read class 
    # exon are summarised based on the union of overlapping unspliced reads
    rcUnsplicedReduced <- getUnsplicedReadClassByReference(
        granges = reads.singleExon, grangesReference = referenceExons,
        confidenceType = "unsplicedNew", stranded = stranded)
    end.ptm <- proc.time()
    if (verbose) message("Finished create single exon transcript models
        (read classes) in ", round((end.ptm - start.ptm)[3] / 60, 1), " mins.")
    exonsByReadClass <- c(rcUnsplicedAnnotation$exonsByReadClass,
        rcUnsplicedReduced$exonsByReadClass)
    return(exonsByReadClass)
}



#' reconstruct read classes using unspliced reads that fall
#' within exons from annotations
#' @importFrom GenomicRanges GRanges relist
#' @importFrom dplyr %>% select group_by summarise .groups mutate cut_group_id
#'     ungroup distinct 
#' @noRd
getUnsplicedReadClassByReference <- function(granges, grangesReference,
    confidenceType = "unspliced", stranded = TRUE) {
    if (is.null(mcols(granges)$id))
        stop("ID column is missing from mcols(granges)")
    hitsWithin <- findOverlaps(granges, grangesReference,
        ignore.strand = !stranded, 
        type = "within", 
        select = "all")
    # find reads that overlap with reference ranges
    hitsDF <- initiateHitsDF(hitsWithin, grangesReference, stranded)
    ## create single exon read class by using the minimum end
    ## and maximum start of all overlapping exons (identical to
    ## minimum equivalent class)
    hitsDF <- hitsDF %>% select(queryHits, chr, start, end, strand) %>%
        group_by(queryHits, chr, strand) %>%
        summarise(start = max(start), end = min(end), .groups = "drop") %>%
        group_by(chr, start, end, strand) %>%
        mutate(readClassId = paste0("rc", confidenceType, ".", 
            cur_group_id())) %>% ungroup()
    readIds <- mcols(granges[hitsDF$queryHits])$id
    hitsDF <- hitsDF %>% select(chr, start, end, strand, readClassId) %>%
        group_by(readClassId) %>% mutate(readCount = n()) %>% 
        ungroup() %>% distinct() %>% 
        mutate(confidenceType = confidenceType, intronStarts = NA,
            intronEnds = NA) %>%
        dplyr::select(chr, start, end, strand, intronStarts, intronEnds, 
            confidenceType, readClassId, readCount)
    exByReadClassUnspliced <- GRanges(
        seqnames = hitsDF$chr,
        ranges = IRanges(start = hitsDF$start, end = hitsDF$end),
        strand = hitsDF$strand)
    exByReadClassUnspliced$exon_rank <- 1
    exByReadClassUnspliced$exon_endRank <- 1
    partitioning <- PartitioningByEnd(seq_along(exByReadClassUnspliced))
    exByReadClassUnspliced <- relist(exByReadClassUnspliced, partitioning)
    names(exByReadClassUnspliced) <- hitsDF$readClassId
    hitsDF <- dplyr::select(hitsDF, chr.rc = chr, strand.rc = strand,
        intronStarts, intronEnds, confidenceType, readCount)
    mcols(exByReadClassUnspliced) <- hitsDF
    return(list(exonsByReadClass = exByReadClassUnspliced, readIds = readIds))
}

#' initiate the hits dataframe
#' @param hitsWithin hitsWithin
#' @param grangesReference grangesReference
#' @param stranded stranded
#' @importFrom dplyr as_tibble 
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
            as.factor(strand(grangesReference)[subjectHits(hitsWithin)])
    }
    return(hitsDF)
}



  assignGeneIdsByReference <- function(grl, annotations) {
    # (1) assign gene Ids based on first intron match to annotations
    geneRanges <- reducedRangesByGenes(annotations)
    ov=findOverlaps(grl, geneRanges)
    geneIds <- rep(NA, length(grl))
    uniqueHits <- which(queryHits(ov) %in% which(countQueryHits(ov)==1))
    geneIds[queryHits(ov)[uniqueHits]] <- names(geneRanges)[subjectHits(ov)[uniqueHits]]
    
    ## next for non unique hits select one gene (maximum overlap? using expand ranges? 
    multiHits <- which(queryHits(ov) %in% which(countQueryHits(ov)>1))
    expandedRanges <- expandRangesList(ranges(grl[queryHits(ov)[multiHits]]),
                                       ranges(geneRanges[subjectHits(ov)[multiHits]]))
    rangeIntersect <- pintersect(expandedRanges, mcols(expandedRanges)$matchRng, resolve.empty = 'start.x')
    intersectById <- tapply(width(rangeIntersect), mcols(expandedRanges)$IdMap, sum)
    
    filteredMultiHits <- as_tibble(ov[multiHits]) %>% mutate(intersectWidth = intersectById) %>% 
      group_by(queryHits) %>% arrange(desc(intersectWidth)) %>% slice(1)
    geneIds[filteredMultiHits$queryHits] <- names(geneRanges)[filteredMultiHits$subjectHits]
    
    return(geneIds)
    }
  
  #for test purpose
  #grl <- IRangesList(IRanges(c(1,90),c(5,95)), IRanges(c(10,20),c(15,25)), IRanges(c(20,30,40), c(25,35,45)),IRanges(c(40,50),c(45,55)), IRanges(c(50,60),c(55,65)), IRanges(c(70,80),c(75,85)),IRanges(c(80,90),c(85,95)), IRanges(c(200,210),c(205,215)), IRanges(c(300),c(305)))
  
  assignGeneIdsNoReference <- function(grl) {
    
    newTxIds <- 1:length(grl)
    newGeneByNewTxId <- rep(NA, length(newTxIds)) # will be filled in
    newTxIdsByExon <- rep(newTxIds, times=elementNROWS(grl))
    
    grSetReduced <- reduce(unlist(grl), with.revmap=T)
    newExonId <- 1:length(grSetReduced)
    revmap=mcols(grSetReduced)$revmap
    newExonIdByMergedExon <- rep(newExonId, times=elementNROWS(revmap))
    newTxIdsByMergedExon <- newTxIdsByExon[unlist(revmap)]
    
    # (1) create initial gene, exon, transcript map
    exonTxMap <- tibble(newExonId=newExonIdByMergedExon, newTxId=newTxIdsByMergedExon) %>% distinct()
    geneExonMap <-tibble(newGeneId=newExonIdByMergedExon, newExonId=newExonIdByMergedExon) %>% distinct()
    geneTxMap <- tibble(newGeneId=newExonIdByMergedExon, newTxId=newTxIdsByMergedExon) %>% distinct()
    
    show(table(is.na(newGeneByNewTxId)))
    # (2) select genes with only unique transcripts (complete assignment)
    exonGeneMap_filter1 <-  geneTxMap %>% group_by(newTxId) %>% mutate(nTx=n()) %>% 
      group_by(newGeneId) %>% filter(sum(nTx)==n()) 
    
    newGeneByNewTxId[exonGeneMap_filter1$newTxId] <- exonGeneMap_filter1$newGeneId
    show(table(is.na(newGeneByNewTxId)))
    
    if(any(is.na(newGeneByNewTxId))) {
      # (3.1) second iteration: non assigned genes
      geneTxMap <- geneTxMap %>% filter(!(newGeneId %in% exonGeneMap_filter1$newGeneId))
      exonTxMap <- exonTxMap %>% filter(!(newTxId %in% exonGeneMap_filter1$newTxId))
      geneExonMap <- geneExonMap %>% filter(!(newGeneId %in% exonGeneMap_filter1$newGeneId))
      
      ## (1) select gene Ids (first exons, can't be merged)
      refGeneExonList <- left_join(exonTxMap, rename(exonTxMap, newExonId.merge=newExonId)) %>% filter(newExonId<newExonId.merge) %>% filter(!(newExonId %in% newExonId.merge)) %>% select(newExonId, newExonId.merge) %>% distinct() %>% left_join(geneExonMap)  %>%  select(newGeneId, newExonId.merge) %>% distinct()
      
      #extend to the right
      refGeneTxMap <- left_join(refGeneExonList,  exonTxMap, by=c('newExonId.merge'='newExonId')) %>% select(newGeneId, newTxId) %>% distinct()
      
      refGeneTxMap <- rbind(geneTxMap %>% filter(newGeneId %in% refGeneTxMap$newGeneId), refGeneTxMap) %>% distinct()
      refGenExonMap <- left_join(refGeneTxMap, exonTxMap) %>% select(newGeneId, newExonId) %>% distinct()
      length_tmp = 0
      
      while(length_tmp<nrow(refGenExonMap)) {
        show('extend exons loop')
        length_tmp=nrow(refGenExonMap)
        refGeneTxMap <- left_join(refGenExonMap,  exonTxMap) %>% select(newGeneId, newTxId) %>% distinct()
        refGenExonMap <- left_join(refGeneTxMap, exonTxMap) %>% select(newGeneId, newExonId) %>% distinct()
      }
      geneTxMap <- refGeneTxMap
      geneExonMap <- refGenExonMap
      
      # combined gene ids
      geneGeneMap <- left_join(geneTxMap, rename(geneTxMap, newGeneId.merge=newGeneId)) %>% filter(newGeneId<newGeneId.merge) %>%filter(!(newGeneId %in% newGeneId.merge)) %>% select(newGeneId, newGeneId.merge) %>% distinct()
      geneTxMap <- rbind(geneTxMap %>%filter(!(newGeneId %in% geneGeneMap$newGeneId.merge)), left_join(geneGeneMap, geneTxMap, by=c('newGeneId.merge'='newGeneId')) %>% select(newGeneId, newTxId)) %>% distinct()
      
      newGeneByNewTxId[geneTxMap$newTxId] <- geneTxMap$newGeneId
      show(table(is.na(newGeneByNewTxId)))
    }
    return(newGeneByNewTxId)
  }
  
    
  assignGeneIds <- function(grl, annotations) {
    geneIds <- assignGeneIdsByReference(grl, annotations) 
    newGeneSet <- which(is.na(geneIds))
    newGeneIds <- assignGeneIdsNoReference(grl[newGeneSet])
    geneIds[newGeneSet] <- newGeneIds
  return(geneIds)
  }
