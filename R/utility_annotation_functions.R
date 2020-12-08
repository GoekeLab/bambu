#' Prepare annotation granges object from GTF file
#' @title Prepare annotation granges object from GTF file into a 
#' GRangesList object
#' @param file a GTF file
#' @return A \code{\link{GRangesList}} object
#' @details Unlike \code{\link{readFromGTF}}, this function finds out the
#' equivalence classes between the transcripts,
#' with \code{\link{mcols}} data having three columns:
#' \itemize{
#'   \item TXNAME specifying prefix for new gene Ids (genePrefix.number),
#'   defaults to empty
#'   \item GENEID indicating whether filter to remove read classes
#'   which are a subset of known transcripts(), defaults to TRUE
#'   \item eqClass specifying minimun read count to consider a read class
#'    valid in a sample, defaults to 2
#'   }
#' @noRd
prepareAnnotationsFromGTF <- function(file) {
    if (missing(file)) {
        stop("A GTF file is required.")
    } else {
        data <- utils::read.delim(file, header = FALSE, comment.char = "#")
        colnames(data) <- c(
            "seqname", "source", "type", "start", "end",
            "score", "strand", "frame", "attribute")
        data <- data[data$type == "exon", ]
        data$strand[data$strand == "."] <- "*"
        data$GENEID <- gsub("gene_id (.*?);.*", "\\1", data$attribute)
        data$TXNAME <- gsub(".*transcript_id (.*?);.*", "\\1", data$attribute)
        geneData <- unique(data[, c("TXNAME", "GENEID")])
        grlist <- GenomicRanges::makeGRangesListFromDataFrame(
            data[, c("seqname", "start", "end", "strand", "TXNAME")],
            split.field = "TXNAME", keep.extra.columns = TRUE)
        grlist <- grlist[IRanges::order(start(grlist))]
        unlistedExons <- unlist(grlist, use.names = FALSE)
        partitioning <- PartitioningByEnd(cumsum(elementNROWS(grlist)),
            names = NULL)
        txIdForReorder <- togroup(PartitioningByWidth(grlist))
        exon_rank <- lapply(elementNROWS(grlist), seq, from = 1)
        exon_rank[which(unlist(unique(strand(grlist))) == "-")] <- lapply(
            exon_rank[which(unlist(unique(strand(grlist))) == "-")], rev
        ) # * assumes positive for exon ranking
        names(exon_rank) <- NULL
        unlistedExons$exon_rank <- unlist(exon_rank)
        unlistedExons <- unlistedExons[order(txIdForReorder,
            unlistedExons$exon_rank)]
        # exonsByTx is always sorted by exon rank, not by strand,
        # make sure that this is the case here
        unlistedExons$exon_endRank <- unlist(lapply(elementNROWS(grlist),
            seq, to = 1), use.names = FALSE)
        unlistedExons <- unlistedExons[order(txIdForReorder,
            start(unlistedExons))]
        mcols(unlistedExons) <- mcols(unlistedExons)[, c("exon_rank",
            "exon_endRank")]
        grlist <- relist(unlistedExons, partitioning)
        # sort the grlist by start position, ranked by exon number
        minEqClasses <- getMinimumEqClassByTx(grlist)
        mcols(grlist) <- DataFrame(geneData[(match(names(grlist),
            geneData$TXNAME)), ])
        mcols(grlist)$eqClass <- minEqClasses$eqClass[match(
            names(grlist), minEqClasses$queryTxId)]
    }
    return(grlist)
}

#' Get minimum equivalent class by Transcript
#' @param exonsByTranscripts exonsByTranscripts
#' @noRd
getMinimumEqClassByTx <- function(exonsByTranscripts) {
    exByTxAnnotated_singleBpStartEnd <-
        cutStartEndFromGrangesList(exonsByTranscripts)
    # estimate overlap only based on junctions
    spliceOverlaps <- findSpliceOverlapsQuick(
        exByTxAnnotated_singleBpStartEnd,
        exByTxAnnotated_singleBpStartEnd
    )
    ## identify transcripts compatible with other (subsets by splice sites)
    spliceOverlaps <- spliceOverlaps[mcols(spliceOverlaps)$compatible == TRUE, ]
    ## select splicing compatible transcript matches

    queryTxId <-
        names(exByTxAnnotated_singleBpStartEnd)[queryHits(spliceOverlaps)]
    subjectTxId <-
        names(exByTxAnnotated_singleBpStartEnd)[subjectHits(spliceOverlaps)]
    subjectTxId <- subjectTxId[order(queryTxId, subjectTxId)]
    queryTxId <- sort(queryTxId)
    eqClass <- unstrsplit(splitAsList(subjectTxId, queryTxId), sep = ".")

    return(tibble(queryTxId = names(eqClass), eqClass = unname(eqClass)))
}

#' calculate distance between first and last exon matches
#' @param candidateList candidateList
#' @param filteredOverlapList filteredOverlapList
#' @noRd
includeOverlapReadClass <- function(candidateList, filteredOverlapList) {
    temp <- left_join(candidateList, filteredOverlapList,
        by = c("subjectHits" = "queryHits")) %>%
        group_by(queryHits) %>%
        filter(!subjectHits.y %in% subjectHits, !is.na(subjectHits.y)) %>%
        ungroup() %>%
        dplyr::select(queryHits, subjectHits.y) %>%
        distinct() %>%
        dplyr::rename(subjectHits = subjectHits.y)
    return(temp)
}
#' Assign New Gene with Gene Ids
#' @param exByTx exByTx
#' @param prefix prefix, defaults to empty
#' @param minoverlap defaults to 5
#' @param ignore.strand defaults to FALSE
#' @noRd
assignNewGeneIds <- function(exByTx, prefix = "", minoverlap = 5,
    ignore.strand = FALSE) {
    if (is.null(names(exByTx))) names(exByTx) <- seq_along(exByTx)
    
    exonSelfOverlaps <- findOverlaps(exByTx, exByTx, select = "all",
        minoverlap = minoverlap, ignore.strand = ignore.strand)
    hitObject <- as_tibble(exonSelfOverlaps) %>% arrange(queryHits, subjectHits)
    candidateList <- hitObject %>% group_by(queryHits) %>%
        filter(queryHits <= min(subjectHits), queryHits != subjectHits) %>%
        ungroup()
    filteredOverlapList <- hitObject %>% filter(queryHits < subjectHits)
    rm(list = c("exonSelfOverlaps", "hitObject"))
    gc(verbose = FALSE)
    length_tmp <- 1
    # loop to include overlapping read classes which are not in order
    while (nrow(candidateList) > length_tmp) {
        length_tmp <- nrow(candidateList)
        temp <- includeOverlapReadClass(candidateList, filteredOverlapList)
        candidateList <- rbind(temp, candidateList)
        while (nrow(temp)) { ## annotate transcripts by new gene id
            temp <- includeOverlapReadClass(candidateList, filteredOverlapList)
            candidateList <- rbind(temp, candidateList)} ## second loop
        tst <- candidateList %>% group_by(subjectHits) %>%
            mutate(subjectCount = n()) %>% group_by(queryHits) %>%
            filter(max(subjectCount) > 1) %>% ungroup()
        temp2 <- inner_join(tst, tst, by = c("subjectHits" = "subjectHits")) %>%
            filter(queryHits.x != queryHits.y) %>%
            mutate(
                queryHits = if_else(queryHits.x > queryHits.y,
                    queryHits.y, queryHits.x),
                subjectHits = if_else(queryHits.x > queryHits.y,
                    queryHits.x, queryHits.y)) %>%
            dplyr::select(queryHits, subjectHits) %>%
            distinct()
        candidateList <- distinct(rbind(temp2, candidateList))
    }
    candidateList <- candidateList %>% filter(!queryHits %in% subjectHits) %>%
        arrange(queryHits, subjectHits)
    idToAdd <-
        (which(!(seq_along(exByTx) %in% unique(candidateList$subjectHits))))
    candidateList <- rbind(candidateList, tibble(
        queryHits = idToAdd, subjectHits = idToAdd
    )) %>% arrange(queryHits, subjectHits) %>%
        mutate(geneId = paste("gene", prefix, ".", queryHits, sep = "")) %>%
        dplyr::select(subjectHits, geneId)
    candidateList$readClassId <- names(exByTx)[candidateList$subjectHits]
    candidateList <- dplyr::select(candidateList, readClassId, geneId)
    return(candidateList)
}

#' generate filtered annotation table
#' @param spliceOverlaps an output from  findSpliceOverlapsByDist()
#' @param primarySecondaryDist default 5
#' @param primarySecondaryDistStartEnd default 5
#' @param exByTx default NULL
#' @param setTMP default NULL
#' @param DistCalculated default FALSE
#' @noRd
genFilteredAnTable <- function(spliceOverlaps, primarySecondaryDist = 5,
    primarySecondaryDistStartEnd = 5, exByTx = NULL, setTMP = NULL,
    DistCalculated = FALSE) {
    ## initiate the table
    if (isFALSE(DistCalculated)) {
        txToAnTable <- as_tibble(spliceOverlaps) %>% group_by(queryHits) %>%
            mutate(dist = uniqueLengthQuery + uniqueLengthSubject) %>%
            mutate(txNumber = n())
    } else {
        txToAnTable <- as_tibble(spliceOverlaps) %>% group_by(queryHits) %>%
            mutate(dist = uniqueLengthQuery + uniqueLengthSubject +
                uniqueStartLengthQuery + uniqueEndLengthQuery) %>%
            mutate(txNumber = n())
    }
    ## change query hits for step 2 and 3
    if (!is.null(exByTx)) {
        txToAnTable$queryHits <-
            (seq_along(exByTx))[-setTMP][txToAnTable$queryHits]
    }
    ## todo: check filters, what happens to reads with only start and end match?
    if (isFALSE(DistCalculated)) {
        txToAnTableFiltered <- txToAnTable %>%
            group_by(queryHits) %>%
            arrange(queryHits, dist) %>%
            filter(dist <= (min(dist) + primarySecondaryDist)) %>%
            filter(queryElementsOutsideMaxDist + 
                subjectElementsOutsideMaxDist == 
                min(queryElementsOutsideMaxDist +
                subjectElementsOutsideMaxDist)) %>% 
                filter((uniqueStartLengthQuery <= primarySecondaryDistStartEnd &
                uniqueEndLengthQuery <= primarySecondaryDistStartEnd) ==
                max(uniqueStartLengthQuery <=
                primarySecondaryDistStartEnd & uniqueEndLengthQuery <=
                primarySecondaryDistStartEnd)) %>%
            mutate(txNumberFiltered = n())
    } else {
        txToAnTableFiltered <- txToAnTable %>%
            group_by(queryHits) %>%
            arrange(queryHits, dist) %>%
            filter(dist <= (min(dist) + primarySecondaryDist)) %>%
            mutate(txNumberFiltered = n())
    }
    return(txToAnTableFiltered)
}
#' Calculate distance from read class to annotation
#' @param exByTx exByTx
#' @param exByTxRef exByTxRef
#' @param maxDist defaults to 35
#' @param primarySecondaryDist defaults to 5
#' @param ignore.strand defaults to FALSE
#' @noRd
calculateDistToAnnotation <- function(exByTx, exByTxRef, maxDist = 35,
    primarySecondaryDist = 5, primarySecondaryDistStartEnd = 5,
    ignore.strand = FALSE) {
    # (1)  find overlaps of read classes with annotated transcripts,
    spliceOverlaps <- findSpliceOverlapsByDist(exByTx, exByTxRef,
        maxDist = maxDist, firstLastSeparate = TRUE,
        dropRangesByMinLength = TRUE, cutStartEnd = TRUE,
        ignore.strand = ignore.strand)
    txToAnTableFiltered <- genFilteredAnTable(spliceOverlaps,
        primarySecondaryDist, DistCalculated = FALSE)
    # (2) calculate splice overlap for any not in the list (new exon >= 35bp)
    setTMP <- unique(txToAnTableFiltered$queryHits)
    spliceOverlaps_rest <- findSpliceOverlapsByDist(exByTx[-setTMP],
        exByTxRef, maxDist = 0, type = "any", firstLastSeparate = TRUE,
        dropRangesByMinLength = FALSE, cutStartEnd = TRUE,
        ignore.strand = ignore.strand)
    txToAnTableRest <-
        genFilteredAnTable(spliceOverlaps_rest, primarySecondaryDist,
        exByTx = exByTx, setTMP = setTMP, DistCalculated = FALSE)
    # (3) find overlaps for remaining reads 
    setTMPRest <- unique(c(txToAnTableRest$queryHits, setTMP))
    txToAnTableRestStartEnd <- NULL
    if (length(exByTx[-setTMPRest])) {
        spliceOverlaps_restStartEnd <-
            findSpliceOverlapsByDist(exByTx[-setTMPRest], exByTxRef,
            maxDist = 0, type = "any", firstLastSeparate = TRUE,
            dropRangesByMinLength = FALSE,
            cutStartEnd = FALSE, ignore.strand = ignore.strand)
        if (length(spliceOverlaps_restStartEnd)) {
            txToAnTableRestStartEnd <-
                genFilteredAnTable(spliceOverlaps_restStartEnd,
                primarySecondaryDist, exByTx = exByTx,
                setTMP = setTMPRest, DistCalculated = TRUE)
        }
    }
    txToAnTableFiltered <- rbind( txToAnTableFiltered,
        txToAnTableRest, txToAnTableRestStartEnd ) %>% ungroup()
    txToAnTableFiltered$readClassId <-
        names(exByTx)[txToAnTableFiltered$queryHits]
    txToAnTableFiltered$annotationTxId <-
        names(exByTxRef)[txToAnTableFiltered$subjectHits]
    return(txToAnTableFiltered)
}

#' Get Empty Read Class From SE
#' @param se summarizedExperiment
#' @param annotationGrangesList defaults to NULL
#' @param max.distScore defaults to 5 
#' @noRd
getEmptyClassFromSE <- function(se = se, annotationGrangesList = NULL) {
    txLength <- data.table(annotationTxId = names(annotationGrangesList),
        txLength = sum(width(annotationGrangesList)))
    distTable <- data.table(metadata(se)$distTable)[, .(readClassId, 
        annotationTxId, readCount, GENEID, dist,equal)]
    distTable <- txLength[distTable, on = "annotationTxId"]
    # filter out multiple geneIDs mapped to the same readClass using rowData(se)
    compatibleData <- as.data.table(as.data.frame(rowData(se)),
        keep.rownames = TRUE)
    setnames(compatibleData, old = c("rn", "geneId"),
        new = c("readClassId", "GENEID"))
    distTable <- distTable[compatibleData[ readClassId %in% 
        unique(distTable$readClassId), .(readClassId, GENEID)],
        on = c("readClassId", "GENEID")]
    distTable[, eqClass := paste(sort(unique(annotationTxId)), collapse = "."),
        by = list(readClassId, GENEID)]
    distTable[, `:=`(equalN = sum(equal*readCount)), 
        by = list(eqClass,annotationTxId)]
    rcTable <- unique(distTable[, .(readClassId, GENEID,
        eqClass, readCount)])
    rcTable[, `:=`(eqClassReadCount = sum(readCount)),
        by = list(eqClass, GENEID)]
    rcTable <- unique(rcTable[, .(eqClass, eqClassReadCount, GENEID)])
    eqClassCountTable <- unique(distTable[, .(annotationTxId, GENEID,
        eqClass, equalN,txLength)][rcTable, 
        on = c("GENEID", "eqClass")])
    setnames(eqClassCountTable, c("annotationTxId"), c("TXNAME"))
    eqClassTable <- as.data.table(mcols(annotationGrangesList)[,
        c("GENEID", "eqClass", "TXNAME")])
    eqClassCountTable <- unique(merge(eqClassCountTable, eqClassTable,
        all = TRUE, on = c("GENEID", "eqClass", "TXNAME")))
    #  new isoforms from eqClassCountTable should be kept
    eqClassCountTable[is.na(eqClassReadCount), eqClassReadCount := 0]
    eqClassCountTable[is.na(equalN), equalN := 0]
    eqClassCountTable[, sum_nobs := sum(eqClassReadCount),
        by = list(GENEID, TXNAME)]
    eqClassCountTable <- unique(eqClassCountTable[sum_nobs > 0, .(
        GENEID,eqClass, eqClassReadCount, TXNAME, equalN,txLength)])
        setnames(eqClassCountTable, old = c("TXNAME", "GENEID", "eqClass",
        "eqClassReadCount", "equalN","txLength"), new = c("tx_id", "gene_id",
        "read_class_id", "nobs","equalN","tx_len"))
    return(eqClassCountTable)
}

#' @noRd 
createMultimappingBaseOnEmptyRC <- function(rcDt, 
    from_symbol = ".", to_symbol = "&"){
    if (!(from_symbol == "&"))
        rcDt$read_class_id <- changeSymbol(rcDt$read_class_id,
        rcDt$tx_id, from_symbol, to_symbol)
    rcDtNew <- rcDt[, list(txNumber = length(unique(tx_id)),
        txNumberExpected = length(unlist(strsplit(read_class_id, "\\&")))),
        by = read_class_id]
    ## For those with txNumber less than txNumberExpected 
    rcDtNew_remap <- rcDtNew[txNumber != txNumberExpected, 
        list(tx_id_new = unlist(strsplit(read_class_id, "\\&"))),
        by = read_class_id]
    rcDt <- rcDtNew_remap[rcDt, on = "read_class_id"]
    rcDt[is.na(tx_id_new), tx_id_new := tx_id]
    rcDt[tx_id_new != tx_id, tx_id := tx_id_new]
    rcDt[, tx_id_new := NULL]
    if (!(from_symbol == "&"))
    rcDt$read_class_id <- changeSymbol(rcDt$read_class_id,
        rcDt$tx_id, from_symbol = to_symbol, to_symbol = from_symbol)
    return(rcDt)
}

#' @noRd
changeSymbol <- function(eqClass, txVec, from_symbol, to_symbol){
    uni_charVec <- unique(substr(txVec,1,1))
     for (uni_char in uni_charVec) {
        eqClass = gsub(paste0("\\",from_symbol,uni_char,""),
        paste0("\\",to_symbol,uni_char,""), eqClass)
     }
    return(eqClass)
}

#' @noRd
getSubsetAnno <- function(rcAnnotations){
    ## for subset transcripts 
    subset_anno <- unique(rcAnnotations[which(subset_status), 
        list(tx_id = unlist(strsplit(eqClass, "\\&"))),
        by = list(GENEID, TXNAME, eqClass)],by = NULL)
    subset_anno <- subset_anno[TXNAME != tx_id]
    setnames(subset_anno, c("eqClass","TXNAME"), 
        c("read_class_id","subset_txid"))
    return(subset_anno)
}

#' @noRd
getSubsetReadClassTable <- function(subset_anno, readClassDt, rcAnnotations){
    readClassDtMatched <- unique(subset_anno[,.(tx_id, subset_txid)], 
        by = NULL)[readClassDt, on = c("tx_id"), allow.cartesian = TRUE]
    readClassDtMatched <- readClassDtMatched[!is.na(subset_txid)]
    readClassDtMatched[, exist := grepl(subset_txid, read_class_id),
        by = list(subset_txid,read_class_id)]
    readClassDtMatched <- readClassDtMatched[which(exist)]
    readClassDtMatched[, total_txNumber := 
        length(unlist(strsplit(read_class_id, "\\&"))), by = read_class_id]
    subset_rcDt <- readClassDtMatched[, list(subset_txNumber = .N), 
        by = list(read_class_id,subset_txid, total_txNumber)] 
    
    ## only keep those with total_txNumber = subset_txNumber + 1 
    subset_rcDt <- subset_rcDt[total_txNumber == subset_txNumber + 1]
    subset_rcDt[, `:=`(total_txNumber = NULL, subset_txNumber = NULL)]
    subset_rcDt <- rbind(subset_rcDt, 
        data.table(read_class_id = rcAnnotations[which(!subset_status)]$eqClass,
        subset_txid = rcAnnotations[which(!subset_status)]$eqClass))
    subset_rcDt[, read_class_id_new := gsub(subset_txid, 
        paste0(subset_txid,"Start"), read_class_id), 
        by = list(read_class_id, subset_txid)]
    return(subset_rcDt)
}

#' @noRd
createFullEquiCounts <- function(subset_rcDt, readClassDt){
    # first merge subset_rcDt with readClassDt
    readClassDtNew <- subset_rcDt[readClassDt, on = "read_class_id"]
    readClassDtNew <- 
        readClassDtNew[!is.na(subset_txid) & (tx_id == subset_txid)]
    readClassDtNew[, `:=`(tx_id = paste0(tx_id, "Start"),
                          nobs = equalN,
                          read_class_id = read_class_id_new)]
    readClassDtNew[, `:=`(subset_txid = NULL,
                          read_class_id_new = NULL)]
    ## whenever edited, remember to create multimapping 
    readClassDtNew <- createMultimappingBaseOnEmptyRC(readClassDtNew, 
        from_symbol = "&")
    # secondly manual edit those are missed by a few bp difference in splicing 
    # sites
    readClassDtNewManual <- readClassDt[!(read_class_id %in% 
        unique(subset_rcDt$read_class_id))][equalN > 0]
    readClassDtNewManual[, `:=`(tx_id = paste0(tx_id, "Start"),
        nobs = equalN,
        read_class_id = gsub(tx_id, paste0(tx_id,"Start"), read_class_id)), 
        by = list(tx_id, read_class_id)]
    readClassDtNewManual <- 
        createMultimappingBaseOnEmptyRC(readClassDtNewManual, from_symbol = "&")
    readClassDtOld <- copy(readClassDt)
    readClassDtOld[, nobs := nobs - equalN]
    readClassDt_final <- do.call("rbind", list(readClassDtNew, readClassDtOld,
        readClassDtNewManual))
    return(readClassDt_final)
}

##============================================================================##
# Steps in modifying readClass
# 1. change concatenating symbol to & to avoid mis-spliting for novel transcript 
#    find the possible first characters in the tx_id
# 2. Modify read class that can be found with a minimal subset transcript 
#    with a full length version of the minimal subset transcript indicated 
#    by original transcript followed by Start (this set by right, should be able
#    to find all shared read class with a subset transcript matching and all
#    transcript with unique splicing )
# 3. Transcripts that are subset of shared read class, if they have a read class
#    for themseleves, modify by adding a full length version as well
# 4. For shared read class that do not have a subset transcript that is 
#    fully compatible with the splicing patterns, do not modify 
# 5. change concatenating symbol back to .

#' Modify readClass to create full length read class
#' @param readClassDt output from \code{getEmptyClassFromSE}
#' @param annotationGrangesList inherits from \code{getEmptyClassFromSE}
#' @noRd
modifyReadClassWtFullLengthTranscript <- function(readClassDt,
    annotationGrangesList){
    readClassDt <- createMultimappingBaseOnEmptyRC(readClassDt)
    rcAnnotations <- data.table(as.data.frame(mcols(annotationGrangesList)))
    rcAnnotations$eqClass <- changeSymbol(eqClass = rcAnnotations$eqClass, 
        txVec = rcAnnotations$TXNAME, from_symbol = ".", to_symbol = "&")
    rcAnnotations[, subset_status := grepl("\\&",eqClass)]
    subset_anno <- getSubsetAnno(rcAnnotations)
    readClassDt$read_class_id <- changeSymbol(readClassDt$read_class_id,
        readClassDt$tx_id, from_symbol = ".", to_symbol = "&")
    subset_rcDt <-
        getSubsetReadClassTable(subset_anno, readClassDt, rcAnnotations)
    readClassDt <- createFullEquiCounts(subset_rcDt, readClassDt)
    readClassDt$read_class_id <- changeSymbol(readClassDt$read_class_id,
        readClassDt$tx_id, from_symbol = "&", to_symbol = ".")
    readClassDt[, `:=`(equalN = NULL)]
   return(readClassDt)
}

##Model the expected ratio of reads at a certain fraction of transcript
##Note here single isoform genes will be used for calculation
##Or should we use uniquely aligned reads?
# calculateExpectedCoverageRatio <- function(rcranges, annotatons){
#     ## calculate coverage for all, as this is pretty fast
#     tx_cvg_numeric <- as(GenomicFeatures::coverageByTranscript(rcranges, 
#         annotations,#[unlist(unique(strand(annotations))) != "*"], 
#         ignore.strand = TRUE), "NumericList")
#     strandInfo <- unlist(unique(strand(annotations)))
#     ## use only single isoform gene and only use a sample of such genes
#     annotation_dt <- data.table(as.data.frame(mcols(annotations)))
#     annotation_dt[, nisoform := length(unique(TXNAME)), by = GENEID]
#     annotation_dt$eqClassNew <- changeSymbol(annotation_dt$eqClass, 
#         annotation_dt$TXNAME, from_symbol = ".", to_symbol = "&")
#     candidate_tx <- annotation_dt[nisoform == 1 & (!grepl("&",
#         eqClassNew)) & (!grepl("unspliced", newTxClass))]$TXNAME
#     candidate_tx <- intersect(candidate_tx, 
#         unique(metadata(readClass)$distTable$annotationTxId))
#     system.time(tx_cvg_df <- BiocParallel::bplapply(as.list(candidate_tx), 
#         calCoverage,
#         tx_cvg_numeric = tx_cvg_numeric, 
#         strandInfo = strandInfo,
#         BPPARAM = bpParameters))
#     tx_cvg <- do.call("rbind",tx_cvg_df)
#     # tx_cvg[, normAveCount := ave_bin_count/max(ave_bin_count), by = tx_id]
#     # tx_cvg_ave <- tx_cvg[!is.na(normAveCount), 
#     #     list(aveNorm = mean(normAveCount)),by = pos_bin]
#     tx_cvg_ave <- tx_cvg[!is.na(normAveCount),
#         list(aveCount= mean(ave_bin_count)),by = pos_bin]
#     ## 
#     ecdf_fun <- approxfun(tx_cvg_ave$pos_bin/100, 
#         tx_cvg_ave$aveCount/max(tx_cvg_ave$aveCount),
#         method = "linear", yleft = 0)
#     return(ecdf_fun)
# }


#' #' @noRd
#' calCoverage <- function(x, tx_cvg_numeric, strandInfo){
#'     dt <- data.table(count = tx_cvg_numeric[[x]],
#'         s = as.character(strandInfo[x]))
#'     dt[, `:=`(pos = ifelse(s == "t", rev(seq_along(dt$s)), seq_along(dt$s)))]
#'     dt[, pos_bin := ceiling(pos/max(pos)*100)]
#'     dt[, ave_bin_count := mean(count), by = pos_bin]
#'     dt[, tx_id := x]
#'     dt <- unique(dt[,.(pos_bin, ave_bin_count, tx_id)], by = NULL)
#'     return(dt)
#' }

## degradation estimation
#' @noRd
estDegradation_factor <- function(se, annotationGrangesList){
    txLength <- data.table(annotationTxId = names(annotationGrangesList),
        txLength = sum(width(annotationGrangesList)))
    distTable <- data.table(metadata(se)$distTable)[, .(readClassId, 
        annotationTxId, readCount, GENEID, dist,equal)]
    distTable <- txLength[distTable, on = "annotationTxId"]
    
    annotation_dt <- data.table(as.data.frame(mcols(annotations)))
    annotation_dt[, nisoform := length(unique(TXNAME)), by = GENEID]
    annotation_dt$eqClassNew <- changeSymbol(annotation_dt$eqClass,
        annotation_dt$TXNAME, from_symbol = ".", to_symbol = "&")
    candidate_tx <- annotation_dt[nisoform == 1 & (!grepl("&",
        eqClassNew)) & (!grepl("unspliced", newTxClass))]$TXNAME
    candidate_tx <- intersect(candidate_tx,
        unique(metadata(readClass)$distTable$annotationTxId))
    
    degradation_data <- unique(distTable[annotationTxId %in% 
        candidate_tx][,.(annotationTxId, txLength, 
        tx_count, full_count, GENEID)])
    
    degradation_rate <- mean((1 - degradation_data2[tx_count > 
        30]$full_count/degradation_data2[tx_count > 30]$tx_count)/ceiling(
        degradation_data2[tx_count > 30]$txLength/1000))
    
    return(degradation_rate)
}
#' From tx ranges to gene ranges
#' @noRd
txRangesToGeneRanges <- function(exByTx, TXNAMEGENEID_Map) {
    # rename names to geneIDs
    names(exByTx) <- as.data.table(TXNAMEGENEID_Map)[match(names(exByTx),
        TXNAME)]$GENEID

    # combine gene exon ranges and reduce overlapping ones
    unlistData <- unlist(exByTx, use.names = TRUE)
    orderUnlistData <- unlistData[order(names(unlistData))]
    orderUnlistData$exon_rank <- NULL
    orderUnlistData$exon_endRank <- NULL

    exByGene <- splitAsList(orderUnlistData, names(orderUnlistData))
    exByGene <- GenomicRanges::reduce(exByGene)

    # add exon_rank and endRank
    unlistData <- unlist(exByGene, use.names = FALSE)
    partitionDesign <- cumsum(elementNROWS(exByGene))
    partitioning <- PartitioningByEnd(partitionDesign, names = NULL)
    geneStrand <- as.character(strand(unlistData))[partitionDesign]
    exon_rank <- lapply(width((partitioning)), seq, from = 1)
    exon_rank[which(geneStrand == "-")] <-
        lapply(exon_rank[which(geneStrand == "-")], rev)
    # * assumes positive for exon ranking
    exon_endRank <- lapply(exon_rank, rev)
    unlistData$exon_rank <- unlist(exon_rank)
    unlistData$exon_endRank <- unlist(exon_endRank)
    exByGene <- relist(unlistData, partitioning)

    return(exByGene)
}
