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

#' This function aims to aggregate RC to equiRCs with the consideration of full
#' or partial alignment status and add empty read class depends on the minimal 
#' eqClass from annotations
#' @noRd
genEquiRCs <- function(readClass, annotations){
    ## aggregate rc's based on their alignment to full or partial transcripts
    distTable <- splitReadClass(readClass)
    ## get total number of read count for each equi read class
    eqClassTable <- getUniCountPerEquiRC(distTable)
    ## merge distTable with eqClassTable 
    equiRCTable <- unique(unique(distTable[, .(tx_id, GENEID, 
        read_class_id)], by = NULL)[eqClassTable, on = c("GENEID",
        "read_class_id"), allow.cartesian = TRUE], by = NULL)
    equiRCTable_final <- addEmptyRC(annotations, equiRCTable)
    setnames(equiRCTable_final, "GENEID", "gene_id")
    return(equiRCTable_final)
}

#' This function formats the distance table obtained from readClass by checking 
#' the distance between readClass and transcripts to create equiValence read 
#' classes
#' @noRd
splitReadClass <- function(readClass){
    rcWidth <- data.table(readClassId = rownames(readClass),
        firstExonWidth = width(unlist(rowRanges(readClass))[unlist(rowRanges(
        readClass))$exon_rank == 1,]),
        totalWidth = sum(width(rowRanges(readClass))))
    distTable <- data.table(metadata(readClass)$distTable)[, .(readClassId, 
        annotationTxId, readCount, GENEID, dist,equal)]
    distTable <- rcWidth[distTable, on = "readClassId"]
    # filter out multiple geneIDs mapped to the same readClass using rowData(se)
    compatibleData <- as.data.table(as.data.frame(rowData(readClass)),
        keep.rownames = TRUE)
    setnames(compatibleData, old = c("rn", "geneId"),
        new = c("readClassId", "GENEID"))
    distTable <- distTable[compatibleData[ readClassId %in% 
        unique(distTable$readClassId), .(readClassId, GENEID)],
        on = c("readClassId", "GENEID")]
    distTable[, tx_id := paste0(annotationTxId, ifelse(equal,"Start",""))]
    distTable[, read_class_id := paste(sort(unique(tx_id)), collapse = "."),
        by = list(readClassId, GENEID)]
    return(distTable)
}

#' get the count and rc_width for each equiRC
#' @noRd
getUniCountPerEquiRC <- function(distTable){
    eqClassTable <- 
        unique(distTable[,.(read_class_id,readClassId, readCount, GENEID,
        firstExonWidth,totalWidth, equal)], by = NULL)
    eqClassTable[, `:=`(rc_width = ifelse(all(!equal), max(totalWidth), 
        max(firstExonWidth))), 
        by = list(read_class_id, GENEID)]
    eqClassTable <- unique(eqClassTable[, .(read_class_id,readClassId, readCount, 
        GENEID, rc_width)], by = NULL)
    eqClassTable[, nobs := sum(readCount), by = list(read_class_id, GENEID)]
    eqClassTable <- unique(eqClassTable[,.(read_class_id, GENEID, nobs,
        rc_width)],by = NULL)
    return(eqClassTable)
}


#' This function adds the empty equiRCs to observed equiRcs to avoid
#' over-estimation of transcripts with unique part but no unique support
#' (Question: in fact, is this still necessary given the unique read support 
#' and unique and partial read counts are estimated specifically for each 
#' transcript?)
#' I.e., consider the scenario where there are there transcripts with only 
#' shared read class observed, if just using the observed shared read class, we 
#' probably will see each of the transcript being similarly distributed, but we 
#' know that of the three transcripts, two of them are longer and have a unique 
#' part, in that case, if we add two empty read class for these two transcripts, 
#' then we would know that given no unique read support found for these two 
#' longer transcripts, then the read count from this shared read class shall 
#' be more assigned to this smaller transcript. 
#' From this angle of view, both the full and partial of the minimal eqRC should
#' be added
#' @noRd
addEmptyRC <- function(annotations, equiRCTable){
    rcAnno <- data.table(as.data.frame(mcols(annotations)))
    rcAnno_partial <- copy(rcAnno)
    setnames(rcAnno_partial, "eqClass","read_class_id")
    rcAnno[, `:=`(read_class_id = gsub(paste0(TXNAME,"$"),
        paste0(TXNAME, "Start"), eqClass), 
        tx_id = paste0(TXNAME, "Start")), by = TXNAME]
    rcAnno[, `:=`(read_class_id = gsub(paste0(TXNAME,"\\."),
        paste0(TXNAME, "Start\\."), eqClass)), by = TXNAME]
    setnames(rcAnno_partial, "TXNAME", "tx_id")
    rcAnnoDt <-
        rbind(unique(rcAnno[, .(tx_id, GENEID, read_class_id)], by = NULL),
        unique(rcAnno_partial[, .(tx_id, GENEID, read_class_id)], by = NULL))
    rcAnnoDt <- createMultimappingBaseOnEmptyRC(rcAnnoDt)
    rcAnnoDt[, minEquiRC := 1] # this empty identifies read class observed only
    equiRCTable_final <- merge(equiRCTable, rcAnnoDt, 
        on = c("tx_id","GENEID","read_class_id"), all = TRUE)
    ## for is.na(nobs) 
    equiRCTable_final[is.na(nobs) ,`:=`(nobs = 0, rc_width = 0)]
    return(equiRCTable_final)
}

#' This function changes the concatenating symbol
#' @noRd
changeSymbol <- function(eqClass, txVec, from_symbol, to_symbol){
    uni_charVec <- unique(substr(txVec,1,1))
     for (uni_char in uni_charVec) {
        eqClass = gsub(paste0("\\",from_symbol,uni_char,""),
        paste0("\\",to_symbol,uni_char,""), eqClass)
     }
    return(eqClass)
}
#' This function creates the multi-mapping of read_class_id to tx_id
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
    ## for those that with txNumber being equal to txNumberExpected
    rcDt[is.na(tx_id_new), tx_id_new := tx_id]
    rcDt[tx_id_new != tx_id, tx_id := tx_id_new]
    rcDt[, tx_id_new := NULL]
    if (!(from_symbol == "&"))
    rcDt$read_class_id <- changeSymbol(rcDt$read_class_id,
        rcDt$tx_id, from_symbol = to_symbol, to_symbol = from_symbol)
    return(rcDt)
}
##Model the expected ratio of reads at a certain fraction of transcript
##Note here single isoform genes will be used for calculation
##Or should we use uniquely aligned reads?
#' @noRd
calculateExpectedCoverageRatio <- function(readClass, annotations, txLength){
    ## calculate coverage for all, as this is pretty fast
    rcranges <- rowRanges(readClass)
    ## cut first and last exon
    unlisted_anno <- unlist(annotations)
    annotations_t <- unlisted_anno[unlisted_anno$exon_rank != 1 & (
        unlisted_anno$exon_endRank != 1)]
    annotations_tList <- split(annotations_t, names(annotations_t))
    tx_cvg_numeric <- as(GenomicFeatures::coverageByTranscript(rcranges[rep(
        rownames(readClass), times = as.integer(assays(readClass)$counts))],
        annotations_tList,#[unlist(unique(strand(annotations))) != "*"],
        ignore.strand = TRUE), "NumericList")
    strandInfo <- unlist(unique(strand(annotations)))
    annotation_dt <- data.table(as.data.frame(mcols(annotations)))
    annotation_dt[, nisoform := length(unique(TXNAME)), by = GENEID]
    annotation_dt$eqClassNew <- changeSymbol(annotation_dt$eqClass,
        annotation_dt$TXNAME, from_symbol = ".", to_symbol = "&")
    candidate_tx <- filterTxForCvg(readClass, annotation_dt[TXNAME %in% 
        names(annotations_tList)])
    tx_cvg_df <- lapply(as.list(candidate_tx),
        calCoverage, tx_cvg_numeric = tx_cvg_numeric,
        strandInfo = strandInfo)
    tx_cvg <- do.call("rbind",tx_cvg_df)
    txLength <- data.table(annotationTxId = names(annotations_tList),
        txLength = sum(width(annotations_tList)))
    tx_cvg <- merge(tx_cvg, txLength, by.x = "tx_id",
        by.y = "annotationTxId", all.x = TRUE)
    tx_cvg[, rel_count := ave_bin_count/max(ave_bin_count), by = tx_id]
    tx_cvg[, adj_pos := txLength*pos_bin/100/1000]#per kb
    tx_cvg[pos_bin > 5 & (pos_bin < 96), 
        coef := lm(rel_count~adj_pos)$coef[2], by = tx_id]
    ## this rate is normalized by length and read 
    d_rate <- median(abs(unique(tx_cvg[!is.na(coef),.(tx_id, coef)])$coef))
    return(d_rate)
}


#' @noRd
filterTxForCvg <- function(readClass, annotation_dt){
    # For transcript with at least 5 full length read count
    # use only single isoform gene and only use a sample of such genes
    distTable <- data.table(as.data.frame(metadata(readClass)$distTable))
    distTable[, `:=`(totalN = sum(readCount),
        equalN = sum(equal*readCount)), by = annotationTxId]
    candidate_tx <- annotation_dt[nisoform == 1 & (!grepl("&",
        eqClassNew)) & (!grepl("unspliced", newTxClass))]$TXNAME
    candidate_tx <- intersect(candidate_tx,
        unique(distTable[equalN >= 5 & (totalN >= 30)]$annotationTxId))
    return(candidate_tx)
}

#' @noRd
calCoverage <- function(x, tx_cvg_numeric, strandInfo){
    dt <- data.table(count = tx_cvg_numeric[[x]],
        s = as.character(strandInfo[x]))
    dt[, `:=`(pos = ifelse(s == "-", rev(seq_along(dt$s)), seq_along(dt$s)))]
    dt[, pos_bin := ceiling(pos/max(pos)*100)]
    dt[, ave_bin_count := mean(count), by = pos_bin]
    dt[, tx_id := x]
    dt <- unique(dt[,.(pos_bin, ave_bin_count, tx_id)], by = NULL)
    return(dt)
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
