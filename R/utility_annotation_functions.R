


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
    rcAnno[, `:=`(read_class_id = gsub(TXNAME,paste0(TXNAME, "Start"),
        eqClass), tx_id = paste0(TXNAME, "Start")), by = TXNAME]
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
