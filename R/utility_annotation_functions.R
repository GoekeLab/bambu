

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

