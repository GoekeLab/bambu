#' Nanopore transcript abundance quantification
#' @title transcript_abundance_quantification
#' @param readClassDt A \code{data.table} with columns
#' @importFrom BiocParallel bpparam bplapply
#' @noRd
abundance_quantification <- function(readClassDt, ncore = 1,
    maxiter = 20000, conv = 10^(-2), minvalue = 10^(-8)) {
    gene_sidList <- unique(readClassDt$gene_sid)
    if (ncore == 1) {
        emResultsList <- lapply(as.list(gene_sidList),
            run_parallel,
            conv = conv,
            minvalue = minvalue, 
            maxiter = maxiter,
            readClassDt = readClassDt
        )
    } else {
        bpParameters <- bpparam()
        bpParameters$workers <- ncore
        emResultsList <- bplapply(as.list(gene_sidList),
            run_parallel,
            conv = conv,
            minvalue = minvalue, 
            maxiter = maxiter,
            readClassDt = readClassDt,
            BPPARAM = bpParameters
        )
    }
    estimates <- do.call("rbind", emResultsList)
    return(estimates)
}


#' function to run in parallel
#' @title run_parallel
#' @param g the serial id of gene
#' @inheritParams abundance_quantification
#' @importFrom methods is
#' @import data.table
#' @noRd
run_parallel <- function(g, conv, minvalue, maxiter, readClassDt) {
    tmp <- unique(readClassDt[gene_sid == g])
    multiMap <- unique(tmp[, .(read_class_sid, multi_align)], 
            by = NULL)[order(read_class_sid)]$multi_align
    n.obs <- unique(tmp[, .(read_class_sid, nobs)], 
            by = NULL)[order(read_class_sid)]$nobs
    K <- as.numeric(sum(n.obs))
    n.obs <- as.numeric(n.obs)/K
    aMatArray <- formatAmat(tmp, multiMap)
    lambda <- sqrt(mean(n.obs))#suggested by Jiang and Salzman
    out <- initialiseOutput(dimnames(aMatArray)[[1]], g, K, n.obs) 
    #The following steps clean up the rc to tx matrix
    #step1:removes transcripts without any read class support
    a_mat <- aMatArray[,,1]
    if (is(a_mat,"numeric")) a_mat <- t(a_mat)
    rids <- which(apply(t(t(a_mat) * n.obs * K),1,sum) != 0)
    #step2: removes read classes without transcript assignment after step1
    a_mat <- aMatArray[rids,,1]
    cids <- which(apply(t(a_mat),1,sum) != 0)
    n.obs <- n.obs[cids]
    aMatArrayNew <- array(NA,dim = c(length(rids), length(cids),4))
    aMatArrayNew <- aMatArray[rids,cids,, drop = FALSE]
    if (is(aMatArrayNew[,,1],"numeric")) {
        aMatArrayUpdated <- K*n.obs*aMatArrayNew
        out[rids, `:=`(counts = aMatArrayUpdated[1],
            FullLengthCounts = aMatArrayUpdated[2],
            PartialLengthCounts = aMatArrayUpdated[3],
            UniqueCounts = aMatArrayUpdated[4])]
    }else{
        est_output <- emWithL1(A = aMatArrayNew, Y = n.obs, K = K,
            lambda = lambda, maxiter = maxiter,
            minvalue = minvalue, conv = conv)
        out <- modifyQuantOut(est_output, rids, cids, out)
    }
    return(out)
}

#' This function generates the wide-format a matrix
#' @import data.table
#' @noRd
formatAmat <- function(tmp, multiMap){
    tmp_wide <- dcast(tmp[order(nobs)],  tx_ori + 
        fullTx ~ read_class_sid, value.var = "aval")
    tmp_wide <- tmp_wide[CJ(tx_ori = unique(tmp_wide$tx_ori),
        fullTx = c(TRUE,FALSE)), on = c("tx_ori", "fullTx")]
    tmp_wide[is.na(tmp_wide)] <- 0
    a_mat_array <- array(NA,dim = c(nrow(tmp_wide)/2,
        ncol(tmp_wide) - 2,4), 
        dimnames = list(unique(tmp$tx_ori), unique(tmp$read_class_sid), NULL))
    a_mat_array[,,2] <- 
        as.matrix(tmp_wide[which(fullTx)][,-seq_len(2),with = FALSE])
    a_mat_array[,,3] <- 
        as.matrix(tmp_wide[which(!fullTx)][,-seq_len(2),with = FALSE])
    a_mat_array[,,4] <- a_mat_array[,,1] <- a_mat_array[,,2] + a_mat_array[,,3]
    a_mat_array[, which(multiMap), 4] <- 0
    return(a_mat_array)
}

#' This function generates a_mat values for all transcripts 
#' @import data.table
#' @noRd
modifyAvaluewithDegradation_rate <- function(tmp, d_rate, d_mode){
        tmp[, multi_align := (length(unique(tx_sid)) > 1),
            by = list(read_class_sid, gene_sid)]
        if (!d_mode) {
            tmp[, aval := 1]
            return(tmp)
        }
        tmp[which(multi_align) , aval := ifelse(fullTx, 1 -
            sum(.SD[which(!fullTx)]$rc_width*d_rate/1000),
            rc_width*d_rate/1000), by = list(gene_sid,tx_ori)]
        if (d_rate == 0) {
            tmp[, par_status := all(!fullTx & multi_align),
            by = list(read_class_sid, gene_sid)]
            tmp[which(par_status), aval := 0.01]
        }
        tmp[, aval := pmax(pmin(aval,1),0)] #d_rate should be contained to 0-1
        tmp[multi_align & fullTx,
            aval := pmin(1,pmax(aval,rc_width*d_rate/1000))]
        tmp[which(!multi_align), aval := 1]
        return(tmp)
}



#' This function initialises the final estimates with default values
#' @import data.table
#' @noRd 
initialiseOutput <- function(matNames, g, K, n.obs){
    return(data.table(tx_sid = matNames,counts = 0,
                FullLengthCounts = 0, PartialLengthCounts = 0,
                UniqueCounts = 0, theta = 0, gene_sid = g, 
                ntotal = as.numeric(K)))
}


#' Calculate degradation rate based on equiRC read counts 
#' @import data.table
#' @noRd
calculateDegradationRate <- function(readClassDt){
    rcCount <- unique(readClassDt[, .(gene_sid,read_class_sid, nobs)])
    rcCountPar <-
        unique(readClassDt[which(!fullTx), .(gene_sid,read_class_sid, nobs)])
    geneCount <- unique(rcCount[, list(nobs = sum(nobs)), by = gene_sid])
    geneCountPar <- unique(rcCountPar[, list(dObs = sum(nobs)), by = gene_sid])
    txLength <- unique(readClassDt[, .(gene_sid, tx_ori, tx_len)])
    geneLength <- 
        unique(txLength[, list(gene_len = max(tx_len)), by = gene_sid])
    geneCountLength <- unique(geneLength[geneCount, on = "gene_sid"])
    geneCountLength <- unique(geneCountPar[geneCountLength, on = "gene_sid"])
    geneCountLength[, d_rate := dObs/nobs]
    if (length(which(geneCountLength$nobs >= 30 & 
        ((geneCountLength$nobs - geneCountLength$dObs) >= 5))) == 0) {
        message("There is not enough read count and full length coverage!
            Hence degradation rate is estimated using all data!")
    } else {
        geneCountLength <- geneCountLength[nobs >= 30 & ((nobs - dObs) >= 5)]
    }
    d_rate <- median(geneCountLength$d_rate * 1000/geneCountLength$gene_len,
        na.rm = TRUE)
    return(c(d_rate, nrow(geneCountLength)))
}


#' Modify default quant output using estimated outputs
#' @import data.table
#' @noRd
modifyQuantOut <- function(est_output, rids, cids, out){
    est_out <- est_output[["theta"]]
    out[rids, `:=`(theta = est_out[1,],
        counts = est_out[2,],
        FullLengthCounts = est_out[3,],
        PartialLengthCounts = est_out[4,],
        UniqueCounts = est_out[5,])]
    return(out)
}


#' This function converts transcript, gene, and read class names to simple
#' integers for more efficient computation
#' @import data.table
#' @noRd
simplifyNames <- function(readClassDt, txVec, geneVec,ori_txvec, readclassVec){
    readClassDt <- as.data.table(readClassDt)
    readClassDt[, gene_sid := match(gene_id, geneVec)]
    readClassDt[, tx_sid := match(tx_id, txVec)]
    readClassDt[, tx_ori := match(gsub("Start","",tx_id),ori_txvec)]
    readClassDt[, read_class_sid := match(read_class_id, readclassVec)]
    readClassDt[, fullTx := grepl("Start",tx_id)]
    readClassDt[, `:=`(tx_id = NULL, gene_id = NULL, read_class_id = NULL)]
    return(readClassDt)
}

#' This function converts transcript and gene ids back to transcript and gene 
#' names 
#' @import data.table
#' @noRd
formatOutput <- function(theta_est, ori_txvec, geneVec){
    theta_est[, `:=`(tx_name = ori_txvec[as.numeric(tx_sid)],
        gene_name = geneVec[gene_sid])]
    theta_est[, `:=`(tx_sid = NULL, gene_sid = NULL)]
    theta_est <- theta_est[, .(tx_name, counts,FullLengthCounts,
        PartialLengthCounts, UniqueCounts, theta)]
    totalCount <- sum(theta_est$counts)
    theta_est[, `:=`(CPM = counts / totalCount * (10^6))]
    return(theta_est)
}


#' Remove duplicate transcript counts originated from multiple genes
#' @import data.table
#' @noRd
removeDuplicates <- function(counts){
    counts_final <- unique(counts[, list(counts = sum(counts),
                            FullLengthCounts = sum(FullLengthCounts),
                            PartialLengthCounts = sum(PartialLengthCounts),
                            UniqueCounts = sum(UniqueCounts),
                            theta = sum(theta),
                            CPM = sum(CPM)), by = tx_name],by = NULL)
    return(counts_final)
}

#' @import data.table
#' @noRd
removeUnObservedGenes <- function(readClassDt){
    uoGenes <- unique(readClassDt[,.I[sum(nobs) == 0], by = gene_sid]$gene_sid)
    if (length(uoGenes) > 0) {
        uo_txGeneDt <- 
            unique(readClassDt[(gene_sid %in% uoGenes),.(tx_ori,gene_sid)])
        readClassDt <- readClassDt[!(gene_sid %in% uoGenes)]
        outList <- data.table(tx_sid = uo_txGeneDt$tx_ori,
            counts = 0, 
            FullLengthCounts = 0,
            PartialLengthCounts = 0,
            UniqueCounts = 0,
            theta = 0,
            gene_sid = uo_txGeneDt$gene_sid, ntotal = 0)
    }else{
        outList <- NULL
    }
    return(list(readClassDt, outList))
}



#' This function aims to aggregate RC to equiRCs with the consideration of full
#' or partial alignment status and add empty read class depends on the minimal 
#' eqClass from annotations
#' @import data.table
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
#' @import data.table
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
#' @import data.table
#' @noRd
getUniCountPerEquiRC <- function(distTable){
    eqClassTable <- 
        unique(distTable[,.(read_class_id,readClassId, readCount, GENEID,
        firstExonWidth,totalWidth, equal)], by = NULL)
    eqClassTable[, `:=`(rc_width = ifelse(all(!equal), max(totalWidth), 
        max(firstExonWidth))), 
        by = list(read_class_id, GENEID)]
    eqClassTable <- unique(eqClassTable[, .(read_class_id,readClassId, 
        readCount, GENEID, rc_width)], by = NULL)
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
#' I.e., consider the scenario where there are transcripts with only 
#' shared read class observed, if just using the observed shared read class, we 
#' probably will see each of the transcript being similarly distributed, but we 
#' know that of the three transcripts, two of them are longer and have a unique 
#' part, in that case, if we add two empty read class for these two transcripts,
#' then we would know that given no unique read support found for these two 
#' longer transcripts, then the read count from this shared read class shall 
#' be more assigned to this smaller transcript. 
#' From this angle of view, both the full and partial of the minimal eqRC should
#' be added
#' @import data.table
#' @noRd
addEmptyRC <- function(annotations, equiRCTable){
    rcAnno <- data.table(as.data.frame(mcols(annotations)))
    rcAnno_partial <- copy(rcAnno)
    setnames(rcAnno_partial, "eqClass","read_class_id")
    rcAnno[, `:=`(read_class_id = gsub(paste0(TXNAME,"$"),
        paste0(TXNAME, "Start"), eqClass), 
        tx_id = paste0(TXNAME, "Start")), by = TXNAME]
    rcAnno[, `:=`(read_class_id = gsub(paste0(TXNAME,"\\."),
        paste0(TXNAME, "Start\\."), read_class_id)), by = TXNAME]
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
#' @import data.table
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


#' @useDynLib bambu, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL


#' @noRd
.onUnload <- function(libpath) {
    library.dynam.unload("bambu", libpath)
}

