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
    multiMap <- unique(tmp[, .(eqClassId, multi_align)], 
                       by = NULL)[order(eqClassId)]$multi_align
    n.obs <- unique(tmp[, .(eqClassId, nobs)], 
                    by = NULL)[order(eqClassId)]$nobs
    K <- as.numeric(sum(n.obs))
    n.obs <- as.numeric(n.obs)/K
    aMatArray <- formatAmat(tmp, multiMap)
    lambda <- sqrt(mean(n.obs))#suggested by Jiang and Salzman
    out <- initialiseOutput(dimnames(aMatArray)[[1]], g, K, n.obs) 
    #The following steps clean up the rc to tx matrix
    #step1:removes transcripts without any read class support
    a_mat <- aMatArray[,,1]
    # in the case, when a_mat dimension is of 1 need to transform 1 more time
    if (is(a_mat,"numeric")) a_mat <- t(a_mat)
    rids <- which(apply(t(t(a_mat) * n.obs * K),1,sum) != 0)
    #step2: removes read classes without transcript assignment after step1
    a_mat <- aMatArray[rids,,1]
    cids <- which(apply(t(a_mat),1,sum) != 0)
    if( is(a_mat, "vector")) cids <- which(a_mat != 0)
    n.obs <- n.obs[cids]
    aMatArrayNew <- array(NA,dim = c(length(rids), length(cids),3))
    aMatArrayNew <- aMatArray[rids,cids,, drop = FALSE]
    if (is(aMatArrayNew[,,1],"numeric")&(nrow(aMatArrayNew)==1)) {
        aMatArrayUpdated <- K*n.obs*aMatArrayNew
        out[rids, `:=`(counts = sum(aMatArrayUpdated[,,1]),
                       FullLengthCounts = sum(aMatArrayUpdated[,,2]),
                       UniqueCounts = sum(aMatArrayUpdated[,,3]))]
    }else{
        est_output <- emWithL1(A = aMatArrayNew, Y = n.obs, K = K,
                               lambda = lambda, maxiter = maxiter,
                               minvalue = minvalue, conv = conv)
        out <- modifyQuantOut(est_output, rids, cids, out)
    }
    end.ptm <- proc.time()
    return(out)
}

#' This function generates the wide-format a matrix
#' @import data.table
#' @noRd
formatAmat <- function(tmp, multiMap){
    tmp_wide <- dcast(tmp[order(nobs)],  txid + 
                          equal ~ eqClassId, value.var = "aval") # this step will automatically sort the table
    tmp_wide <- tmp_wide[CJ(txid = unique(tmp_wide$txid),
                            equal = c(TRUE,FALSE)), on = c("txid", "equal")]
    tmp_wide[is.na(tmp_wide)] <- 0
    a_mat_array <- array(NA,dim = c(nrow(tmp_wide)/2,
                                    ncol(tmp_wide) - 2,3), 
                         dimnames = list(sort(unique(tmp$txid)), sort(unique(tmp$eqClassId)), NULL))
    a_mat_array[,,2] <- 
        as.matrix(tmp_wide[which(equal)][,-seq_len(2),with = FALSE])
    a_mat_array[,,3] <- a_mat_array[,,1] <- a_mat_array[,,2] + as.matrix(tmp_wide[which(!equal)][,-seq_len(2),with = FALSE])
    a_mat_array[, which(multiMap), 3] <- 0
    return(a_mat_array)
}

#' This function generates a_mat values for all transcripts 
#' @import data.table
#' @noRd
modifyAvaluewithDegradation_rate <- function(tmp, d_rate, d_mode){
    tmp[, multi_align := (length(unique(txid)) > 1),
        by = list(eqClassId, gene_sid)]
    if (!d_mode) {
        tmp[, aval := 1]
        return(tmp)
    }
    tmp[which(multi_align) , aval := ifelse(equal, 1 -
        sum(.SD[which(!equal)]$rcWidth*d_rate/1000),
                                            rcWidth*d_rate/1000), by = list(gene_sid,txid)]
    if (d_rate == 0) {
        tmp[, par_status := all(!equal & multi_align),
            by = list(eqClassById, gene_sid)]
        tmp[which(par_status), aval := 0.01]
    }
    tmp[, aval := pmax(pmin(aval,1),0)] #d_rate should be contained to 0-1
    tmp[multi_align & equal,
        aval := pmin(1,pmax(aval,rcWidth*d_rate/1000))]
    tmp[which(!multi_align), aval := 1]
    return(tmp)
}



#' This function initialises the final estimates with default values
#' @import data.table
#' @noRd 
initialiseOutput <- function(matNames, g, K, n.obs){
    return(data.table(txid = sort(as.numeric(matNames)),counts = 0,
                      FullLengthCounts = 0,
                      UniqueCounts = 0, gene_sid = g, 
                      ntotal = as.numeric(K)))
}


#' Calculate degradation rate based on equiRC read counts 
#' @import data.table
#' @noRd
calculateDegradationRate <- function(readClassDt){
    rcCount <- unique(readClassDt[, .(gene_sid,eqClassId, nobs)])
    rcCountPar <-
        unique(readClassDt[which(!equal), .(gene_sid,eqClassId, nobs)])
    geneCount <- unique(rcCount[, list(nobs = sum(nobs)), by = gene_sid])
    geneCountPar <- unique(rcCountPar[, list(dObs = sum(nobs)), by = gene_sid])
    txLength <- unique(readClassDt[, .(gene_sid, txid, txlen)])
    geneLength <- 
        unique(txLength[, list(gene_len = max(txlen)), by = gene_sid])
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
    out[rids, `:=`(
        counts = est_out[1,],
        FullLengthCounts = est_out[2,],
        UniqueCounts = est_out[3,])]
    return(out)
}


#' This function converts transcript, gene, and read class names to simple
#' integers for more efficient computation
#' @import data.table
#' @noRd
simplifyNames <- function(readClassDt, geneVec){
    readClassDt <- as.data.table(readClassDt)
    readClassDt[, gene_sid := match(GENEID, geneVec)]
    readClassDt[, `:=`(GENEID = NULL, eqClassById = NULL)]
    return(readClassDt)
}

#' This function converts transcript and gene ids back to transcript and gene 
#' names 
#' @noRd
formatOutput <- function(theta_est, txVec, geneVec){
    theta_est[, `:=`(gene_name = geneVec[gene_sid])]
    theta_est[, `:=`(gene_sid = NULL)]
    theta_est <- theta_est[, .(txid, counts,FullLengthCounts,
UniqueCounts)]
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
                                         UniqueCounts = sum(UniqueCounts),
                                         CPM = sum(CPM)), by = txid],by = NULL)
    return(counts_final)
}

#' @import data.table
#' @noRd
removeUnObservedGenes <- function(readClassDt){
    uoGenes <- unique(readClassDt[,.I[sum(nobs) == 0], by = gene_sid]$gene_sid)
    if (length(uoGenes) > 0) {
        uo_txGeneDt <- 
            unique(readClassDt[(gene_sid %in% uoGenes),.(txid,gene_sid)])
        readClassDt <- readClassDt[!(gene_sid %in% uoGenes)]
        outList <- data.table(txid = uo_txGeneDt$txid,
                              counts = 0, 
                              FullLengthCounts = 0,
                              UniqueCounts = 0,
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
genEquiRCs <- function(readClassDist, annotations, verbose){
    start.ptm <- proc.time()
    ## aggregate rc's based on their alignment to full or partial transcripts
    distTable <- genEquiRCsBasedOnObservedReads(readClassDist)
    eqClassCount <- getUniCountPerEquiRC(distTable)
    eqClassTable <- addEmptyRC(eqClassCount, annotations)
    # create equiRC id 
    eqClassTable <- eqClassTable %>% 
        group_by(eqClassById) %>%
        mutate(eqClassId = cur_group_id())
    return(data.table(eqClassTable))
}

#' This function formats the distance table obtained from readClass by checking 
#' the distance between readClass and transcripts to create equiValence read 
#' classes
#' @import data.table
#' @noRd
genEquiRCsBasedOnObservedReads <- function(readClass){
    unlisted_rowranges <- unlist(rowRanges(readClass))
    rcWidth <- data.table(readClassId = rownames(readClass),
                          firstExonWidth =
                              width(unlisted_rowranges[unlisted_rowranges$exon_rank == 1,]),
                          totalWidth = sum(width(rowRanges(readClass))))
    distTable <- data.table(as.data.frame(metadata(readClass)$distTable))[!grepl("unidentified", annotationTxId), .(readClassId, 
                                                               annotationTxId, readCount, GENEID, dist,equal,txid)]
    distTable <- rcWidth[distTable, on = "readClassId"]
    # filter out multiple geneIDs mapped to the same readClass using rowData(se)
    compatibleData <- as.data.table(as.data.frame(rowData(readClass)),
                                    keep.rownames = TRUE)
    setnames(compatibleData, old = c("rn", "geneId"),
             new = c("readClassId", "GENEID"))
    distTable <- distTable[compatibleData[ readClassId %in% 
                                               unique(distTable$readClassId), .(readClassId, GENEID)],
                           on = c("readClassId", "GENEID")]
    #here, each transcript should be assigned to one gene only based on isore.estimateDistanceToAnnotation function
    ##this step is very slow, consider to use integers instead of tx_ids
    distTable[order(readClassId,GENEID,txid), 
              eqClassById :=.(list(sort(unique(ifelse(equal, -1*txid,txid))))),
              by = list(readClassId, GENEID)]
    return(distTable)
}

#' get the count and rc_width for each equiRC
#' @import data.table
#' @noRd
getUniCountPerEquiRC <- function(distTable){
    eqClassCount <- distTable %>% 
        group_by(eqClassById) %>%
        mutate(anyEqual = any(equal)) %>%
        select(eqClassById, firstExonWidth,totalWidth, readCount,GENEID,anyEqual) %>% #eqClassByIdTemp,
        distinct() %>%
        mutate(nobs = sum(readCount),
               rcWidth = ifelse(anyEqual, max(totalWidth), 
                               max(firstExonWidth))) %>%
        select(eqClassById,GENEID,nobs,rcWidth) %>% #eqClassByIdTemp,
        ungroup()  %>%
        distinct()
    return(eqClassCount)
}


# add minimal equiRC 
#' @import data.table
#' @noRd
addEmptyRC <- function(eqClassCount, annotations){
    minEquiRC <- as.data.frame(mcols(annotations)[,c("eqClassById","GENEID","txid")])
    minEquiRC$eqClassById <- unAsIs(minEquiRC$eqClassById)
    #colnames(minEquiRC)[1] <- "eqClassByIdTemp"
   
    minEquiRCTemp <- minEquiRC  %>% 
        mutate(txidTemp = eqClassById) %>% 
        unnest(c(txidTemp)) %>%
        mutate(equal = ifelse(txid == txidTemp, TRUE, FALSE)) %>%
        mutate(txidTemp = ifelse(txid == txidTemp, -1*as.numeric(txidTemp),txidTemp)) %>%
        group_by(eqClassById, txid) %>%
        mutate(eqClassByIdTemp = list(sort(unique(txidTemp)))) %>%
        ungroup() %>% 
        mutate(txid = abs(txidTemp), eqClassById = eqClassByIdTemp, eqClassByIdTemp = NULL, txidTemp = NULL) %>%
        distinct()
    
    minEquiRC <- minEquiRC %>% 
        mutate(txidTemp = eqClassById) %>% 
        unnest(c(txidTemp)) %>%
        mutate(minRC = 1, equal = FALSE, txid = txidTemp, txidTemp = NULL)
    
    minEquiRC <- bind_rows(minEquiRC, minEquiRCTemp)
    eqClassCount <- createEqClassToTxMapping(eqClassCount)
    eqClassCountJoin <- full_join(eqClassCount, minEquiRC, by = c("eqClassById","GENEID","txid","equal"))
    eqClassCountJoin[is.na(eqClassCountJoin)] <- 0
    eqClassCount_final <- eqClassCountJoin %>% 
        group_by(eqClassById) %>%
        mutate(nobs = max(nobs),
               rcWidth = max(rcWidth),
               minRC = max(minRC)) %>%
        ungroup() %>%
        distinct()
    return(eqClassCount_final)
}

#' Function to get rid of AsIs class so that group_by can be used on eqClassById column in mcols(annotations)
#' credit to https://stackoverflow.com/questions/12865218/getting-rid-of-asis-class-attribute
#' @noRd
unAsIs <- function(X) {
    if("AsIs" %in% class(X)) {
        class(X) <- class(X)[-match("AsIs", class(X))]
    }
    X
}


#' Create eqClass to tx mapping based on eqClassById
#' @import tidyr 
#' @noRd
createEqClassToTxMapping <- function(eqClassTable){
    eqClassTable_unnest <- eqClassTable %>% 
        mutate(txid = eqClassById) %>% 
        unnest(c(txid)) %>%
        mutate(equal = ifelse(txid<0,TRUE,FALSE)) %>%
        mutate(txid = abs(txid))
    return(eqClassTable_unnest)
}


#' modifiy incompatible read classes assignment
#' @import data.table
#' @noRd
modifyIncompatibleAssignment <- function(distTable){
    distTable <- data.table(as.data.frame(distTable))
    distTable[,`:=`(anyCompatible = any(compatible), 
                    anyEqual = any(equal)),
              by = readClassId]
    ## for uncompatible reads, change annotationTxId to gene_id_unidentified
    distTable[which(!anyCompatible), 
              annotationTxId := paste0(GENEID, 
                                       "_unidentified")]
    distTable[,`:=`(anyCompatible = NULL,
                    anyEqual = NULL)]
    distTable <- unique(DataFrame(setDF(distTable)))
    return(distTable)
}


 

#' Combine count se object while preserving the metadata objects
#' @noRd
combineCountSes <- function(countsSe, trackReads = FALSE, returnDistTable = FALSE){
    sampleNames = sapply(countsSe, FUN = function(x){colnames(x)})
    if(trackReads){
        readToTranscriptMaps = lapply(countsSe, FUN = function(se){metadata(se)$readToTranscriptMap})
        names(readToTranscriptMaps) = sampleNames
        countsSe = lapply(countsSe, FUN = function(se){
            metadata(se)$readToTranscriptMap=NULL
            return(se)})
    }
    if(returnDistTable){
        distTables = lapply(countsSe, FUN = function(se){metadata(se)$distTable})
        names(distTables) = sampleNames
        countsSe = lapply(countsSe, FUN = function(se){
            metadata(se)$distTable=NULL
            return(se)})
    }
    # combine incompatible counts
    incompatibleCounts = Reduce(merge_wrapper, lapply(countsSe, FUN = function(se){metadata(se)$incompatibleCounts}))
    countsSe = lapply(countsSe, FUN = function(se){
        metadata(se)$incompatibleCounts=NULL
        return(se)})
    countsSe <- do.call(SummarizedExperiment::cbind, countsSe)
    if(trackReads) metadata(countsSe)$readToTranscriptMaps = readToTranscriptMaps
    if(returnDistTable) metadata(countsSe)$distTables = distTables
    metadata(countsSe)$incompatibleCounts = incompatibleCounts
    return(countsSe)
}

# Quick wrapper function (https://stackoverflow.com/questions/13273833/merging-multiple-data-tables)
#' @noRd 
merge_wrapper <- function(x,y){
    merge.data.table(x,y,by = "GENEID",all=TRUE)
}

#' @useDynLib bambu, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL


#' @noRd
.onUnload <- function(libpath) {
    library.dynam.unload("bambu", libpath)
}

