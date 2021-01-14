#' Nanopore transcript abundance quantification
#' @title transcript_abundance_quantification
#' @param method A string variable indicates the whether a one-step or two-step
#' approach will be used. See \code{Details}
#' for details on one-step and two-step approach.
#' @param readClassDt A \code{data.table} with columns
#' @importFrom BiocParallel bplapply
#' @noRd
abundance_quantification <- function(readClassDt, ncore = 1, bias = TRUE,
    maxiter = 20000, conv = 10^(-2), minvalue = 10^(-8)) {
    gene_sidList <- unique(readClassDt$gene_sid)
    if (ncore == 1) {
        emResultsList <- lapply(as.list(gene_sidList),
            run_parallel,
            conv = conv,
            minvalue = minvalue, 
            bias = bias,
            maxiter = maxiter,
            readClassDt = readClassDt
        )
    } else {
        bpParameters <- BiocParallel::bpparam()
        bpParameters$workers <- ncore

        emResultsList <- BiocParallel::bplapply(as.list(gene_sidList),
            run_parallel,
            conv = conv,
            minvalue = minvalue, 
            bias = bias,
            maxiter = maxiter,
            readClassDt = readClassDt,
            BPPARAM = bpParameters
        )
    }
    estimates <- list(
        do.call("rbind", lapply(
            seq_along(emResultsList),
            function(x) emResultsList[[x]][[1]]
        )),
        do.call("rbind", lapply(
            seq_along(emResultsList),
            function(x) emResultsList[[x]][[2]]
        ))
    )
    return(estimates)
}


#' function to run in parallel
#' @title run_parallel
#' @param g the serial id of gene
#' @inheritParams abundance_quantification
#' @noRd
run_parallel <-
    function(g, conv, bias, minvalue, maxiter, readClassDt) {
    tmp <- unique(readClassDt[gene_sid == g])
    multiMap <- unique(tmp[, .(read_class_sid, multi_align)], 
            by = NULL)[order(read_class_sid)]$multi_align
    n.obs <- unique(tmp[, .(read_class_sid, nobs)], 
            by = NULL)[order(read_class_sid)]$nobs
    K <- as.numeric(sum(n.obs))
    n.obs <- as.numeric(n.obs)/K
    aMatList <- formatAmat(tmp, multiMap)
    a_mat <- aMatList[["combined"]]
    lambda <- sqrt(mean(n.obs))#suggested by Jiang and Salzman
    out <- initialiseOutput(a_mat, g, K, n.obs) 
    if (K == 0) return(out)
    #The following steps clean up the rc to tx matrix
    #step1:removes transcripts without any read class support
    rids <- which(apply(t(t(a_mat) * n.obs * K),1,sum) != 0)
    a_mat <- a_mat[rids,]
    #step2: removes read classes without transcript assignment after step1
    if (is(a_mat, "numeric")) {
        cids <- 1
    }else{
        cids <- which(apply(a_mat,2,sum) != 0)
        a_mat <- a_mat[,cids]
    }
    n.obs <- n.obs[cids]
    aMatList[["full"]] <- aMatList[["full"]][rids, cids]
    aMatList[["partial"]] <- aMatList[["partial"]][rids, cids]
    aMatList[["unique"]] <- aMatList[["unique"]][rids, cids]
    if (is.null(nrow(a_mat))) {
        out[[1]][rids]$counts <- K*n.obs*a_mat
        out[[1]][rids]$FullLengthCounts <- K*n.obs*aMatList$full
        out[[1]][rids]$PartialLengthCounts <- K*n.obs*aMatList$partial
        out[[1]][rids]$UniqueCounts <- K*n.obs*aMatList$unique
    }else{
        est_output <- emWithL1(X = as.matrix(a_mat), Y = n.obs,
            lambda = lambda, d = bias, maxiter = maxiter,
            minvalue = minvalue, conv = conv)
        out <- updateOutput(est_output, a_mat, rids, 
            cids,out,aMatList, K, n.obs)
    }
    return(out)
}

# This function generates the wide-format a matrix
#' @noRd
formatAmat <- function(tmp, multiMap){
        tmp_wide <- dcast(tmp[order(nobs)],  tx_ori + 
            fullTx ~ read_class_sid, value.var = "aval")
        tmp_wide <- tmp_wide[CJ(tx_ori = unique(tmp_wide$tx_ori),
            fullTx = c(TRUE,FALSE)), on = c("tx_ori", "fullTx")]
        tmp_wide[is.na(tmp_wide)] <- 0
        a_mat_full <- setDF(tmp_wide[which(fullTx)][,-c(1:2),with = FALSE])
        a_mat_partial <- setDF(tmp_wide[which(!fullTx)][,-c(1:2),with = FALSE])
        a_mat_unique <- a_mat <- a_mat_full + a_mat_partial
        a_mat_unique[, which(multiMap)] <- 0
        rownames(a_mat) <- rownames(a_mat_partial) <- rownames(a_mat_full) <-
            rownames(a_mat_unique) <- unique(tmp_wide$tx_ori)
        return(list(combined = a_mat, 
                    full = a_mat_full, 
                    partial = a_mat_partial,
                    unique = a_mat_unique))
}

# This function generates a_mat values for all transcripts 
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



# This function initialises the final estimates with default values
#' @noRd 
initialiseOutput <- function(a_mat, g, K, n.obs){
    return(list(data.table(tx_sid = rownames(a_mat),counts = 0,
                FullLengthCounts = 0, PartialLengthCounts = 0,
                UniqueCounts = 0, theta = 0, gene_sid = g, 
                ntotal = as.numeric(K)),
                data.table(read_class_sid = colnames(a_mat),biasEstimates = 0,
                gene_sid = g,ntotal = as.numeric(K))))#pre-define output
}

# This function generates the estimates for overall, full, partial, and unique
# estimates
#' @noRd
updateOutput <- function(est_output, a_mat, rids, cids,
    out, aMatList, K, n.obs){
    est_output[["counts"]] <- 
        getFullandPartialEstimates(a_mat, a_mat, est_output, n.obs)
    est_output[["FullLengthCounts"]] <- 
        getFullandPartialEstimates(aMatList[["full"]], a_mat, est_output, n.obs)
    est_output[["PartialLengthCounts"]] <- getFullandPartialEstimates(
        aMatList[["partial"]], a_mat, est_output, n.obs)
    est_output[["UniqueCounts"]] <- getFullandPartialEstimates(
        a_mat = aMatList[["unique"]], a_mat_sum = a_mat, est_output, n.obs)
    out <- modifyQuantOut(est_output, a_mat, rids, cids, out, K)
    return(out)
}

#' Calculate degradation rate based on equiRC read counts 
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
        ((geneCountLength$nobs - geneCountLength$dObs) >= 5))) == 0){
        warning("There is not enough read count and full length coverage!
            Hence degradation rate is estimated using all data!")
    } else {
         geneCountLength <- geneCountLength[nobs >= 30 & ((nobs - dObs) >= 5)]
    }
    d_rate <- median(geneCountLength$d_rate * 1000/geneCountLength$gene_len,
        na.rm = TRUE)
    return(list(d_rate, nrow(geneCountLength)))
}


# Modify default quant output using estimated outputs
#' @noRd
modifyQuantOut <- function(est_output, a_mat, rids, cids, out,K){
    out[[1]][rids]$theta <- est_output[["theta"]]
    out[[1]][rids]$counts <- K * est_output[["counts"]]
    out[[1]][rids]$FullLengthCounts <- K * est_output[["FullLengthCounts"]]
    out[[1]][rids]$PartialLengthCounts <- 
        K * est_output[["PartialLengthCounts"]]
    out[[1]][rids]$UniqueCounts <- K * est_output[["UniqueCounts"]]
    out[[2]][cids]$biasEstimates <- as.numeric(t(est_output[["b"]]))
    return(out)
}

# This function uses the E-step in the EM algorithm to compute the expected
# relative read count attributed to each transcript from each read class, note 
# there are a_mat, one is the overall sampling matrix, i.e., a_mat_sum
#' @param a_mat This is the sampling matrix of the target estimates
#' @param a_mat_sum This is the sampling matrix of the overall estimates, it 
#' will be mainly used in normalizing the count.
#' @noRd
getFullandPartialEstimates <- function(a_mat, a_mat_sum, est_output, 
    n.obs){
    atheta <- a_mat * est_output[["theta"]]
    atheta_sum <- a_mat_sum * est_output[["theta"]]
    return(colSums(n.obs * t(atheta) / colSums(atheta_sum), na.rm = TRUE))
}

# This function converts transcript, gene, and read class names to simple
# integers for more efficient computation
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
#' @noRd
formatOutput <- function(outList, ori_txvec, geneVec){
    theta_est <- outList[[1]]
    theta_est[, `:=`(tx_name = ori_txvec[as.numeric(tx_sid)],
        gene_name = geneVec[gene_sid])]
    theta_est[, `:=`(tx_sid = NULL, gene_sid = NULL)]
    theta_est <- theta_est[, .(tx_name, counts,FullLengthCounts,
        PartialLengthCounts, UniqueCounts, theta)]
    totalCount <- sum(theta_est$counts)
    theta_est[, `:=`(CPM = counts / totalCount * (10^6))]
    return(theta_est)
}


## Remove duplicate transcript counts originated from multiple genes
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

#' @useDynLib bambu, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @noRd
.onUnload <- function(libpath) {
    library.dynam.unload("bambu", libpath)
}

