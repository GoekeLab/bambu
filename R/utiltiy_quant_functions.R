#' Nanopore transcript abundance quantification
#' @title transcript_abundance_quantification
#' @param method A string variable indicates the whether a one-step or two-step
#' approach will be used. See \code{Details}
#' for details on one-step and two-step approach.
#' @param readClassDt A \code{data.table} with columns
#' @importFrom BiocParallel bplapply
#' @noRd
abundance_quantification <- function(readClassDt, ncore = 1,
    bias = TRUE, maxiter = 20000, conv = 10^(-8), d_rate = 0.1) {
    gene_sidList <- unique(readClassDt$gene_sid)
    if (ncore == 1) {
        emResultsList <- lapply(as.list(gene_sidList),
            run_parallel,
            conv = conv,
            bias = bias,
            maxiter = maxiter,
            d_rate = d_rate,
            readClassDt = readClassDt
        )
    } else {
        bpParameters <- BiocParallel::bpparam()
        bpParameters$workers <- ncore

        emResultsList <- BiocParallel::bplapply(as.list(gene_sidList),
            run_parallel,
            conv = conv,
            bias = bias,
            maxiter = maxiter,
            d_rate = d_rate,
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
run_parallel <- function(g, conv, bias, maxiter, d_rate, readClassDt) {
    tmp <- unique(readClassDt[gene_sid == g])
    if ((nrow(tmp) == 1)) {
        out <- list(data.table(tx_sid = tmp$tx_sid,estimates = tmp$nobs,
                gene_sid = tmp$gene_sid,ntotal = sum(tmp$nobs)),
                data.table(read_class_sid = tmp$read_class_sid,
                counts = 0, FullLengthCounts = 0, PartialLengthCounts = 0,
                UniqueCounts = 0, n = tmp$nobs,gene_sid = g,
                ntotal = sum(tmp$nobs)))
    } else {
        n.obs <- unique(tmp[, .(read_class_sid, nobs)], 
            by = NULL)[order(read_class_sid)]$nobs
        K <- sum(n.obs)
        n.obs <- n.obs/K
        aMatList <- formatAmat(tmp, d_rate)
        a_mat <- aMatList[["combined"]]
        lambda <- sqrt(max(n.obs)) ##using the Jiang and Salzman suggested value
        out <- initialiseOutput(a_mat, g, K, n.obs) 
        ## The following two steps clean up the rc to tx matrix so that it won't 
        ## give NA estimates
        ## note: this step removes transcripts without any read class support
        removed_txs <- which(apply(t(t(a_mat) * n.obs * K),1,sum) == 0)
        if (length(removed_txs) > 0) a_mat <- a_mat[-removed_txs,]
        ## note: this tep removes read classes without any transcript assignment 
        ## after the first step
        removed_rcs <- which(apply(a_mat,2,sum) == 0)
        if (length(removed_rcs) > 0) { 
            a_mat <- a_mat[,-removed_rcs]
            n.obs <- n.obs[-removed_rcs]
        }
        aMatList[["full"]] <- 
            updateAmat(aMatList[["full"]], removed_txs, removed_rcs)
        aMatList[["partial"]] <- 
            updateAmat(aMatList[["partial"]], removed_txs, removed_rcs)
        aMatList[["unique"]] <- 
            updateAmat(aMatList[["unique"]], removed_txs, removed_rcs)
        if (is.null(nrow(a_mat))) {
            out[[1]][-as.numeric(removed_txs)]$counts <- K*n.obs
        }else{
            est_output <- emWithL1(X = as.matrix(a_mat), Y = n.obs, 
                lambda = lambda, d = bias, maxiter = maxiter, conv = conv)
            out <- updateOutput(est_output, a_mat, bias, removed_rcs, 
                removed_txs,out,K,aMatList,n.obs)
        }
    }
    return(out)
}

#' @noRd
formatAmat <- function(tmp, d_rate){
        tmp[, multi_align := (length(unique(tx_sid)) > 1), by = read_class_sid]
        tmp[which(multi_align) , aval := ifelse(fullTx, 1 - 
            sum(.SD[which(!fullTx)]$rc_width*d_rate/1000),
            rc_width*d_rate/1000), by = tx_ori]
        tmp[multi_align & fullTx, aval := aval + rc_width*d_rate/1000]
        tmp[, ave_status := any(fullTx & multi_align), by = tx_ori]
        tmp[which(ave_status), aval := aval/sum(aval,na.rm = TRUE), by = tx_ori]
        tmp[which(!multi_align), aval := 1]
        tmp_wide <- dcast(tmp[order(nobs)],  tx_ori + 
            fullTx ~ read_class_sid, value.var = "aval")
        tmp_wide <- tmp_wide[CJ(tx_ori = unique(tmp_wide$tx_ori),
            fullTx = c(TRUE,FALSE)), on = c("tx_ori", "fullTx")]
        tmp_wide[is.na(tmp_wide)] <- 0
        a_mat_full <- setDF(tmp_wide[which(fullTx)][,-c(1:2),with = FALSE])
        a_mat_partial <- setDF(tmp_wide[which(!fullTx)][,-c(1:2),with = FALSE])
        a_mat_unique <- a_mat <- a_mat_full + a_mat_partial
        a_mat_unique[,apply(a_mat_unique > 0,2,sum) > 1] <- 0
        rownames(a_mat) <- rownames(a_mat_partial) <- rownames(a_mat_full) <-
            rownames(a_mat_unique) <- unique(tmp_wide$tx_ori)
        return(list(combined = a_mat, 
                    full = a_mat_full, 
                    partial = a_mat_partial,
                    unique = a_mat_unique))
}

#' @noRd
updateAmat <- function(a_mat, removed_txs, removed_rcs){
    if (length(removed_txs) > 0)  a_mat <- a_mat[-removed_txs,]
    if (length(removed_rcs) > 0)  a_mat <- a_mat[,-removed_rcs]
    return(a_mat)
}

#' @noRd 
initialiseOutput <- function(a_mat, g, K, n.obs){
    return(list(data.table(tx_sid = rownames(a_mat),counts = 0,
                FullLengthCounts = 0, PartialLengthCounts = 0,
                UniqueCounts = 0, gene_sid = g,ntotal = K),
                data.table(read_class_sid = colnames(a_mat),biasEstimates = 0,
                n = n.obs, gene_sid = g,ntotal = K)))#pre-define output
}

updateOutput <- function(est_output, a_mat, bias, 
    removed_rcs, removed_txs,out,K,aMatList,n.obs){
    est_output[["counts"]] <- 
        getFullandPartialEstimates(a_mat, a_mat, est_output, n.obs)
    est_output[["FullLengthCounts"]] <- 
        getFullandPartialEstimates(aMatList[["full"]], a_mat, est_output, n.obs)
    est_output[["PartialLengthCounts"]] <- getFullandPartialEstimates(
        aMatList[["partial"]], a_mat, est_output, n.obs)
    est_output[["UniqueCounts"]] <- getFullandPartialEstimates(
        aMatList[["unique"]], a_mat, est_output, n.obs)
    out <- modifyQuantOut(est_output, a_mat, bias, 
        removed_rcs, removed_txs,out,K)
}

# Modify default quant output using estimated outputs
#' @noRd
modifyQuantOut <- function(est_output, a_mat,
    bias, removed_rcs, removed_txs, out,K){
    if (bias) {
        b_est <- as.numeric(t(est_output[["b"]]))
        if (length(removed_rcs) > 0) {
            out[[2]][-removed_rcs]$biasEstimates <- b_est
        }else{
            out[[2]]$biasEstimates <- b_est
        }
    } 
    #t_est <- as.numeric(t(K * apply(a_mat * est_output[["theta"]], 1, sum)))
    t_est <- K * est_output[["counts"]]
    t_est_full <- K * est_output[["FullLengthCounts"]]
    t_est_partial <- K * est_output[["PartialLengthCounts"]]
    t_est_unique <- K * est_output[["UniqueCounts"]]
    if (length(removed_txs) > 0) {
        out[[1]][-removed_txs]$counts <- t_est
        out[[1]][-removed_txs]$FullLengthCounts <- t_est_full
        out[[1]][-removed_txs]$PartialLengthCounts <- t_est_partial
        out[[1]][-removed_txs]$UniqueCounts <- t_est_unique
    }else{
        out[[1]]$counts <- t_est
        out[[1]]$FullLengthCounts <- t_est_full
        out[[1]]$PartialLengthCounts <- t_est_partial
        out[[1]]$UniqueCounts <- t_est_unique
    }
    return(out)
}

#' @noRd
getFullandPartialEstimates <- function(a_mat, a_mat_sum, est_output, 
    n.obs){
    atheta <- a_mat*est_output[["theta"]]
    atheta_sum <- a_mat_sum*est_output[["theta"]]
    return(colSums(n.obs*t(atheta)/colSums(atheta_sum), na.rm = TRUE))
}

#' @noRd
simplifyNames <- function(readClassDt, txVec, geneVec,ori_txvec){
    readclassVec <- unique(readClassDt$read_class_id)
    readClassDt <- as.data.table(readClassDt)
    readClassDt[, gene_sid := match(gene_id, geneVec)]
    readClassDt[, tx_sid := match(tx_id, txVec)]
    readClassDt[, tx_ori := match(gsub("Start","",tx_id),ori_txvec)]
    readClassDt[, read_class_sid := match(read_class_id, readclassVec)]
    readClassDt[, fullTx := grepl("Start",tx_id)]
    readClassDt[, `:=`(tx_id = NULL, gene_id = NULL, read_class_id = NULL)]
    return(readClassDt)
}

#' Format output 
#' @noRd
formatOutput <- function(outList, ori_txvec, geneVec){
    theta_est <- outList[[1]]
    theta_est[, `:=`(tx_name = ori_txvec[as.numeric(tx_sid)],
        gene_name = geneVec[gene_sid])]
    theta_est[, `:=`(tx_sid = NULL, gene_sid = NULL)]
    theta_est <- theta_est[, .(tx_name, counts,FullLengthCounts,
        PartialLengthCounts, UniqueCounts)]
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
                            CPM = sum(CPM)), by = tx_name],by = NULL)
    return(counts_final)
}
