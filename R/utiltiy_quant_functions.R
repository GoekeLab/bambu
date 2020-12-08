#' Nanopore transcript abundance quantification
#' @title transcript_abundance_quantification
#' @param method A string variable indicates the whether a one-step or two-step
#' approach will be used. See \code{Details}
#' for details on one-step and two-step approach.
#' @param readClassDt A \code{data.table} with columns
#' @importFrom BiocParallel bplapply
#' @noRd
abundance_quantification <- function(readClassDt, ncore = 1,
    bias = TRUE, maxiter = 20000, conv = 10^(-8)) {
    gene_sidList <- unique(readClassDt$gene_sid)
    if (ncore == 1) {
        emResultsList <- lapply(as.list(gene_sidList),
            run_parallel,
            conv = conv,
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
run_parallel <- function(g, conv, bias, maxiter, readClassDt) {
    tmp <- unique(readClassDt[gene_sid == g])
    if ((nrow(tmp) == 1)) {
        out <- list(data.table(tx_sid = tmp$tx_sid,estimates = tmp$nobs,
                gene_sid = tmp$gene_sid,ntotal = sum(tmp$nobs)),
                data.table(read_class_sid = tmp$read_class_sid,
                estimates = 0,n = tmp$nobs,gene_sid = g,ntotal = sum(tmp$nobs)))
    } else {
        n.obs <- unique(tmp[, .(read_class_sid, nobs)], by = NULL)$nobs
        K <- sum(n.obs)
        n.obs <- n.obs/K
        tmp_wide <- dcast(tmp[order(nobs)], tx_sid ~ read_class_sid,
            fun.aggregate = length, value.var = "nobs")
        a_mat <- tmp_wide[, -1, with = FALSE]
        setDF(a_mat)
        rownames(a_mat) <- tmp_wide$tx_sid
        lambda <- sqrt(max(n.obs)) ##using the Jiang and Salzman suggested value
        out <- list(data.table(tx_sid = rownames(a_mat),estimates = 0,
                gene_sid = g,ntotal = K),
            data.table(read_class_sid = colnames(a_mat),estimates = 0,
                n = n.obs, gene_sid = g,ntotal = K))#pre-define output
        removed_txs <- which(apply(t(t(a_mat) * n.obs*K),1,sum) == 0)
        if (length(removed_txs) > 0) a_mat <- a_mat[-removed_txs,]
        removed_rcs <- which(apply(a_mat,2,sum) == 0)
        if (length(removed_rcs) > 0) { 
            a_mat <- a_mat[,-removed_rcs]
            n.obs <- n.obs[-removed_rcs]
        }
        est_output <- emWithL1(X = as.matrix(a_mat), Y = n.obs, lambda = lambda,
            d = bias, maxiter = maxiter, conv = conv)
        est_output$theta <- K*apply(a_mat * est_output$theta,1,sum)
        out <- modifyQuantOut(est_output, bias, removed_rcs, removed_txs,out)
    }
    return(out)
}

# Modify default quant output using estimated outputs
#' @noRd
modifyQuantOut <- function(est_output, bias, removed_rcs, removed_txs,out){
    if (bias) {
        b_est <- as.numeric(t(est_output[["b"]]))
        if (length(removed_rcs) > 0) {
            out[[2]][-removed_rcs]$estimates <- b_est
        }else{
            out[[2]]$estimates <- b_est
        }
    } 
    t_est <- as.numeric(t(est_output[["theta"]]))
    if (length(removed_txs) > 0) {
        out[[1]][-removed_txs]$estimates <- t_est
    }else{
        out[[1]]$estimates <- t_est
    }
    return(out)
}

