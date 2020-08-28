#' Nanopore transcript abundance quantification
#' @title transcript_abundance_quantification
#' @param method A string variable indicates the whether a one-step or two-step approach will be used. See \code{Details}
#' for details on one-step and two-step approach.
#' @param readClassDt A \code{data.table} with columns
#' @importFrom BiocParallel bplapply
#' @noRd
abundance_quantification <- function(readClassDt,ncore = 1,
                                     bias = TRUE,
                                     maxiter = 20000,
                                     conv = 10^(-8)){
  gene_sidList <- unique(readClassDt$gene_sid)
 
    if(ncore == 1){
      
      emResultsList <- lapply(as.list(gene_sidList),
                                              run_parallel,
                                              conv = conv,
                                              bias = bias,
                                              maxiter = maxiter,
                                              readClassDt = readClassDt)
    }else{
      bpParameters <- BiocParallel::bpparam()
      bpParameters$workers <- ncore
      
      emResultsList <- BiocParallel::bplapply(as.list(gene_sidList),
                                              run_parallel,
                                              conv = conv,
                                              bias = bias,
                                              maxiter = maxiter,
                                              readClassDt = readClassDt,
                                              BPPARAM=bpParameters)
    }
    
  


  estimates <- list(do.call('rbind',lapply(seq_along(emResultsList), function(x) emResultsList[[x]][[1]])),
                    do.call('rbind',lapply(seq_along(emResultsList), function(x) emResultsList[[x]][[2]])))

  return(estimates)
}

run_parallel <- function(g,conv,bias,maxiter, readClassDt){
  tmp <- readClassDt[gene_sid==g]
  if((nrow(tmp)==1)){
    out <- list(data.table(tx_sid = tmp$tx_sid,
                               estimates = tmp$nobs,
                               gene_sid = tmp$gene_sid,
                               ntotal = sum(tmp$nobs)),
                    data.table(read_class_sid = tmp$read_class_sid,
                               estimates = 0,
                               n = tmp$nobs,
                               gene_sid = g,
                               ntotal = sum(tmp$nobs)))
  }else{
    tmp_wide <- dcast(tmp[order(nobs)], tx_sid~read_class_sid, fun.aggregate = length,
                      value.var = 'nobs')
    a_mat <- tmp_wide[,-1,with=FALSE]
    setDF(a_mat)
    #a_mat <- a_mat/rowSums(a_mat)
    rownames(a_mat) <- tmp_wide$tx_sid

    n.obs <- unique(tmp[,.(read_class_sid,nobs)],by=NULL)$nobs


    ## using the Jiang and Salzman suggested value
    lambda <- sqrt(max(n.obs))#,25)
    est_output <- emWithL1(X = as.matrix(a_mat),
                          Y = n.obs,
                          lambda = lambda,
                          d = bias,
                          maxiter = maxiter,
                          conv = conv)
    t_est <- as.numeric(t(est_output[["theta"]]))

    if(bias){
      b_est <- as.numeric(t(est_output[["b"]]))
    }else{
      b_est <- rep(0,ncol(a_mat))
    }
    nobs <- sum(n.obs)
    out <- list(data.table(tx_sid = rownames(a_mat),
                               estimates = t_est,
                               gene_sid = g,
                               ntotal = nobs),
                    data.table(read_class_sid = colnames(a_mat),
                               estimates = b_est,
                               n = n.obs,
                               gene_sid = g,
                               ntotal = nobs))
  }
  return(out)
}
