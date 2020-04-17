#' Nanopore transcript abundance quantification
#' @title transcript_abundance_quantification
#' @param method A string variable indicates the whether a one-step or two-step approach will be used. See \code{Details}
#' for details on one-step and two-step approach.
#' @param read_classDT A \code{data.table} with columns
#' @examples
#' @noRd
abundance_quantification <- function(read_classDT,mc.cores = 1,
                                     bias_correction = TRUE,
                                     maxiter = 20000,
                                     conv.control = 10^(-8)){
  gene_sidList <- unique(read_classDT$gene_sid)
  n_gene <- length(gene_sidList)
  n_chunk <- ceiling(n_gene/2000)
  splitIndex <- split(1:n_gene, rep(1:n_chunk, ceiling(n_gene/n_chunk))[1:n_gene])
  pb <- progress::progress_bar$new(
    format = "  Genes [:bar] :percent in :elapsed",
    total = n_chunk, clear = FALSE, width= 60)
  pb$tick(0)

  outData <- lapply(as.list(1:n_chunk), function(i){
    gene_sid_chunk <- gene_sidList[splitIndex[[i]]]
    outData <- parallel::mclapply(gene_sid_chunk,run_parallel,
                                  conv.control = conv.control,
                                  bias_correction = bias_correction,
                                  maxiter = maxiter, read_classDT = read_classDT,
                                  mc.set.seed = TRUE,mc.preschedule = TRUE,
                                  mc.silent = FALSE, mc.cores = mc.cores)

    outData <- list(do.call('rbind',lapply(1:length(outData), function(x) outData[[x]][[1]])),
                    do.call('rbind',lapply(1:length(outData), function(x) outData[[x]][[2]])))
    pb$tick()
    return(outData)
  })

  estimates <- list(do.call('rbind',lapply(1:length(outData), function(x) outData[[x]][[1]])),
                    do.call('rbind',lapply(1:length(outData), function(x) outData[[x]][[2]])))

  return(estimates)
}

run_parallel <- function(g,conv.control,bias_correction,maxiter, read_classDT){
  tmp <- read_classDT[gene_sid==g]
  if((nrow(tmp)==1)){
    outData <- list(data.table(tx_sid = tmp$tx_sid,
                               estimates = tmp$nobs,
                               gene_sid = tmp$gene_sid,
                               ntotal = sum(tmp$nobs)),
                    data.table(read_class_sid = tmp$read_class_sid,
                               estimates = 0,
                               n = tmp$nobs,
                               gene_sid = g,
                               ntotal = sum(tmp$nobs)))
  }else{
    tmp_wide <- dcast(tmp[order(nobs)], tx_sid~read_class_sid, fun.agg = length,
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
                          d = bias_correction,
                          maxiter = maxiter,
                          conv = conv.control)
    t_est <- as.numeric(t(est_output[["theta"]]))

    if(bias_correction){
      b_est <- as.numeric(t(est_output[["b"]]))
    }else{
      b_est <- rep(0,ncol(a_mat))
    }
    nobs <- sum(n.obs)
    outData <- list(data.table(tx_sid = rownames(a_mat),
                               estimates = t_est,
                               gene_sid = g,
                               ntotal = nobs),
                    data.table(read_class_sid = colnames(a_mat),
                               estimates = b_est,
                               n = n.obs,
                               gene_sid = g,
                               ntotal = nobs))
  }
  return(outData)
}
