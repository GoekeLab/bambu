#' Nanopore transcript abundance quantification
#' @title transcript_abundance_quantification
#' @param method A string variable indicates the whether a one-step or two-step approach will be used. See \code{Details}
#' for details on one-step and two-step approach.
#' @param read_classDT A \code{data.table} with columns
#' @examples
#' \dontrun{
#' system.file("extdata", "example_data.rds", package = "bamboo")
#' abundance_quantification(example_data[[1]])
#' abundance_quantification(example_data[[2]])
#' abundance_quantification(example_data[[3]])
#' abundance_quantification(example_data[[4]])
#' }
abundance_quantification <- function(read_classDT,mc.cores = 1,
                                     method = 'two-step',
                                     conv.control = 10^(-3)){
  gene_sidList <- unique(read_classDT$gene_sid)
  n_gene <- length(gene_sidList)
  n_chunk <- ceiling(n_gene/2000)
  splitIndex <- split(1:n_gene, rep(1:n_chunk, ceiling(n_gene/n_chunk))[1:n_gene])
  outData <- lapply(as.list(1:n_chunk), function(i){
    gene_sid_chunk <- gene_sidList[splitIndex[[i]]]
    #print(i)
    outData <- parallel::mclapply(gene_sid_chunk,run_parallel,method = method, conv.control = conv.control, read_classDT = read_classDT,
                        mc.cores = mc.cores, mc.preschedule = TRUE)
    outDataList <- list(do.call('rbind',lapply(1:length(outData), function(x) outData[[x]][[1]])),
                        do.call('rbind',lapply(1:length(outData), function(x) outData[[x]][[2]])))
    return(outDataList)
  })
  estimates <- list(do.call('rbind',lapply(1:length(outData), function(x) outData[[x]][[1]])),
                    do.call('rbind',lapply(1:length(outData), function(x) outData[[x]][[2]])))

  return(estimates)
}

run_parallel <- function(g,method,conv.control,read_classDT){
  tmp <- read_classDT[gene_sid==g]
  if((nrow(tmp)==1)|(sum(tmp$nobs)==0)){
    outData <- list(data.table(tx_sid = tmp$tx_sid,
                               estimates = tmp$nobs,
                               gene_sid = tmp$gene_sid,
                               ntotal = tmp$nobs),
                    NULL)
  }else{
    tmp_wide <- dcast(tmp, tx_sid~read_class_sid, fun.agg = length,
                      value.var = 'nobs')
    a_mat <- tmp_wide[,-1,with=FALSE]
    setDF(a_mat)
    #a_mat <- a_mat/rowSums(a_mat)
    rownames(a_mat) <- tmp_wide$tx_sid

    n.obs <- unique(tmp[,.(read_class_sid,nobs)],by=NULL)$nobs

    ## using the Jiang and Salzman suggested value
    lambda <- sqrt(max(n.obs))#,25)

    ## two - step approach
    ##===========================
    ## Step 1 estimate with L1 penalty
    ## EM with L1
    theta.hat <- emWithL1(X = as.matrix(a_mat),
                          Y = n.obs,
                          lambda = lambda,
                          conv = conv.control)

    b_est <- as.numeric(t(theta.hat$b))
    ##===========================
    ## Step 2 estimate without L1 penalty
    ## EM without L1
    if(method == 'two-step'){

      theta.hat <-emWithoutL1(X = as.matrix(a_mat),
                               Y = n.obs,
                               par = b_est,#rep(0.5,ncol(a_mat)),#
                               lambda = lambda,
                               conv = conv.control)

    }
    t_est <- as.numeric(t(theta.hat$theta))
    nobs <- sum(n.obs)
    outData <- list(data.table(tx_sid = rownames(a_mat),
                               estimates = t_est,
                               gene_sid = g,
                               ntotal = nobs),
                    data.table(read_class_sid = colnames(a_mat),
                               estimates = b_est,
                               gene_sid = g,
                               ntotal = nobs))
  }
  return(outData)
}
