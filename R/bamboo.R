bamboo <- function(bam.file = NULL, se = NULL, dt = NULL,txdb = NULL, txdbTablesList = NULL, algo.control = NULL, ir.control = NULL, fa.file = NULL){

  if(!is.null(dt)){
    return(bamboo.quantDT(dt = dt,algo.control = algo.control))
  }

  if(is.null(txdbTablesList)){
    if(!is.null(txdb)) # txdb object by reading txdb file
    {
      if(class(txdb) != 'TxDb'){
        stop("txdb object is missing.")
      }
      txdbTablesList <- prepareAnnotations(txdb) ## note: check which annotation tables are reused multiple times, and inlcude only those. Optimise required annotations if possible
    }else{
      stop("txdb object is missing.")
    }
  }

    if(!is.null(se)){
      return(bamboo.quantSE(se = se,txdb = txdb, txdbTablesList = txdbTablesList, algo.control = algo.control))
          }
    if(!is.null(bam.file)){
      return(bamboo.quantISORE(bam.file = bam.file,algo.control = algo.control, fa.file=fa.file,
                               txdb=txdb,
                               txdbTablesList=txdbTablesList, ir.control = ir.control))
    }

  stop("At least bam.file or summarizedExperiment output from isore need to be provided.")
}

bamboo.quantDT <- function(dt = dt,algo.control = NULL){
  if(is.null(dt)){
    stop("Input object is missing.")
  }else if(any(!(c('gene_id','tx_id','read_class_id','nobs') %in% colnames(dt)))){
    stop("Columns gene_id, tx_id, read_class_id, nobs, are missing from object.")
  }

  ## check quantification parameters
  default.algo.control <- list(ncore = parallel::detectCores(),
                               method = "two-step",
                               convcontrol = 10^(-3))

  if(is.null(algo.control)){
    algo.control <- default.algo.control
  }else{
    if(is.null(algo.control[["ncore"]])){
      algo.control[["ncore"]] <- default.algo.control[["ncore"]]
    }else if(as.numeric(algo.control[["ncore"]])>default.algo.control[["ncore"]]){
      algo.control[["ncore"]] <- default.algo.control[["ncore"]]
    }
    if(is.null(algo.control[["method"]])){
      algo.control[["method"]] <- default.algo.control[["method"]]
    }else if(!(algo.control[["method"]] %in% c("one-step","two-step"))){
      algo.control[["method"]] <- default.algo.control[["method"]]
    }
    if(is.null(algo.control[["convcontrol"]])){
      algo.control[["convcontrol"]] <- default.algo.control[["convcontrol"]]
    }
  }

  ##----step2: match to simple numbers to increase claculation efficiency
  geneVec <- unique(dt$gene_id)
  txVec <- unique(dt$tx_id)
  readclassVec <- unique(dt$read_class_id)
  dt[, gene_sid:=match(gene_id, geneVec)]
  dt[, tx_sid:=match(tx_id, txVec)]
  dt[, read_class_sid:=match(read_class_id, readclassVec)]

  dt[,`:=`(tx_id = NULL, gene_id = NULL, read_class_id = NULL)]

  ##----step3: aggregate read class
  temp <- aggReadClass(dt)
  dt <- temp[[1]]
  eqClassVec <- temp[[2]]

  ##----step4: quantification
  start.time <- proc.time()
  outList <- abundance_quantification(dt,
                                      mc.cores = algo.control[["ncore"]],
                                      method = algo.control[["method"]],
                                      conv.control = algo.control[["convcontrol"]])
  end.time <- proc.time()
  cat(paste0('Finished EM estimation in ', round((end.time-start.time)[3]/60,1), ' mins', ' \n'))


  theta_est <- outList[[1]]
  theta_est[, `:=`(tx_name = txVec[as.numeric(tx_sid)],
                   gene_name = geneVec[gene_sid] )]
  theta_est[,`:=`(tx_sid=NULL, gene_sid=NULL)]


  b_est <- outList[[2]]
  b_est[, `:=`(gene_name = geneVec[gene_sid], eqClass = eqClassVec[as.numeric(read_class_sid)])]
  b_est[, `:=`(gene_sid = NULL,read_class_sid=NULL)]




  theta_est <- theta_est[,.(tx_name, estimates)]
  theta_est[,`:=`(CPM = estimates/sum(estimates)*(10^6))]




  est.list <- list(counts = theta_est, metadata = b_est)
  return(est.list)
}





bamboo.quantSE <- function(se = se,txdb = NULL, txdbTablesList = NULL, algo.control = NULL){
  if(is.null(txdbTablesList)){
    if(!is.null(txdb)) # txdb object by reading txdb file
    {
      if(class(txdb) != 'TxDb'){
        stop("txdb object is missing.")
      }
      txdbTablesList <- prepareAnnotations(txdb) ## note: check which annotation tables are reused multiple times, and inlcude only those. Optimise required annotations if possible
    }else{
      stop("txdb object is missing.")
    }
  }

  dt <- getEmptyClassFromSE(se, txdbTablesList)

  ## To do:
  ## task1: optional: to implement filtering function
  ## task2: to implement for multiple samples, when multiple samples are provided, run txdbtableslist for one time
  est <- bamboo.quantDT(dt,algo.control = algo.control)
  counts <- est$counts
  ColNames <- colnames(se)
  counts <- merge(counts,data.table(tx_name = names(txdbTablesList$exonsByTx)), all = TRUE,  on = 'tx_name')
  counts[is.na(estimates),`:=`(estimates = 0, CPM = 0) ]
  counts <- setDF(counts)
  seOutput <- SummarizedExperiment::SummarizedExperiment(assays = SimpleList(estimates = matrix(counts$estimates,ncol = length(ColNames), dimnames = list(counts$tx_name, ColNames)),
                                                                       normEstimates = matrix(counts$CPM, ncol =  length(ColNames), dimnames = list(counts$tx_name, ColNames))),
                                   rowRanges = txdbTablesList$exonsByTx[counts$tx_name],
                                   colData = colData(se),
                                   metadata = est$metadata) # transcript annotation with read class information
  return(seOutput)
}



bamboo.quantISORE <- function(bam.file = bam.file, algo.control = NULL, fa.file=NULL, txdb=NULL, txdbTablesList=NULL, ir.control = NULL, yieldSize = NULL, quickMode = FALSE){
  if(is.null(fa.file)){
    stop("Genome fa file is missing!")
  }
  if(is.null(ir.control)){
    ir.control <- list(stranded = FALSE,
                       protocol = NULL,
                       prefix = '',
                       minimumReadSupport = 2,
                       minimumTxFraction = 0.02)
  }else{
    if(is.null(ir.control[['stranded']])){
      ir.control[['stranded']] <- FALSE
    }
    if(is.null(ir.control[['prefix']])){
      ir.control[['prefix']] <- ''
    }
    if(is.null(ir.control[['minimumReadSupport']])){
      ir.control[['minimumReadSupport']] <- 2
    }
    if(ir.control[['minimumTxFraction']]){
      ir.control[['minimumTxFraction']] <- 0.02
    }
  }

  if(is.null(txdbTablesList)){
    if(!is.null(txdb)) # txdb object by reading txdb file
    {
      if(class(txdb) != 'TxDb'){
        stop("txdb object is missing.")
      }
      txdbTablesList <- prepareAnnotations(txdb) ## note: check which annotation tables are reused multiple times, and inlcude only those. Optimise required annotations if possible
    }else{
      stop("txdb object is missing.")
    }
  }

  ## create BamFileList object from character ##
  if(class(bam.file)=='BamFile') {
    if(!is.null(yieldSize)) {
      yieldSize(bam.file) <- yieldSize
    } else {
      yieldSize <- yieldSize(bam.file)
    }
    bam.file <- Rsamtools::BamFileList(bam.file)
  }else if(class(bam.file)=='BamFileList') {
    if(!is.null(yieldSize)) {
      yieldSize(bam.file) <- yieldSize
    } else {
      yieldSize <- min(yieldSize(bam.file))
    }
  }else if(any(!grepl('\\.bam$',bam.file))){
    stop("Bam file is missing from arguments.")
  }else{
    if(is.null(yieldSize)) {
      yieldSize <- NA
    }
    bam.file <- Rsamtools::BamFileList(bam.file, yieldSize = yieldSize)
  }

  if(extendAnnotations==FALSE) {
    for(bam.file.index in seq_along(bam.file)){
      start.time <- proc.time()
      se  <- isore(bamFile = bam.file[[bam.file.index]],
                   txdb = NULL,
                   txdbTablesList = txdbTablesList,
                   genomeFA = fa.file,
                   stranded = ir.control[['stranded']],
                   protocol = ir.control[['protocol']],
                   prefix = ir.control[['prefix']],
                   minimumReadSupport= ir.control[['minimumReadSupport']],
                   minimumTxFraction = ir.control[['minimumTxFraction']],
                   quickMode= quickMode)
      end.time <- proc.time()
      cat(paste0('Finished build transcript models in ', round((end.time-start.time)[3]/60,1), ' mins', ' \n'))

      se.quant <- bamboo.quantSE(se = se,txdb = NULL, txdbTablesList = txdbTablesList, algo.control = algo.control)

      if(bam.file.index==1){
        seOutput <- se.quant  # create se object
      }else {
        seOutput <- SummarizedExperiment::cbind(seOutput,se.quant)  # combine se object
      }
    }
  }else{
    seList=list()
    for(bam.file.index in seq_along(bam.file)){
      start.time <- proc.time()
      seList[[bam.file.index]]  <- isore.constructReadClasses(bamFile = bam.file[[bam.file.index]],
                   txdbTablesList = txdbTablesList,
                   genomeFA = fa.file,
                   stranded = ir.control[['stranded']],
                   protocol = ir.control[['protocol']],
                   prefix = ir.control[['prefix']],
                   minimumReadSupport= ir.control[['minimumReadSupport']],
                   minimumTxFraction = ir.control[['minimumTxFraction']],
                   quickMode= quickMode)
      end.time <- proc.time()
      cat(paste0('Finished build transcript models in ', round((end.time-start.time)[3]/60,1), ' mins', ' \n'))
    }
  }
  return(seOutput)
}




