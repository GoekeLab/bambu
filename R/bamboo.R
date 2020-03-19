bamboo <- function(bam.file = NULL, se = NULL, dt = NULL,txdb = NULL, annotationGrangesList = NULL, algo.control = NULL, ir.control = NULL, fa.file = NULL, extendAnnotations = FALSE){

  if(!is.null(dt)){
    return(bamboo.quantDT(dt = dt,algo.control = algo.control))
  }

  if(is.null(annotationGrangesList)){
    if(!is.null(txdb)) # txdb object by reading txdb file
    {
      if(class(txdb) != 'TxDb'){
        stop("txdb object is missing.")
      }
      annotationGrangesList <- prepareAnnotations(txdb) ## note: check which annotation tables are reused multiple times, and inlcude only those. Optimise required annotations if possible
    }else{
      stop("txdb object is missing.")
    }
  }

    if(!is.null(se)){
      return(bamboo.quantSE(se = se, annotationGrangesList = annotationGrangesList, algo.control = algo.control))
          }
    if(!is.null(bam.file)){
      return(bamboo.quantISORE(bam.file = bam.file,
                               algo.control = algo.control,
                               fa.file = fa.file,
                               annotationGrangesList = annotationGrangesList,
                               ir.control = ir.control,
                               extendAnnotations = extendAnnotations))
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





bamboo.quantSE <- function(se, annotationGrangesList, algo.control = NULL){
  # if(is.null(txdbTablesList)){
  #   if(!is.null(txdb)) # txdb object by reading txdb file
  #   {
  #     if(class(txdb) != 'TxDb'){
  #       stop("txdb object is missing.")
  #     }
  #     txdbTablesList <- prepareAnnotations(txdb) ## note: check which annotation tables are reused multiple times, and inlcude only those. Optimise required annotations if possible
  #   }else{
  #     stop("txdb object is missing.")
  #   }
  # }
  # annotationGrangesList <- txdbTablesList$exonsByTx
  # mcols(annotationGrangesList) <- dplyr::select(txdbTablesList$txIdToGeneIdTable, TXNAME, GENEID, eqClass)

  dt <- getEmptyClassFromSE(se, annotationGrangesList)

  ## To do:
  ## task1: optional: to implement filtering function
  ## task2: to implement for multiple samples, when multiple samples are provided, run txdbtableslist for one time
  est <- bamboo.quantDT(dt,algo.control = algo.control)
  counts <- est$counts
  ColNames <- colnames(se)
  counts <- merge(counts,data.table(tx_name = names(annotationGrangesList)), all = TRUE,  on = 'tx_name')
  counts[is.na(estimates),`:=`(estimates = 0, CPM = 0) ]
  counts <- setDF(counts)
  seOutput <- SummarizedExperiment::SummarizedExperiment(assays = SimpleList(estimates = matrix(counts$estimates,ncol = length(ColNames), dimnames = list(counts$tx_name, ColNames)),
                                                                       normEstimates = matrix(counts$CPM, ncol =  length(ColNames), dimnames = list(counts$tx_name, ColNames))),
                                   rowRanges = annotationGrangesList[counts$tx_name],
                                   colData = colData(se),
                                   metadata = est$metadata) # transcript annotation with read class information
  return(seOutput)
}



bamboo.quantISORE <- function(bam.file = bam.file, algo.control = NULL, fa.file=NULL, annotationGrangesList , ir.control = NULL, yieldSize = NULL, quickMode = FALSE, extendAnnotations=FALSE, outputReadClassToFolder = NULL){

  if(is.null(fa.file)){
    stop("GenomeFA file is missing.")
  }else if(class(fa.file) != 'FaFile'){
    if(!grepl('.fa',fa.file)){
      stop("GenomeFA file is missing.")
    }else{
      fa.file <- Rsamtools::FaFile(fa.file)
    }
  }

  if(!is.null(outputReadClassToFolder)) {
    if(!dir.exists(outputReadClassToFolder)) {
      stop("output folder does not exist")
    }
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

  # if(is.null(txdbTablesList)){
  #   if(!is.null(txdb)) # txdb object by reading txdb file
  #   {
  #     if(class(txdb) != 'TxDb'){
  #       stop("txdb object is missing.")
  #     }
  #     txdbTablesList <- prepareAnnotations(txdb) ## note: check which annotation tables are reused multiple times, and inlcude only those. Optimise required annotations if possible
  #   }else{
  #     stop("txdb object is missing.")
  #   }
  # }

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

  ## remove
  #annotationGrangesList <- txdbTablesList$exonsByTx
 # mcols(annotationGrangesList) <- dplyr::select(txdbTablesList$txIdToGeneIdTable, TXNAME, GENEID, eqClass)

  bam.file.basenames <- tools::file_path_sans_ext(BiocGenerics::basename(bam.file))
  seOutput = NULL
  if(extendAnnotations==FALSE) {
    for(bam.file.index in seq_along(bam.file)){
      start.time <- proc.time()
      readGrgList <- isore.preprocessBam(bam.file[[bam.file.index]])
      se  <- isore.constructReadClasses(readGrgList = readGrgList,
                                                              runName =bam.file.basenames[bam.file.index],
                                                              annotationGrangesList = annotationGrangesList,
                                                              genomeFA = fa.file,
                                                              stranded = ir.control[['stranded']],
                                                              protocol = ir.control[['protocol']],
                                                              prefix = ir.control[['prefix']],
                                                              minimumReadSupport= ir.control[['minimumReadSupport']],
                                                              minimumTxFraction = ir.control[['minimumTxFraction']],
                                                              quickMode= quickMode)
      end.time <- proc.time()
      cat(paste0('Finished build read classes models in ', round((end.time-start.time)[3]/60,1), ' mins', ' \n'))
      start.time <- proc.time()

      seWithDist <- isore.estimateDistanceToAnnotations(se, annotationGrangesList, stranded=ir.control[['stranded']])
      end.time <- proc.time()
      cat(paste0('Finished calculate distance to transcripts in ', round((end.time-start.time)[3]/60,1), ' mins', ' \n'))
      start.time <- proc.time()

      se.quant <- bamboo.quantSE(se = seWithDist, annotationGrangesList = annotationGrangesList, algo.control = algo.control)
      if(bam.file.index==1){
        seOutput <- se.quant  # create se object
      }else {
        seOutput <- SummarizedExperiment::cbind(seOutput,se.quant)  # combine se object
      }
      end.time <- proc.time()
      cat(paste0('Finished transcript abundance quantification in ', round((end.time-start.time)[3]/60,1), ' mins', ' \n'))
    }
  }else if (!is.null(outputReadClassToFolder)){  # if data is written to output directory
    readClassFiles <- fs::path(outputReadClassToFolder,paste0(bam.file.basenames,'_readClassSe'), ext='rds')
    combinedTxCandidates <- NULL
    for(bam.file.index in seq_along(bam.file)){  # first loop to reconstruct read classes
      start.time <- proc.time()
      readGrgList <- isore.preprocessBam(bam.file[[bam.file.index]])
      se <- isore.constructReadClasses(readGrgList = readGrgList,
                                       runName = bam.file.basenames[bam.file.index],
                                       annotationGrangesList = annotationGrangesList,
                                       genomeFA = fa.file,
                                       stranded = ir.control[['stranded']],
                                       protocol = ir.control[['protocol']],
                                       prefix = ir.control[['prefix']],
                                       minimumReadSupport= ir.control[['minimumReadSupport']],
                                       minimumTxFraction = ir.control[['minimumTxFraction']],
                                       quickMode = quickMode)
      end.time <- proc.time()
      rm(readGrgList)
      cat(paste0('Finished build transcript models in ', round((end.time-start.time)[3]/60,1), ' mins', ' \n'))
      #readClassFile <-  fs::path(outputReadClassToFolder,paste0(bam.file.basenames[bam.file.index], '_readClassSe'), ext = 'rds')
      if(file.exists(readClassFiles[bam.file.index])){
        warning(paste(readClassFiles[bam.file.index], 'exists, will be overwritten'))
      }
      saveRDS(se, file=readClassFiles[bam.file.index])
      combinedTxCandidates <- isore.combineTranscriptCandidates(se, readClassSeRef = combinedTxCandidates)
      rm(se)
      gc()
    }

    # add new transcripts
   # annotationGrangesList <- txdbTablesList$exonsByTx
   # mcols(annotationGrangesList) <- dplyr::select(txdbTablesList$txIdToGeneIdTable, TXNAME, GENEID, referenceTXNAME, eqClass, minEqClassSize)
    extendedAnnotationGRangesList = isore.extendAnnotations(se=combinedTxCandidates, annotationGrangesList=annotationGrangesList) ## missing

    for(bam.file.index in seq_along(bam.file)){  # second loop after adding new gene annotations
      se <- readRDS(file=readClassFiles[bam.file.index])
      seWithDist <- isore.estimateDistanceToAnnotations(se, extendedAnnotationGRangesList, stranded=ir.control[['stranded']])
      se.quant <- bamboo.quantSE(se = seWithDist, annotationGrangesList=annotationGrangesList, algo.control = algo.control) ## NOTE: replace txdbTableList with new annotation table list
      if(bam.file.index==1){
        seOutput <- se.quant  # create se object
      }else {
        seOutput <- SummarizedExperiment::cbind(seOutput,se.quant)  # combine se object
      }
    }

  }else { # if computation is done in memory in a single session
    seList = list()
    combinedTxCandidates = NULL
    for(bam.file.index in seq_along(bam.file)){  # first loop to reconstruct read classes
      start.time <- proc.time()
      readGrgList <- isore.preprocessBam(bam.file[[bam.file.index]])
      seList[[bam.file.index]]  <- isore.constructReadClasses(readGrgList = readGrgList,
                                                              runName = bam.file.basenames[bam.file.index],
                                                              annotationGrangesList = annotationGrangesList,
                                                              genomeFA = fa.file,
                                                              stranded = ir.control[['stranded']],
                                                              protocol = ir.control[['protocol']],
                                                              prefix = ir.control[['prefix']],
                                                              minimumReadSupport= ir.control[['minimumReadSupport']],
                                                              minimumTxFraction = ir.control[['minimumTxFraction']],
                                                              quickMode= quickMode)
      end.time <- proc.time()
      cat(paste0('Finished build transcript models in ', round((end.time-start.time)[3]/60,1), ' mins', ' \n'))
      combinedTxCandidates <- isore.combineTranscriptCandidates(seList[[bam.file.index]], readClassSeRef = combinedTxCandidates)
    }

    #annotationGrangesList <- txdbTablesList$exonsByTx
   # mcols(annotationGrangesList) <- dplyr::select(txdbTablesList$txIdToGeneIdTable, TXNAME, GENEID)

    extendedAnnotationGRangesList = isore.extendAnnotations(se=combinedTxCandidates, annotationGrangesList=annotationGrangesList) ## missing

    for(bam.file.index in seq_along(bam.file)){  # second loop after adding new gene annotations


      seWithDist <- isore.estimateDistanceToAnnotations(seList[[bam.file.index]], extendedAnnotationGRangesList, stranded=ir.control[['stranded']]) ## NOTE: replace txdbTableList with new annotation table list
      se.quant <- bamboo.quantSE(se = seWithDist, annotationGrangesList = annotationGrangesList, algo.control = algo.control) ## NOTE: replace txdbTableList with new annotation table list
      if(bam.file.index==1){
        seOutput <- se.quant  # create se object
      }else {
        seOutput <- SummarizedExperiment::cbind(seOutput,se.quant)  # combine se object
      }
    }
  }
  return(seOutput)
}




