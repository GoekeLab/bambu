bamboo <- function(bam.file = NULL, se = NULL, readclass.file = NULL, outputReadClassDir = NULL, txdb = NULL, annotationGrangesList = NULL, genomeSequence = NULL, algo.control = NULL, yieldSize = NULL, ir.control = NULL, extendAnnotations = FALSE){


  #===# Check annotation inputs #===#
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

  ## When SE object from bamboo.quantISORE is provided ##

  if(!is.null(se)){
    return(bamboo.quantSE(se = se, annotationGrangesList = annotationGrangesList, algo.control = algo.control))
  }

  if(!is.null(bam.file) | (!is.null(readclass.file))){
    #===# set default controlling parameters for isoform reconstruction  #===#
    ir.control.default <- list(stranded = FALSE,
                               prefix = '',
                               remove.subsetTx = TRUE, # filter to remove read classes which are a subset of known transcripts.
                               min.readCount = 2,  # minimun read count to consider a read class valid in a sample
                               min.readFractionByGene = 0.05,  ## minimum relative read count per gene, highly expressed genes will have many high read count low relative abundance transcripts that can be filtered
                               min.sampleNumber = 1,  # minimum sample number with minimum read count
                               min.exonDistance = 35,  # minum distance to known transcript to be considered valid as new
                               min.exonOverlap = 10, # minimum number of bases shared with annotation to be assigned to the same gene id
                               prefix='')  ## prefix for new gene Ids (genePrefix.number)
    if(!is.null(ir.control)){
      for(i in names(ir.control)) {
        ir.control.default[[i]] <- ir.control[[i]]
      }
    }
    ir.control <- ir.control.default

    #===# Check whether outputReadClassDir is provided  #===#
    if(!is.null(outputReadClassDir)) {
      if(!dir.exists(outputReadClassDir)) {
        stop("output folder does not exist")
      }
    }

    if(!is.null(readclass.file)){
      if(!all(grepl(".rds", readclass.file))){
        stop("Read class files should be provided in rds format.")
      }
    }

    ## When only directory to readClass SE are provided
    if(is.null(bam.file)){
      return(bamboo.combineQuantify(readclass.file = readclass.file,
                                    annotationGrangesList = annotationGrangesList,
                                    ir.control = ir.control,
                                    algo.control = algo.control,
                                    extendAnnotations = extendAnnotations))
    }else{
      #===# Check genome.fa file  #===#
      if(is.null(genomeSequence)){
        stop("genomeSequence is missing, please provde fasta file or BSgenome name")
      }



      #===# create BamFileList object from character #===#
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

      #===# When there are a lot number of samples, nSample > 10
      if(length(bam.file)>10 &(is.null(outputReadClassDir))){
        outputReadClassDir <- paste0(getwd(),"/tmpReadClassFolder")
        tempfile(pattern = paste0(bam.file.basenames,"_readClassFolder"),
                 tmpdir =  outputReadClassDir,
                 fileext = ".rds")
      }


      ## When bam.file and are provided ##
      if(is.null(outputReadClassDir)){
        return(bamboo.quantISORE(bam.file = bam.file,
                                 algo.control = algo.control,
                                 genomeSequence = genomeSequence,
                                 annotationGrangesList = annotationGrangesList,
                                 ir.control = ir.control,
                                 extendAnnotations = extendAnnotations))
      }else{
        return(bamboo.preprocess(bam.file = bam.file,
                                 algo.control = algo.control,
                                 genomeSequence = genomeSequence,
                                 annotationGrangesList = annotationGrangesList,
                                 ir.control = ir.control,
                                 extendAnnotations = extendAnnotations,
                                 outputReadClassDir = outputReadClassDir))
      }
    }


  }
  stop("At least bam.file, summarizedExperiment output from isore, or directory to saved readClass objects need to be provided.")
}

bamboo.quantDT <- function(dt = dt,algo.control = NULL){
  if(is.null(dt)){
    stop("Input object is missing.")
  }else if(any(!(c('gene_id','tx_id','read_class_id','nobs') %in% colnames(dt)))){
    stop("Columns gene_id, tx_id, read_class_id, nobs, are missing from object.")
  }

  ## check quantification parameters
  default.algo.control <- list(ncore = 1,#parallel::detectCores(),
                               bias_correction = FALSE,
                               maxiter = 10000,
                               convcontrol = 10^(-4))

  if(is.null(algo.control)){
    algo.control <- default.algo.control
  }else{
    if(is.null(algo.control[["ncore"]])){
      algo.control[["ncore"]] <- default.algo.control[["ncore"]]
    }else if(as.numeric(algo.control[["ncore"]])>parallel::detectCores()){
      algo.control[["ncore"]] <- default.algo.control[["ncore"]]
    }
    if(is.null(algo.control[["bias_correction"]])){
      algo.control[["bias_correction"]] <- default.algo.control[["bias_correction"]]
    }
    if(is.null(algo.control[["maxiter"]])){
      algo.control[["maxiter"]] <- default.algo.control[["maxiter"]]
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
                                      bias_correction = algo.control[["bias_correction"]],
                                      maxiter = algo.control[["maxiter"]],
                                      conv.control = algo.control[["convcontrol"]])
  end.time <- proc.time()
  cat(paste0('Finished EM estimation in ', round((end.time-start.time)[3]/60,1), ' mins', ' \n'))


  theta_est <- outList[[1]]
  theta_est[, `:=`(tx_name = txVec[as.numeric(tx_sid)],
                   gene_name = geneVec[gene_sid] )]
  theta_est[,`:=`(tx_sid=NULL, gene_sid=NULL)]
  theta_est <- theta_est[,.(tx_name, estimates)]
  theta_est[,`:=`(CPM = estimates/sum(estimates)*(10^6))]

  b_est <- outList[[2]]
  b_est[, `:=`(gene_name = geneVec[gene_sid], eqClass = eqClassVec[as.numeric(read_class_sid)])]
  b_est[, `:=`(gene_sid = NULL,read_class_sid=NULL)]

  est.list <- list(counts = theta_est, metadata = b_est)
  return(est.list)
}





bamboo.quantSE <- function(se, annotationGrangesList , algo.control = NULL){

  dt <- getEmptyClassFromSE(se, annotationGrangesList)

  est <- bamboo.quantDT(dt,algo.control = algo.control)

  counts <- est$counts
  ColNames <- colnames(se)
  if(length(setdiff(counts$tx_name,names(annotationGrangesList)))>0){
    stop("The provided annotation is incomplete!")
  }
  counts <- counts[data.table(tx_name = names(annotationGrangesList)),  on = 'tx_name']

  counts[is.na(estimates),`:=`(estimates = 0, CPM = 0) ]
  counts <- setDF(counts)
  seOutput <- SummarizedExperiment::SummarizedExperiment(assays = SimpleList(estimates = matrix(counts$estimates,ncol = length(ColNames), dimnames = list(counts$tx_name, ColNames)),
                                                                             normEstimates = matrix(counts$CPM, ncol =  length(ColNames), dimnames = list(counts$tx_name, ColNames))),
                                                         rowRanges = annotationGrangesList[counts$tx_name],
                                                         colData = colData(se),
                                                         metadata = est$metadata) # transcript annotation with read class information

  return(seOutput)
}




bamboo.quantISORE <- function(bam.file = bam.file,annotationGrangesList, genomeSequence = NULL, algo.control = NULL,  ir.control = NULL,  quickMode = FALSE, extendAnnotations=FALSE, outputReadClassDir = NULL){

  bam.file.basenames <- tools::file_path_sans_ext(BiocGenerics::basename(bam.file))
  seOutput = NULL
  if(extendAnnotations==FALSE){
    for(bam.file.index in seq_along(bam.file)){
      start.time <- proc.time()
      readGrgList <- prepareDataFromBam(bam.file[[bam.file.index]])
      se  <- isore.constructReadClasses(readGrgList = readGrgList,
                                        runName =bam.file.basenames[bam.file.index],
                                        annotationGrangesList = annotationGrangesList,
                                        genomeSequence = genomeSequence,
                                        stranded = ir.control[['stranded']],
                                        quickMode= quickMode)
      rm(readGrgList)
      gc()
      end.time <- proc.time()
      cat(paste0('Finished build read classes models in ', round((end.time-start.time)[3]/60,1), ' mins', ' \n'))

      start.time <- proc.time()
      seWithDist <- isore.estimateDistanceToAnnotations(se, annotationGrangesList, min.exonDistance = ir.control[['min.exonDistance']])
      rm(se)
      gc()
      end.time <- proc.time()
      cat(paste0('Finished calculate distance to transcripts in ', round((end.time-start.time)[3]/60,1), ' mins', ' \n'))

      start.time <- proc.time()
      se.quant <- bamboo.quantSE(se = seWithDist, annotationGrangesList = annotationGrangesList, algo.control = algo.control)
      if(bam.file.index==1){
        seOutput <- se.quant  # create se object
      }else {
        seOutput <- SummarizedExperiment::cbind(seOutput,se.quant)  # combine se object
      }
      rm(list=c("se.quant","seWithDist"))
      gc()
      end.time <- proc.time()
      cat(paste0('Finished transcript abundance quantification in ', round((end.time-start.time)[3]/60,1), ' mins', ' \n'))
    }
  }else { # if computation is done in memory in a single session
    seList = list()
    combinedTxCandidates = NULL
    for(bam.file.index in seq_along(bam.file)){  # first loop to reconstruct read classes
      start.time <- proc.time()

      readGrgList <- prepareDataFromBam(bam.file[[bam.file.index]])
      seList[[bam.file.index]]  <- isore.constructReadClasses(readGrgList = readGrgList,
                                                              runName = bam.file.basenames[bam.file.index],
                                                              annotationGrangesList = annotationGrangesList,
                                                              genomeSequence = genomeSequence,
                                                              stranded = ir.control[['stranded']],
                                                              quickMode= quickMode)


      rm(readGrgList)
      gc()
      end.time <- proc.time()
      cat(paste0('Finished build transcript models in ', round((end.time-start.time)[3]/60,1), ' mins', ' \n'))
      combinedTxCandidates <- isore.combineTranscriptCandidates(seList[[bam.file.index]], readClassSeRef = combinedTxCandidates)
    }
    start.time <- proc.time()
    extendedAnnotationGRangesList = isore.extendAnnotations(se=combinedTxCandidates,
                                                            annotationGrangesList=annotationGrangesList,
                                                            remove.subsetTx = ir.control[['remove.subsetTx']],
                                                            min.readCount = ir.control[['min.readCount']],
                                                            min.readFractionByGene = ir.control[['min.sampleNumber']],
                                                            min.sampleNumber = ir.control[['min.sampleNumber']],
                                                            min.exonDistance = ir.control[['min.exonDistance']],
                                                            min.exonOverlap = ir.control[['min.exonOverlap']],
                                                            prefix = ir.control[['prefix']])
    rm(annotationGrangesList)
    gc()
    end.time <- proc.time()
    cat(paste0('Finished extending annotations in ', round((end.time-start.time)[3]/60,1), ' mins', ' \n'))

    for(bam.file.index in seq_along(bam.file)){  # second loop after adding new gene annotations
      seWithDist <- isore.estimateDistanceToAnnotations(seList[[bam.file.index]], extendedAnnotationGRangesList, min.exonDistance = ir.control[['min.exonDistance']])
      se.quant <- bamboo.quantSE(se = seWithDist, annotationGrangesList = extendedAnnotationGRangesList, algo.control = algo.control)
      if(bam.file.index==1){
        seOutput <- se.quant  # create se object
      }else {
        seOutput <- SummarizedExperiment::cbind(seOutput,se.quant)  # combine se object
      }
    }
  }
  return(seOutput)
}


bamboo.preprocess <- function(bam.file = bam.file, annotationGrangesList, genomeSequence = NULL, algo.control = NULL,  ir.control = NULL,  quickMode = FALSE, extendAnnotations=FALSE, outputReadClassDir = NULL){

  bam.file.basenames <- tools::file_path_sans_ext(BiocGenerics::basename(bam.file))
  seOutput = NULL

  readClassFiles <- fs::path(outputReadClassDir,paste0(bam.file.basenames,'_readClassSe'), ext='rds')
  noprint <- lapply(seq_along(bam.file), function(bam.file.index){  # first loop to reconstruct read classes
    start.time <- proc.time()

    readGrgList <- prepareDataFromBam(bam.file[[bam.file.index]])
    se <- isore.constructReadClasses(readGrgList = readGrgList,
                                     runName = bam.file.basenames[bam.file.index],
                                     annotationGrangesList = annotationGrangesList,
                                     genomeSequence = genomeSequence,
                                     stranded = ir.control[['stranded']],
                                     quickMode = quickMode)
    rm(readGrgList)
    gc()

    end.time <- proc.time()
    cat(paste0('Finished build transcript models in ', round((end.time-start.time)[3]/60,1), ' mins', ' \n'))
    if(file.exists(readClassFiles[bam.file.index])){
      warning(paste(readClassFiles[bam.file.index], 'exists, will be overwritten'))
    }
    saveRDS(se, file=readClassFiles[bam.file.index])
    rm(se)
    gc()
  })

  seOutput <- bamboo.combineQuantify(readclass.file = readClassFiles,
                                     annotationGrangesList = annotationGrangesList,
                                     ir.control = ir.control,
                                     algo.control = algo.control,
                                     extendAnnotations = extendAnnotations)

  return(seOutput)
}

#' Combine readClass objects and perform quantification
bamboo.combineQuantify <- function(readclass.file, annotationGrangesList, ir.control, algo.control, extendAnnotations){

  seOutput <- NULL
  if(extendAnnotations==FALSE){
    for(readclass.file.index in seq_along(readclass.file)){  # second loop after adding new gene annotations
      se <- readRDS(file=readclass.file[readclass.file.index])
      seWithDist <- isore.estimateDistanceToAnnotations(se, annotationGrangesList, min.exonDistance = ir.control[['min.exonDistance']])
      se.quant <- bamboo.quantSE(se = seWithDist, annotationGrangesList, algo.control = algo.control) ## NOTE: replace txdbTableList with new annotation table list
      if(readclass.file.index==1){
        seOutput <- se.quant  # create se object
      }else {
        seOutput <- SummarizedExperiment::cbind(seOutput,se.quant)  # combine se object
      }
    }
  }else{
    start.time <- proc.time()
    combinedTxCandidates <- NULL
    for(readclass.file.index in seq_along(readclass.file)){  # second loop after adding new gene annotations
      se <- readRDS(file=readclass.file[readclass.file.index])
      combinedTxCandidates <- isore.combineTranscriptCandidates(se, readClassSeRef = combinedTxCandidates)
      rm(se)
      gc()
    }
    end.time <- proc.time()
    cat(paste0('Finished combining transcript candidates across samples in ', round((end.time-start.time)[3]/60,1), ' mins', ' \n'))



    start.time <- proc.time()
    extendedAnnotationGRangesList = isore.extendAnnotations(se=combinedTxCandidates,
                                                            annotationGrangesList=annotationGrangesList,
                                                            remove.subsetTx = ir.control[['remove.subsetTx']],
                                                            min.readCount = ir.control[['min.readCount']],
                                                            min.readFractionByGene = ir.control[['min.readFractionByGene']],
                                                            min.sampleNumber = ir.control[['min.sampleNumber']],
                                                            min.exonDistance = ir.control[['min.exonDistance']],
                                                            min.exonOverlap = ir.control[['min.exonOverlap']],
                                                            prefix = ir.control[['prefix']])
    rm(list = c("combinedTxCandidates","annotationGrangesList"))
    gc()
    end.time <- proc.time()
    cat(paste0('Finished extending annotations in ', round((end.time-start.time)[3]/60,1), ' mins', ' \n'))


    for(readclass.file.index in seq_along(readclass.file)){  # second loop after adding new gene annotations
      se <- readRDS(file=readclass.file[readclass.file.index])
      seWithDist <- isore.estimateDistanceToAnnotations(se, annotationGrangesList = extendedAnnotationGRangesList, min.exonDistance = ir.control[['min.exonDistance']])
      se.quant <- bamboo.quantSE(se = seWithDist, annotationGrangesList =  extendedAnnotationGRangesList, algo.control = algo.control) ## NOTE: replace txdbTableList with new annotation table list
      if(readclass.file.index==1){
        seOutput <- se.quant  # create se object
      }else {
        seOutput <- SummarizedExperiment::cbind(seOutput,se.quant)  # combine se object
      }
      rm(list = c("seWithDist","se.quant"))
      gc()
    }
  }
  return(seOutput)
}



