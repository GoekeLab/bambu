#' Main function
#' @title long read isoform reconstruction and quantification
#' @description This function takes bam file of genomic alignments and performs isoform recontruction and gene and transcript expression quantification.
#' It also allows saving of read class files of alignments, extending provided annotations, and quantification based on extended annotations.
#' When multiple samples are provided, extended annotations will be combined across samples to allow comparison.
#' @param reads A string or a vector of strings specifying the paths of bam files for genomic alignments, or a \code{\link{BamFile}} object or a \code{\link{BamFileList}}  object (see \code{\link{Rsamtools}}).
#' @param readclass.file A string or a vector of strings specifying the read class files that are saved during previous run of \code{\link{bambu}}.
#' @param outputReadClassDir A string variable specifying the path to where read class files will be saved.
#' @param annotations A \code{\link{TxDb}} object or A GRangesList object obtained by \code{\link{prepareAnnotations}} or \code{\link{prepareAnnotationsFromGTF}}.
#' @param extendAnnotations A logical variable indicating whether annotations are to be extended for quantification.
#' @param genomeSequence A fasta file or a BSGenome object.
#' @param algo.control A list of controlling parameters for quantification algorithm estimation process:
#' \itemize{
#'   \item ncore specifying number of cores used when parallel processing is used, defaults to 1.
#'   \item maxiter specifying maximum number of run interations, defaults to 10000.
#'   \item bias_correction specifying whether to correct for bias, defaults to FALSE.
#'   \item convcontrol specifying the covergence trheshold control, defaults to 0.0001.
#' }
#' @param yieldSize see \code{\link{Rsamtools}}.
#' @param ir.control A list of controlling parameters for isoform reconstruction process:
#' \itemize{
#'   \item whether stranded, defaults to FALSE
#'   \item prefix specifying prefix for new gene Ids (genePrefix.number), defaults to empty
#'   \item remove.subsetTx indicating whether filter to remove read classes which are a subset of known transcripts(), defaults to TRUE
#'   \item min.readCount specifying minimun read count to consider a read class valid in a sample, defaults to 2
#'   \item min.readFractionByGene specifying minimum relative read count per gene, highly expressed genes will have many high read count low relative abundance transcripts that can be filtered, defaults to 0.05
#'   \item min.sampleNumber specifying minimum sample number with minimum read count, defaults to 1
#'   \item min.exonDistance specifying minum distance to known transcript to be considered valid as new, defaults to 35
#'   \item min.exonOverlap specifying minimum number of bases shared with annotation to be assigned to the same gene id, defaults 10 base pairs
#' }
#' @param verbose A logical variable indicating whether processing messages will be printed.
#' @details
#' @return A list of two SummarizedExperiment object for transcript expression and gene expression.
#' @examples
#' \dontrun{
#' ## =====================
#' ## More stringent new gene/isoform discovery: new isoforms are identified with at least 5 read count in 1 sample
#' ## Increase EM convergence threshold to 10^(-6)
#' seOutput <- bambu(reads, annotationGrangesList,
#' genomeSequence, ir.control = list(min.readCount=5),
#' algo.control = list(convcontrol = 10^(-6))
#' }
#' @export
bambu <- function(reads = NULL, readclass.file = NULL, outputReadClassDir = NULL, annotations = NULL, genomeSequence = NULL, algo.control = NULL, yieldSize = NULL, ir.control = NULL, extendAnnotations = TRUE, verbose = FALSE){


  #===# Check annotation inputs #===#
  if(!is.null(annotations)){
      if(class(annotations) == 'TxDb'){
        txdb <- annotations
        annotationGrangesList <- prepareAnnotations(txdb) 
      }else if(class(annotations) == "CompressedGRangesList"){
        annotationGrangesList <- annotations
        ## check if annotationGrangesList is as expected
        if(!all(c("TXNAME","GENEID","eqClass") %in% colnames(mcols(annotationGrangesList)))){
         stop("The annotations is not properly prepared.\nPlease prepareAnnnotations using prepareAnnotations or prepareAnnotationsFromGTF functions.")
        }
      }else{
        stop("The annotations is not a GRangesList object.")
      }
        }else{
      stop("Annotations is missing.")
    }


  ## When SE object from bambu.quantISORE is provided ##
  if(!is.null(reads) | (!is.null(readclass.file))){
    #===# set default controlling parameters for isoform reconstruction  #===#
    ir.control.default <- list(stranded = FALSE,
                               remove.subsetTx = TRUE, #
                               min.readCount = 2,  #
                               min.readFractionByGene = 0.05,  ##
                               min.sampleNumber = 1,  #
                               min.exonDistance = 35,  #
                               min.exonOverlap = 10, #
                               prefix='')  ##
    if(!is.null(ir.control)){
      for(i in names(ir.control)) {
        ir.control.default[[i]] <- ir.control[[i]]
      }
    }
    ir.control <- ir.control.default

    #===# Check whether provided outputReadClassDir exists  #===#
    if(!is.null(outputReadClassDir)) {
      if(!dir.exists(outputReadClassDir)) {
        stop("output folder does not exist")
      }
    }

    #===# Check whether provided readclass files are all in rds format #===#
    if(!is.null(readclass.file)){
      if(!all(grepl(".rds", readclass.file))){
        stop("Read class files should be provided in rds format.")
      }
    }

    ## When only directory to readClass SE are provided
    if(is.null(reads)){
      return(bambu.combineQuantify(readclass.file = readclass.file,
                                    annotationGrangesList = annotationGrangesList,
                                    ir.control = ir.control,
                                    algo.control = algo.control,
                                    extendAnnotations = extendAnnotations,
                                    verbose = verbose))
    }else{
      bam.file <- reads

      #===# Check genome.fa file  #===#
      if(is.null(genomeSequence)){
        stop("genomeSequence is missing, please provde fasta file or BSgenome name")
      }



      #===# create BamFileList object from character #===#
      if(class(bam.file)=='BamFile') {
        if(!is.null(yieldSize)) {
          Rsamtools::yieldSize(bam.file) <- yieldSize
        } else {
          yieldSize <- Rsamtools::yieldSize(bam.file)
        }
        bam.file<- Rsamtools::BamFileList(bam.file)
      }else if(class(bam.file)=='BamFileList') {
        if(!is.null(yieldSize)) {
          Rsamtools::yieldSize(bam.file) <- yieldSize
        } else {
          yieldSize <- min(Rsamtools::yieldSize(bam.file))
        }
      }else if(any(!grepl('\\.bam$',bam.file))){
        stop("Bam file is missing from arguments.")
      }else{
        if(is.null(yieldSize)) {
          yieldSize <- NA
        }

        bam.file<- Rsamtools::BamFileList(bam.file, yieldSize = yieldSize)
      }

      #===# When there are a lot number of samples, nSample > 10
      if(length(bam.file)>10 &(is.null(outputReadClassDir))){
        outputReadClassDir <- paste0(getwd(),"/tmpReadClassFolder")
        warning(paste0("There are more than 10 samples, read class files will be saved to ",outputReadClassDir, " for more efficient process!"))
        tempfile(pattern = paste0(bam.file.basenames,"_readClassFolder"),
                 tmpdir =  outputReadClassDir,
                 fileext = ".rds")
      }


      ## When bam files are provided ##
      if(is.null(outputReadClassDir)){
        return(bambu.quantISORE(bam.file= bam.file,
                                 algo.control = algo.control,
                                 genomeSequence = genomeSequence,
                                 annotationGrangesList = annotationGrangesList,
                                 ir.control = ir.control,
                                 extendAnnotations = extendAnnotations,
                                 verbose = verbose))
      }else{
        return(bambu.preprocess(bam.file= bam.file,
                                 algo.control = algo.control,
                                 genomeSequence = genomeSequence,
                                 annotationGrangesList = annotationGrangesList,
                                 ir.control = ir.control,
                                 extendAnnotations = extendAnnotations,
                                 outputReadClassDir = outputReadClassDir,
                                 verbose = verbose))
      }
    }


  }
  stop("At least bam file, summarizedExperiment output from isore, or directory to saved readClass objects need to be provided.")
}


#' Process data.table object
#' @param dt A data.table object
#' @inheritParams bambu
#' @noRd
bambu.quantDT <- function(dt = dt,algo.control = NULL, verbose = FALSE){
  if(is.null(dt)){
    stop("Input object is missing.")
  }else if(any(!(c('gene_id','tx_id','read_class_id','nobs') %in% colnames(dt)))){
    stop("Columns gene_id, tx_id, read_class_id, nobs, are missing from object.")
  }

  ## check quantification parameters
  algo.control.default <- list(ncore = 1,#parallel::detectCores(),
                               bias_correction = FALSE,
                               maxiter = 10000,
                               convcontrol = 10^(-4))

  if(!is.null(algo.control)){
    for(i in names(algo.control)) {
      algo.control.default[[i]] <- algo.control[[i]]
    }
  }
  algo.control <- algo.control.default


  ##----step2: match to simple numbers to increase claculation efficiency
  geneVec <- unique(dt$gene_id)
  txVec <- unique(dt$tx_id)
  readclassVec <- unique(dt$read_class_id)
  dt <- as.data.table(dt)
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
  if(verbose)   message('Finished EM estimation in ', round((end.time-start.time)[3]/60,1), ' mins.')


  theta_est <- outList[[1]]
  theta_est[, `:=`(tx_name = txVec[as.numeric(tx_sid)],
                   gene_name = geneVec[gene_sid] )]
  theta_est[,`:=`(tx_sid=NULL, gene_sid=NULL)]
  theta_est <- theta_est[,.(tx_name, estimates)]

  theta_est[,`:=`(CPM = estimates/sum(estimates)*(10^6))]



  b_est <- outList[[2]]
  b_est[, `:=`(gene_name = geneVec[gene_sid], eqClass = eqClassVec[as.numeric(read_class_sid)])]
  b_est[, `:=`(gene_sid = NULL,read_class_sid=NULL)]

  est.list <- list(counts = theta_est,
                   metadata = b_est)
  return(est.list)
}




#' Process SummarizedExperiment object
#' @param se A summarizedExperiment object
#' @inheritParams bambu
#' @noRd
bambu.quantSE <- function(se, annotationGrangesList , algo.control = NULL, verbose = FALSE){

  dt <- getEmptyClassFromSE(se, annotationGrangesList)

  est <- bambu.quantDT(dt,algo.control = algo.control, verbose = verbose)

  counts <- est$counts
  ColNames <- colnames(se)
  if(length(setdiff(counts$tx_name,names(annotationGrangesList)))>0){
    stop("The provided annotation is incomplete!")
  }
  counts <- counts[data.table(tx_name = names(annotationGrangesList)),  on = 'tx_name']

  counts[is.na(estimates),`:=`(estimates = 0, CPM = 0) ]
  counts <- setDF(counts)
  seOutput <- SummarizedExperiment::SummarizedExperiment(assays = SimpleList(counts = matrix(counts$estimates,ncol = length(ColNames), dimnames = list(counts$tx_name, ColNames)),
                                                                             CPM = matrix(counts$CPM, ncol =  length(ColNames), dimnames = list(counts$tx_name, ColNames))),
                                                         rowRanges = annotationGrangesList[counts$tx_name],
                                                         colData = colData(se),
                                                         metadata = est$metadata) # transcript annotation with read class information

  return(seOutput)
}


#' Process bam files without saving to folders.
#' @inheritParams bambu
#' @noRd
bambu.quantISORE <- function(bam.file = bam.file, annotationGrangesList, genomeSequence = NULL, algo.control = NULL,  ir.control = NULL, extendAnnotations=FALSE, outputReadClassDir = NULL, verbose = FALSE){

  bam.file.basenames <- tools::file_path_sans_ext(BiocGenerics::basename(bam.file))
  seOutput = NULL
  if(extendAnnotations==FALSE){
    for(bam.file.index in seq_along(bam.file)){
      start.time <- proc.time()
      readGrgList <- prepareDataFromBam(bam.file[[bam.file.index]], verbose = verbose)
      seqlevelsStyle(readGrgList) <- seqlevelsStyle(annotationGrangesList)[1]
      
      se  <- isore.constructReadClasses(readGrgList = readGrgList,
                                        runName =bam.file.basenames[bam.file.index],
                                        annotationGrangesList = annotationGrangesList,
                                        genomeSequence = genomeSequence,
                                        stranded = ir.control[['stranded']],
                                        verbose = verbose)
      rm(readGrgList)
      gc(verbose = FALSE)
      end.time <- proc.time()
      if(verbose)   message('Finished build read classes models in ', round((end.time-start.time)[3]/60,1), ' mins.')

      start.time <- proc.time()
      seWithDist <- isore.estimateDistanceToAnnotations(se, annotationGrangesList, min.exonDistance = ir.control[['min.exonDistance']], verbose = verbose)
      rm(se)
      gc(verbose = FALSE)
      end.time <- proc.time()
      if(verbose)   message('Finished calculate distance to transcripts in ', round((end.time-start.time)[3]/60,1), ' mins.')

      start.time <- proc.time()
      se.quant <- bambu.quantSE(se = seWithDist, annotationGrangesList = annotationGrangesList, algo.control = algo.control, verbose = verbose)
      if(bam.file.index==1){
        seOutput <- se.quant  # create se object
      }else {
        seOutput <- SummarizedExperiment::cbind(seOutput,se.quant)  # combine se object
      }
      rm(list=c("se.quant","seWithDist"))
      gc(verbose = FALSE)
      end.time <- proc.time()
      if(verbose)   message('Finished transcript abundance quantification in ', round((end.time-start.time)[3]/60,1), ' mins.')
    }


  }else { # if computation is done in memory in a single session
    seList = list()
    combinedTxCandidates = NULL
    for(bam.file.index in seq_along(bam.file)){  # first loop to reconstruct read classes
      start.time <- proc.time()

      readGrgList <- prepareDataFromBam(bam.file[[bam.file.index]], verbose = verbose)
      seqlevelsStyle(readGrgList) <- seqlevelsStyle(annotationGrangesList)[1]
      
      seList[[bam.file.index]]  <- isore.constructReadClasses(readGrgList = readGrgList,
                                                              runName = bam.file.basenames[bam.file.index],
                                                              annotationGrangesList = annotationGrangesList,
                                                              genomeSequence = genomeSequence,
                                                              stranded = ir.control[['stranded']],
                                                              verbose = verbose)


      rm(readGrgList)
      gc(verbose = FALSE)
      end.time <- proc.time()
      if(verbose)   message('Finished build transcript models in ', round((end.time-start.time)[3]/60,1), ' mins.')
      combinedTxCandidates <- isore.combineTranscriptCandidates(seList[[bam.file.index]], readClassSeRef = combinedTxCandidates, verbose = verbose)
    }
    start.time <- proc.time()
    extendedAnnotationGRangesList = isore.extendAnnotations(se=combinedTxCandidates,
                                                            annotationGrangesList=annotationGrangesList,
                                                            remove.subsetTx = ir.control[['remove.subsetTx']],
                                                            min.readCount = ir.control[['min.readCount']],
                                                            min.readFractionByGene = ir.control[['min.readFractionByGene']],
                                                            min.sampleNumber = ir.control[['min.sampleNumber']],
                                                            min.exonDistance = ir.control[['min.exonDistance']],
                                                            min.exonOverlap = ir.control[['min.exonOverlap']],
                                                            prefix = ir.control[['prefix']],
                                                            verbose = verbose)
    rm(annotationGrangesList)
    gc(verbose = FALSE)
    end.time <- proc.time()
    if(verbose)   message('Finished extending annotations in ', round((end.time-start.time)[3]/60,1), ' mins.')

    for(bam.file.index in seq_along(bam.file)){  # second loop after adding new gene annotations
      start.time <- proc.time()
      seWithDist <- isore.estimateDistanceToAnnotations(seList[[bam.file.index]], extendedAnnotationGRangesList, min.exonDistance = ir.control[['min.exonDistance']], verbose = verbose)
      end.time <- proc.time()
      if(verbose)   message('Finished calculate distance to transcripts in ', round((end.time-start.time)[3]/60,1), ' mins.')

      se.quant <- bambu.quantSE(se = seWithDist, annotationGrangesList = extendedAnnotationGRangesList, algo.control = algo.control, verbose = verbose)
      if(bam.file.index==1){
        seOutput <- se.quant  # create se object
      }else {
        seOutput <- SummarizedExperiment::cbind(seOutput,se.quant)  # combine se object
      }
      rm(list = c("seWithDist","se.quant"))
      gc(verbose = FALSE)
    }

  }

  return(seOutput)
}

#' Preprocess bam files and save read class files
#' @inheritParams bambu
#' @noRd
bambu.preprocess <- function(bam.file= bam.file, annotationGrangesList, genomeSequence = NULL, algo.control = NULL,  ir.control = NULL, extendAnnotations=FALSE, outputReadClassDir = NULL, verbose = FALSE){

  bam.file.basenames <- tools::file_path_sans_ext(BiocGenerics::basename(bam.file))
  seOutput = NULL

  readClassFiles <- fs::path(outputReadClassDir,paste0(bam.file.basenames,'_readClassSe'), ext='rds')
  noprint <- lapply(seq_along(bam.file), function(bam.file.index){  # first loop to reconstruct read classes
    start.time <- proc.time()

    readGrgList <- prepareDataFromBam(bam.file[[bam.file.index]], verbose = verbose)
    seqlevelsStyle(readGrgList) <- seqlevelsStyle(annotationGrangesList)[1]
    
    se <- isore.constructReadClasses(readGrgList = readGrgList,
                                     runName = bam.file.basenames[bam.file.index],
                                     annotationGrangesList = annotationGrangesList,
                                     genomeSequence = genomeSequence,
                                     stranded = ir.control[['stranded']],
                                     verbose = verbose)
    rm(readGrgList)
    gc(verbose = FALSE)

    end.time <- proc.time()
    if(verbose)   message('Finished build transcript models in ', round((end.time-start.time)[3]/60,1), ' mins.')
    if(file.exists(readClassFiles[bam.file.index])){
      warning(paste(readClassFiles[bam.file.index], 'exists, will be overwritten'))
    }
    saveRDS(se, file=readClassFiles[bam.file.index])
    rm(se)
    gc(verbose = FALSE)
  })

  seOutput <- bambu.combineQuantify(readclass.file = readClassFiles,
                                     annotationGrangesList = annotationGrangesList,
                                     ir.control = ir.control,
                                     algo.control = algo.control,
                                     extendAnnotations = extendAnnotations,
                                     verbose = verbose)

  return(seOutput)
}

#' Combine readClass objects and perform quantification
#' @inheritParams bambu
#' @noRd
bambu.combineQuantify <- function(readclass.file, annotationGrangesList, ir.control, algo.control, extendAnnotations, verbose = FALSE){

  seOutput <- NULL
  if(extendAnnotations==FALSE){
    for(readclass.file.index in seq_along(readclass.file)){  # second loop after adding new gene annotations

      start.time <- proc.time()
      se <- readRDS(file=readclass.file[readclass.file.index])
      seqlevelsStyle(se) <- seqlevelsStyle(annotationGrangesList)[1]
      seWithDist <- isore.estimateDistanceToAnnotations(se, annotationGrangesList, min.exonDistance = ir.control[['min.exonDistance']], verbose = verbose)
      end.time <- proc.time()
      if(verbose)   message('Finished calculate distance to transcripts in ', round((end.time-start.time)[3]/60,1), ' mins.')

      se.quant <- bambu.quantSE(se = seWithDist, annotationGrangesList, algo.control = algo.control, verbose = verbose) ## NOTE: replace txdbTableList with new annotation table list
      if(readclass.file.index==1){
        seOutput <- se.quant  # create se object
      }else {
        seOutput <- SummarizedExperiment::cbind(seOutput,se.quant)  # combine se object
      }
      rm(list = c("seWithDist","se.quant"))
      gc(verbose = FALSE)
    }

  }else{
    start.time <- proc.time()
    combinedTxCandidates <- NULL
    for(readclass.file.index in seq_along(readclass.file)){  # second loop after adding new gene annotations
      se <- readRDS(file=readclass.file[readclass.file.index])
      seqlevelsStyle(se) <- seqlevelsStyle(annotationGrangesList)[1]
      combinedTxCandidates <- isore.combineTranscriptCandidates(se, readClassSeRef = combinedTxCandidates, verbose = verbose)
      rm(se)
      gc(verbose = FALSE)
    }
    end.time <- proc.time()
    if(verbose)   message('Finished combining transcript candidates across samples in ', round((end.time-start.time)[3]/60,1), ' mins.')



    start.time <- proc.time()
    extendedAnnotationGRangesList = isore.extendAnnotations(se=combinedTxCandidates,
                                                            annotationGrangesList=annotationGrangesList,
                                                            remove.subsetTx = ir.control[['remove.subsetTx']],
                                                            min.readCount = ir.control[['min.readCount']],
                                                            min.readFractionByGene = ir.control[['min.readFractionByGene']],
                                                            min.sampleNumber = ir.control[['min.sampleNumber']],
                                                            min.exonDistance = ir.control[['min.exonDistance']],
                                                            min.exonOverlap = ir.control[['min.exonOverlap']],
                                                            prefix = ir.control[['prefix']],
                                                            verbose = verbose)
    rm(list = c("combinedTxCandidates","annotationGrangesList"))
    gc(verbose = FALSE)
    end.time <- proc.time()
    if(verbose)   message('Finished extending annotations in ', round((end.time-start.time)[3]/60,1), ' mins.')


    for(readclass.file.index in seq_along(readclass.file)){  # second loop after adding new gene annotations
      start.time <- proc.time()
      se <- readRDS(file=readclass.file[readclass.file.index])
      seqlevelsStyle(se) <- seqlevelsStyle(annotationGrangesList)[1]
      
      seWithDist <- isore.estimateDistanceToAnnotations(se, annotationGrangesList = extendedAnnotationGRangesList, min.exonDistance = ir.control[['min.exonDistance']], verbose = verbose)
      end.time <- proc.time()
      if(verbose)   message('Finished calculate distance to transcripts in ', round((end.time-start.time)[3]/60,1), ' mins.')

      se.quant <- bambu.quantSE(se = seWithDist, annotationGrangesList =  extendedAnnotationGRangesList, algo.control = algo.control, verbose = verbose) ## NOTE: replace txdbTableList with new annotation table list
      if(readclass.file.index==1){
        seOutput <- se.quant  # create se object
      }else {
        seOutput <- SummarizedExperiment::cbind(seOutput,se.quant)  # combine se object
      }
      rm(list = c("seWithDist","se.quant"))
      gc(verbose = FALSE)
    }


  }
  return(seOutput)
}


