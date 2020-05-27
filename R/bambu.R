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
bambu <- function(reads = NULL, readclass.file = NULL, outputReadClassDir = NULL, annotations = NULL, genomeSequence = NULL, algo.control = NULL, yieldSize = NULL, ir.control = NULL, extendAnnotations = TRUE, stranded = FALSE, ncore = 1, verbose = FALSE){


  #===# Check annotation inputs #===#
  if(!is.null(annotations)){
      if(class(annotations) == 'TxDb'){
        annotations <- prepareAnnotations(annotations)
      }else if(class(annotations) == "CompressedGRangesList"){
        ## check if annotations is as expected
        if(!all(c("TXNAME","GENEID","eqClass") %in% colnames(mcols(annotations)))){
         stop("The annotations is not properly prepared.\nPlease prepareAnnnotations using prepareAnnotations or prepareAnnotationsFromGTF functions.")
        }
      }else{
        stop("The annotations is not a GRangesList object.")
      }
        }else{
      stop("Annotations is missing.")
    }


  ## When SE object from bambu.quantISORE is provided ##
  if(!is.null(reads) & (!is.null(readclass.file))){
    
    stop("At least bam file or path to readClass file needs to be provided.")
  }
  
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

    ## check quantification parameters
    algo.control.default <- list(bias_correction = FALSE,
                                 maxiter = 10000,
                                 convcontrol = 10^(-4))
    
    if(!is.null(algo.control)){
      for(i in names(algo.control)) {
        algo.control.default[[i]] <- algo.control[[i]]
      }
    }
    algo.control <- algo.control.default
    
    rm.readClassSe <- FALSE  # indicator to remove temporary read class files
    
    bpParameters <- BiocParallel::bpparam()
    #===# set parallel options: If more CPUs than samples available, use parallel computing on each sample, otherwise use parallel to distribute samples (more efficient)
    if(length(reads)<=(0.5*ncore)) {
      bpParameters$workers <- 1
    } else {
      bpParameters$workers <- ncore
      ncore <- 1
    }
    
    if(!is.null(reads)){  # calculate readClass objects
      
      #===# create BamFileList object from character #===#
      if(class(reads)=='BamFile') {
        if(!is.null(yieldSize)) {
          Rsamtools::yieldSize(reads) <- yieldSize
        } else {
          yieldSize <- Rsamtools::yieldSize(reads)
        }
        reads<- Rsamtools::BamFileList(reads)
        names(reads) <- tools::file_path_sans_ext(BiocGenerics::basename(reads))
      }else if(class(reads)=='BamFileList') {
        if(!is.null(yieldSize)) {
          Rsamtools::yieldSize(reads) <- yieldSize
        } else {
          yieldSize <- min(Rsamtools::yieldSize(reads))
        }
      }else if(any(!grepl('\\.bam$',reads))){
        stop("Bam file is missing from arguments.")
      }else{
        if(is.null(yieldSize)) {
          yieldSize <- NA
        }
        reads<- Rsamtools::BamFileList(reads, yieldSize = yieldSize)
        names(reads) <- tools::file_path_sans_ext(BiocGenerics::basename(reads))
      }
      
      #===# When more than 10 samples are provided, files will be written to a temporary directory

      if(length(reads)>10 &(is.null(outputReadClassDir))){
        outputReadClassDir <- tempdir()
        warning(paste0("There are more than 10 samples, read class files will be temporarily saved to ",outputReadClassDir, " for more efficient processing"))
        rm.readClassSe <- TRUE # remove temporary read class files from system
      }
      
   
        
      
      
      readClassList <- BiocParallel::bplapply(names(reads), function(bamFileName){
                              bambu.constructReadClass(
                              bam.file= reads[bamFileName],
                              outputReadClassDir=outputReadClassDir,
                              genomeSequence = genomeSequence,
                              annotations = annotations,
                              stranded=stranded,
                              ncore = ncore,
                              verbose = verbose)}, 
                              BPPARAM=bpParameters)
    } else {
      readClassList <- readclass.file
    }
    

    if(extendAnnotations) {
      annotations <- bambu.extendAnnotations(readClassList, annotations, ir.control, verbose=verbose)
      gc(verbose = FALSE)
    }
    countsSe <- BiocParallel::bplapply(readClassList,
                                       bambu.quantify,
                                       annotations=annotations,
                                       min.exonDistance= ir.control[['min.exonDistance']],
                                       algo.control = algo.control,
                                       ncore = ncore,
                                       verbose = verbose, 
                                       BPPARAM=bpParameters)
    countsSe <- do.call(SummarizedExperiment::cbind, countsSe)
    rowRanges(countsSe) <- annotations
    
    #===# Clean up temp directory
      if(rm.readClassSe){
        file.remove(unlist(readClassList))
      }
     return(countsSe)
}

#' Extend annotations
#' @inheritParams bambu
#' @noRd
bambu.extendAnnotations <- function(readClassList, annotations, ir.control, verbose = FALSE){
  combinedTxCandidates = NULL
  for(readClassIndex in seq_along(readClassList)){
    readClass <- readClassList[[readClassIndex]]
    if(is.character(readClass)){
      readClass <- readRDS(file=readClass)
      seqlevelsStyle(readClass) <- seqlevelsStyle(annotations)[1]
    }
    combinedTxCandidates <- isore.combineTranscriptCandidates(readClass, readClassSeRef = combinedTxCandidates, verbose = verbose)
  }
  annotations <- isore.extendAnnotations(se=combinedTxCandidates,
                                         annotationGrangesList=annotations,
                                         remove.subsetTx = ir.control[['remove.subsetTx']],
                                         min.readCount = ir.control[['min.readCount']],
                                         min.readFractionByGene = ir.control[['min.readFractionByGene']],
                                         min.sampleNumber = ir.control[['min.sampleNumber']],
                                         min.exonDistance = ir.control[['min.exonDistance']],
                                         min.exonOverlap = ir.control[['min.exonOverlap']],
                                         prefix = ir.control[['prefix']],
                                         verbose = verbose)
  return(annotations)
}

#' Perform quantification
#' @inheritParams bambu
#' @noRd
bambu.quantify <- function(readClass, annotations, algo.control, min.exonDistance=35, ncore = 1, verbose = FALSE){
  if(is.character(readClass)){
    readClass <- readRDS(file=readClass)
    seqlevelsStyle(readClass) <- seqlevelsStyle(annotations)[1]
  }
  
  readClass <- isore.estimateDistanceToAnnotations(readClass, annotations, min.exonDistance = min.exonDistance, verbose = verbose)
  dt <- getEmptyClassFromSE(readClass, annotations)
  colNameRC <- colnames(readClass)
  colDataRC <- colData(readClass)
  rm(readClass)
  gc(verbose=FALSE)
  counts <- bambu.quantDT(dt,algo.control = algo.control, ncore = ncore, verbose = verbose)
  if(length(setdiff(counts$tx_name,names(annotations)))>0){
    stop("The provided annotation is incomplete")
  }
  counts <- counts[data.table(tx_name = names(annotations)),  on = 'tx_name']
  counts[is.na(estimates),`:=`(estimates = 0, CPM = 0) ]
  
  seOutput <- SummarizedExperiment::SummarizedExperiment(assays = SimpleList(counts = matrix(counts$estimates,ncol = 1, dimnames = list(NULL, colNameRC)),
                                                                             CPM = matrix(counts$CPM, ncol =  1, dimnames = list(NULL, colNameRC))),
                                                         colData = colDataRC)
  return(seOutput)
}

#' Preprocess bam files and save read class files
#' @inheritParams bambu
#' @noRd
bambu.constructReadClass <- function(bam.file, genomeSequence, annotations, outputReadClassDir=NULL, stranded=FALSE, ncore=1, verbose=FALSE){
  

  readGrgList <- prepareDataFromBam(bam.file[[1]], cores=ncore, verbose = verbose)
  seqlevelsStyle(readGrgList) <- seqlevelsStyle(annotations)[1]
  se <- isore.constructReadClasses(readGrgList = readGrgList,
                                   runName = names(bam.file)[1],
                                   annotationGrangesList = annotations,
                                   genomeSequence = genomeSequence,
                                   stranded = stranded,
                                   ncore = ncore,
                                   verbose = verbose)
  seqlevels(se) <- unique(c(seqlevels(se), seqlevels(annotations)))
  if(!is.null(outputReadClassDir)){
    readClassFile <- fs::path(outputReadClassDir,paste0(names(bam.file),'_readClassSe'), ext='rds')
    if(file.exists(readClassFile)){
      show(paste(readClassFile, 'exists, will be overwritten')) #warning is not printed, use show in addition
      warning(paste(readClassFile, 'exists, will be overwritten'))
    }
    saveRDS(se, file=readClassFile)
    se <- readClassFile
  }
  return(se)
}
  
  

#' Process data.table object
#' @param dt A data.table object
#' @inheritParams bambu
#' @noRd
bambu.quantDT <- function(dt = dt,algo.control = NULL,ncore = 1, verbose = FALSE){
  if(is.null(dt)){
    stop("Input object is missing.")
  }else if(any(!(c('gene_id','tx_id','read_class_id','nobs') %in% colnames(dt)))){
    stop("Columns gene_id, tx_id, read_class_id, nobs, are missing from object.")
  }

  ## check quantification parameters
  algo.control.default <- list(bias_correction = FALSE,
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
                                      ncore = ncore,
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



  #b_est <- outList[[2]]
  #b_est[, `:=`(gene_name = geneVec[gene_sid], eqClass = eqClassVec[as.numeric(read_class_sid)])]
  #b_est[, `:=`(gene_sid = NULL,read_class_sid=NULL)]

  #est.list <- list(counts = theta_est,
  #                 metadata = b_est)
  #return(est.list)
  return(theta_est)
}


