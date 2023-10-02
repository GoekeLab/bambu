#' Main function
#' @title long read isoform reconstruction and quantification
#' @description This function takes bam file of genomic alignments and performs
#' isoform recontruction and gene and transcript expression quantification.
#' It also allows saving of read class files of alignments, extending provided
#' annotations, and quantification based on extended annotations. When multiple
#' samples are provided, extended annotations will be combined across samples to
#' allow comparison.
#' @param reads A string or a vector of strings specifying the paths of bam
#' files for genomic alignments, or a \code{BamFile} object or a
#' \code{BamFileList}  object (see \code{Rsamtools}). Alternatively 
#' string or a vector of strings specifying the read class files 
#' that are saved during previous run of \code{\link{bambu}}.
#' @param annotations A path to a .gtf file or a \code{TxDb} object 
#' or a GRangesList object obtained by \code{\link{prepareAnnotations}}.
#' @param genome A path to a fasta file or a BSGenome object.
#' @param NDR specifying the maximum NDR rate to novel transcript
#' output from detected read classes, defaults to an automatic recommendation
#' @param opt.discovery A list of controlling parameters for isoform
#' reconstruction process:
#' \describe{
#'     \item{remove.subsetTx}{indicating whether filter to remove read classes
#'     which are a subset of known transcripts(), defaults to TRUE}
#'     \item{min.readCount}{specifying minimum read count to consider a read
#'     class valid in a sample, defaults to 2}
#'     \item{min.readFractionByGene}{specifying minimum relative read count per
#'     gene, highly expressed genes will have many high read count low relative
#'     abundance transcripts that can be filtered, defaults to 0.05}
#'     \item{min.sampleNumber}{specifying minimum sample number with minimum read
#'     count, defaults to 1}
#'     \item{min.exonDistance}{specifying minum distance to known transcript 
#'     to be considered valid as new, defaults to 35bp}
#'     \item{min.exonOverlap}{specifying minimum number of bases shared with
#'     annotation to be assigned to the same gene id, defaults to 10bp}
#'     \item{min.primarySecondaryDist}{specifying the minimum number of distance 
#'     threshold, defaults to 5bp}
#'     \item{min.primarySecondaryDistStartEnd1}{specifying the minimum number 
#'     of distance threshold, used for extending annotation, defaults to 5bp}
#'     \item{min.primarySecondaryDistStartEnd2}{specifying the minimum number 
#'     of distance threshold, used for estimating distance to annotation, 
#'     defaults to 5bp}
#'     \item{min.txScore.multiExon}{specifying the minimum transcript level 
#'     threshold for multi-exon transcripts during sample combining, 
#'     defaults to 0}
#'     \item{min.txScore.singleExon}{specifying the minimum transcript level 
#'     threshold for single-exon transcripts during sample combining, defaults 
#'     to 1}
#'     \item{fitReadClassModel}{ A boolean specifying if Bambu should attempt
#'     to train a transcript discovery model for all samples. Defaults to TRUE}
#'     \item{defaultModels}{A model object obtained by code{\link{trainBambu}}
#'     or when returnModel is TRUE}
#'     \item{returnModel}{A boolean specifying if the trained model is output
#'     with the readclass files. Defaults to FALSE}
#'     \item{baselineFDR}{A number between 0 - 1, specifying the false discovery
#'     rate used during NDR recomendation. Defaults to 0.1}
#'     \item{min.readFractionByEqClass}{indicating the minimum relative read
#'     count of a subset transcript compared to all superset transcripts 
#'     (ie the relative read count within the minimum equivalent class). This 
#'     filter is applied on the set of annotations across all samples using the 
#'     total read count, this is not a per-sample filter. Please use with 
#'     caution. defaults to 0}
#'     \item{prefix}{specifying prefix for new gene Ids (genePrefix.number),
#'     defaults to "Bambu"}
#' }
#' @param opt.em A list of controlling parameters for quantification
#' algorithm estimation process:
#' \describe{
#'     \item{maxiter}{specifying maximum number of run iterations,
#'     defaults to 10000}
#'     \item{degradationBias}{correcting for degradation bias, defaults to TRUE}
#'     \item{conv}{specifying the covergence threshold control, 
#'     defaults to 0.0001}
#'     \item{minvalue}{specifying the minvalue for convergence consideration, 
#'     defaults to 0.00000001}
#'     \item{sig.digit}{specifying the maximum significant digits of the reported estimates}
#' }
#' @param rcOutDir A string variable specifying the path to where
#' read class files will be saved.
#' @param discovery A logical variable indicating whether annotations
#' are to be extended. Defaults to TRUE
#' @param quant A logical variable indicating whether quantification will 
#' be performed. If false the output type will change. Defaults to TRUE
#' @param stranded A boolean for strandedness, defaults to FALSE.
#' @param ncore specifying number of cores used when parallel processing 
#' is used, defaults to 1.
#' @param yieldSize see \code{Rsamtools}.
#' @param trackReads When TRUE read names will be tracked and output as
#' metadata in the final output as readToTranscriptMaps detailing. 
#' the assignment of reads to transcripts. The output is a list with 
#' an entry for each sample.
#' @param returnDistTable When TRUE the calculated distance table between
#' read classes and annotations will be output as metadata as 
#' distTables. The output is a list with an entry for each sample.
#' @param lowMemory Read classes will be processed by chromosomes when lowMemory 
#' is specified. This option provides an efficient way to process big samples.
#' @param fusionMode A logical variable indicating whether run in fusion mode
#' @param verbose A logical variable indicating whether processing messages will
#' be printed.
#' @details
#' @return \code{bambu} will output different results depending on whether
#' \emph{quant} mode is on. By default, \emph{quant} is set to TRUE, so 
#' \code{bambu} will generate a \emph{SummarizedExperiment} object that contains
#' the transcript expression estimates. Transcript expression estimates can be 
#' accessed by \emph{counts()}, including the following variables
#' \describe{
#'     \item{counts}{expression estimates}
#'     \item{CPM}{sequencing depth normalized estimates}
#'     \item{fullLengthCounts}{estimates of read counts mapped as full length 
#'     reads for each transcript}
#'     \item{uniqueCounts}{counts of reads that are uniquely mapped to each 
#'     transcript}
#' }
#' Output annotations that are usually the annotations with/without novel 
#' transcripts/genes added, depending on whether \emph{discovery} mode is on
#' can be accessed by \emph{rowRanges()}
#' Transcript to gene map can be accessed by \emph{rowData()}, with 
#' \emph{eqClass} that defining equivalent class for each transcript
#' 
#' In the case when \emph{quant} is set to FALSE, i.e., only transcript 
#' discovery is performed, \code{bambu} will report the \emph{grangeslist} of 
#' the extended annotations.
#' @importFrom BiocParallel bplapply
#' @importFrom SummarizedExperiment cbind
#' @examples
#' ## =====================
#' test.bam <- system.file("extdata",
#'     "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam",
#'     package = "bambu")
#' fa.file <- system.file("extdata", 
#'     "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", 
#'     package = "bambu")
#' gr <- readRDS(system.file("extdata", 
#'     "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds",
#'     package = "bambu"))
#' se <- bambu(reads = test.bam, annotations = gr, 
#'     genome = fa.file,  discovery = TRUE, quant = TRUE)
#' @export
bambu <- function(reads, annotations = NULL, genome = NULL, NDR = NULL,
                  opt.discovery = NULL, opt.em = NULL, rcOutDir = NULL, discovery = TRUE, 
                  quant = TRUE, stranded = FALSE,  ncore = 1, yieldSize = NULL,  
                  trackReads = FALSE, returnDistTable = FALSE, lowMemory = FALSE, 
                  fusionMode = FALSE, verbose = FALSE, demultiplexed = FALSE, readGrgListFile = NULL) {
  
  ### add ###
  library(Matrix)
  
  if(demultiplexed == TRUE) {trackReads = TRUE; lowMemory = TRUE}
  ### add ### 
  
  if(is.null(annotations)) { annotations = GRangesList()
  } else annotations <- checkInputs(annotations, reads,
                                    readClass.outputDir = rcOutDir, genomeSequence = genome)
  isoreParameters <- setIsoreParameters(isoreParameters = opt.discovery)
  #below line is to be compatible with earlier version of running bambu
  if(!is.null(isoreParameters$max.txNDR)) NDR = isoreParameters$max.txNDR
  
  emParameters <- setEmParameters(emParameters = opt.em)
  
  ### add ###
  if (is.null(readGrgListFile)){
    bpParameters <- setBiocParallelParameters(reads, ncore, verbose) 
    
  } else{
    bpParameters <- setBiocParallelParameters(readGrgListFile, ncore, verbose) 
  }
  ### add ### 
  
  rm.readClassSe <- FALSE
  readClassList = reads
  isRDSs = all(sapply(reads, class)=="RangedSummarizedExperiment")
  isBamFiles = !isRDSs
  if(!isRDSs) isBamFiles = ifelse(!is(reads, "BamFileList"), all(grepl(".bam$", reads)), FALSE)
  if (isBamFiles | is(reads, "BamFileList")) {
    if (length(reads) > 10 & (is.null(rcOutDir))) {
      rcOutDir <- tempdir() #>=10 samples, save to temp folder
      message("There are more than 10 samples, read class files
                will be temporarily saved to ", rcOutDir,
              " for more efficient processing")
      rm.readClassSe <- TRUE # remove temporary read class files 
    }
    message("--- Start generating read class files ---")
    
    readClassList <- bambu.processReads(reads, annotations,
       genomeSequence = genome,
       readClass.outputDir = rcOutDir, yieldSize = yieldSize,
       bpParameters, stranded, verbose,
       isoreParameters, trackReads = trackReads, fusionMode = fusionMode,
       lowMemory = lowMemory, demultiplexed = demultiplexed, readGrgListFile = readGrgListFile)
  }
  if (is.null(readGrgListFile)) warnings = handleWarnings(readClassList, verbose)
  if (!discovery & !quant) return(readClassList)
  if (discovery) {
    message("--- Start extending annotations ---")
    annotations <- bambu.extendAnnotations(readClassList, annotations, NDR,
                                           isoreParameters, stranded, bpParameters, fusionMode, verbose)
    metadata(annotations)$warnings = warnings
    
    if (!quant) return(annotations)
  }
  
  if (quant) {
    
    message("--- Start isoform quantification ---")
    if(length(annotations)==0) stop("No valid annotations, if running
                                de novo please try less stringent parameters")
    
    ### add ### 
    if (!is.null(readGrgListFile)){
  
      library(BiocFileCache) # temporary
      bfc <- BiocFileCache(rcOutDir, ask = FALSE)
      info <- bfcinfo(bfc)
      
      countsSeCompressed <- bplapply(info$rpath, bambu.quantify,
                                     annotations = annotations, isoreParameters = isoreParameters,
                                     emParameters = emParameters, trackReads = trackReads, 
                                     returnDistTable = returnDistTable, verbose = verbose, 
                                     readGrgListFile = readGrgListFile, 
                                     BPPARAM = bpParameters)
      
      ### countsData
      countsData <- c("incompatibleCounts", "counts", "CPM", "fullLengthCounts", "uniqueCounts")
      countsDataMat <- list()
      
      for (k in seq_along(countsData)){
        countsVecList <- lapply(seq_along(countsSeCompressed), function(j){countsSeCompressed[[j]][[countsData[k]]]})
        
        countsMat <- sparseMatrix(i = unlist(lapply(seq_along(countsVecList), function(j){countsVecList[[j]]@i})),
                                  j = unlist(lapply(seq_along(countsVecList), function(j){rep(j, length(countsVecList[[j]]@i))})),
                                  x = unlist(lapply(seq_along(countsVecList), function(j){countsVecList[[j]]@x})),
                                  dims = c(length(countsVecList[[1]]), length(countsVecList)))
        
        colnames(countsMat) <- unlist(lapply(seq_along(countsSeCompressed), function(j){countsSeCompressed[[j]]$colnames}))
        
        if (countsData[k] == "incompatibleCounts"){
          countsMat <- data.table(as.data.frame(as.matrix(countsMat)) %>%
                                    mutate(GENEID = unique(mcols(annotations)$GENEID)) %>%
                                    select(GENEID, everything()))
        }
        
        countsDataMat[[countsData[k]]] <- countsMat 
      }
      
      countsSe <- SummarizedExperiment(assays = SimpleList(counts = countsDataMat$counts, 
                                                           CPM = countsDataMat$CPM, 
                                                           fullLengthCounts = countsDataMat$fullLengthCounts, 
                                                           uniqueCounts = countsDataMat$uniqueCounts), 
                                       colData = as(bind_rows(lapply(lapply(seq_along(countsSeCompressed), 
                                                 function(j){countsSeCompressed[[j]]$colData}), as.data.frame)), "DataFrame"))
      
      metadata(countsSe)$incompatibleCounts <- countsDataMat$incompatibleCounts
      rowRanges(countsSe) <- annotations
      metadata(countsSe)$warnings = warnings
        
    } else {
    ### add ###
      countsSe <- bplapply(readClassList, bambu.quantify,
                           annotations = annotations, isoreParameters = isoreParameters,
                           emParameters = emParameters, trackReads = trackReads, 
                           returnDistTable = returnDistTable, verbose = verbose, 
                           readGrgListFile = readGrgListFile, 
                           BPPARAM = bpParameters)
      
      countsSe <- combineCountSes(countsSe, trackReads, returnDistTable)
      rowRanges(countsSe) <- annotations
      metadata(countsSe)$warnings = warnings
    }
    
    ### add ###
    #if (!is.null(readGrgListFile)) unlink(unique(dirname(readGrgListFile)), recursive = TRUE)
    ### add ###
    
    #if (rm.readClassSe) file.remove(unlist(readClassList))
    message("--- Finished running Bambu ---")
  
    return(countsSe)
  }
}  
