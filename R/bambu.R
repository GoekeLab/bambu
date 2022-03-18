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
#' \code{BamFileList}  object (see \code{Rsamtools}).
#' @param rcFile A string or a vector of strings specifying the read
#' class files that are saved during previous run of \code{\link{bambu}}.
#' @param rcOutDir A string variable specifying the path to where
#' read class files will be saved.
#' @param annotations A \code{TxDb} object or A GRangesList object
#' obtained by \code{\link{prepareAnnotations}}.
#' @param genome A fasta file or a BSGenome object.
#' @param stranded A boolean for strandedness, defaults to FALSE.
#' @param ncore specifying number of cores used when parallel processing 
#' is used, defaults to 1.
#' @param NDR specifying the maximum NDR rate to novel transcript
#'     output from detected read classes, defaults to 0.1
#' @param yieldSize see \code{Rsamtools}.
#' @param opt.discovery A list of controlling parameters for isoform
#' reconstruction process:
#' \describe{
#'     \item{prefix}{specifying prefix for new gene Ids (genePrefix.number),
#'     defaults to empty}
#'     \item{remove.subsetTx}{indicating whether filter to remove read classes
#'     which are a subset of known transcripts(), defaults to TRUE}
#'     \item{min.readCount}{specifying minimun read count to consider a read
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
#' }
#' @param trackReads When TRUE read names will be tracked and output as
#' metadata in the final output as readToTranscriptMaps detailing. 
#' the assignment of reads to transcripts. The output is a list with 
#' an entry for each sample.
#' @param outputDistTable When TRUE the calculated distance table between
#' read classes and annotations will be output as metadata as 
#' distTables. The output is a list with an entry for each sample.
#' @param discovery A logical variable indicating whether annotations
#' are to be extended
#' @param quant A logical variable indicating whether quantification will 
#' be performed
#' @param verbose A logical variable indicating whether processing messages will
#' be printed.
#' @param lowMemory Read classes will be processed by chromosomes when lowMemory 
#' is specified. This option provides an efficient way to process big samples.
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
#'     \item{partialLengthCounts}{estimates of read counts mapped as partial 
#'     length reads for each transcript}
#'     \item{uniqueCounts}{counts of reads that are uniquely mapped to each 
#'     transcript}
#'     \item{theta}{raw estimates}
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
bambu <- function(reads = NULL, rcFile = NULL, rcOutDir = NULL,
    annotations = NULL, genome = NULL, stranded = FALSE, ncore = 1, NDR = 0.1,
    yieldSize = NULL, opt.discovery = NULL, opt.em = NULL, trackReads = FALSE, 
    returnDistTable = FALSE, discovery = TRUE, quant = TRUE, verbose = FALSE, 
    lowMemory = FALSE) {
    if (!(discovery+quant)) stop("At least 1 of discovery and quant must be 
    TRUE. Rerun with either 1 or both parameters as TRUE")
    if(is.null(annotations)) { annotations = GRangesList()
    } else annotations <- checkInputs(annotations, reads, readClass.file = rcFile,
            readClass.outputDir = rcOutDir, genomeSequence = genome)
    if(!is.null(reads)) genomeSequence <- checkInputSequence(genome)
    isoreParameters <- setIsoreParameters(isoreParameters = opt.discovery)

    #below line is to be compatible with earlier version of running bambu
    if(!is.null(isoreParameters$max.txNDR)) NDR = isoreParameters$max.txNDR

    emParameters <- setEmParameters(emParameters = opt.em)
    bpParameters <- setBiocParallelParameters(reads, readClass.file = rcFile,
        ncore, verbose)
    if (bpParameters$workers > 1) ncore <- 1
    rm.readClassSe <- FALSE
    if (!is.null(reads)) {
        if (length(reads) > 10 & (is.null(rcOutDir))) {
            rcOutDir <- tempdir() #>=10 samples, save to temp folder
            message("There are more than 10 samples, read class files
                will be temporarily saved to ", rcOutDir,
                " for more efficient processing")
            rm.readClassSe <- TRUE # remove temporary read class files 
        }
        readClassList <- bambu.processReads(reads, annotations, 
            genomeSequence = genomeSequence, 
            readClass.outputDir = rcOutDir, yieldSize, 
            bpParameters, stranded, verbose,
            isoreParameters, trackReads = trackReads)
    } else { 
        if(is.list(rcFile)) {readClassList <- rcFile}
        else {readClassList <- as.list(rcFile)}}
    if (discovery) {
        annotations <- bambu.extendAnnotations(readClassList, annotations, NDR,
            isoreParameters, stranded, bpParameters, verbose = verbose)
        if (!verbose) message("Finished extending annotations.")
        if (!quant){
            return(annotations=annotations)
        }
    }
    if (quant) {
        if (!verbose) message("Start isoform quantification")
        if(length(annotations)==0) stop("No valid annotations, if running
                                de novo please try less stringent parameters")
        countsSe <- bplapply(readClassList, bambu.quantify,
            annotations = annotations, isoreParameters = isoreParameters,
            emParameters = emParameters, ncore = ncore, verbose = verbose, 
            BPPARAM = bpParameters)
        if(trackReads){ 
            readToTranscriptMaps = bplapply(Map(list,readClassList,countsSe), generateReadToTranscriptMap,
                annotations, BPPARAM = bpParameters)
        }
        if(returnDistTable){
            distTables = lapply(countsSe, FUN = function(se){metadata(se)$distTable})}
        countsSe = lapply(countsSe, FUN = function(se){
            metadata(se)$distTable=NULL
            return(se)})
        countsSe <- do.call(SummarizedExperiment::cbind, countsSe)
        if(returnDistTable) metadata(countsSe)$distTables = distTables
        rowRanges(countsSe) <- annotations
        if(trackReads) metadata(countsSe)$readToTranscriptMaps = readToTranscriptMaps
        if (!verbose) message("Finished isoform quantification.")
        if (rm.readClassSe) file.remove(unlist(readClassList))
        return(countsSe)
    }
}
