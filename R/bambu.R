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
#' @param yieldSize see \code{Rsamtools}.
#' @param opt.discovery A list of controlling parameters for isoform
#' reconstruction process:
#' \itemize{
#'     \item prefix specifying prefix for new gene Ids (genePrefix.number),
#'     defaults to empty
#'     \item remove.subsetTx indicating whether filter to remove read classes
#'     which are a subset of known transcripts(), defaults to TRUE
#'     \item min.readCount specifying minimun read count to consider a read
#'     class valid in a sample, defaults to 2
#'     \item min.readFractionByGene specifying minimum relative read count per
#'     gene, highly expressed genes will have many high read count low relative
#'     abundance transcripts that can be filtered, defaults to 0.05
#'     \item min.sampleNumber specifying minimum sample number with minimum read
#'     count, defaults to 1
#'     \item min.exonDistance specifying minum distance to known transcript 
#'     to be considered valid as new, defaults to 35
#'     \item min.exonOverlap specifying minimum number of bases shared with
#'     annotation to be assigned to the same gene id, defaults 10 base pairs
#'     \item min.primarySecondaryDist specifying the minimum number of distance 
#'     threshold
#'     \item min.primarySecondaryDistStartEnd1 specifying the minimum number 
#'     of distance threshold, used for extending annotation
#'     \item min.primarySecondaryDistStartEnd2 specifying the minimum number 
#'     of distance threshold, used for estimating distance to annotation
#' }
#' @param opt.em A list of controlling parameters for quantification
#' algorithm estimation process:
#' \itemize{
#'     \item maxiter specifying maximum number of run interations,
#'     defaults to 10000.
#'     \item degradationBias correcting for degradation bias, defaults to TRUE.
#'     \item conv specifying the covergence trheshold control,
#'     defaults to 0.0001.
#'     \item minvalue specifying the minvalue for convergence consideration
#' }
#' @param discovery A logical variable indicating whether annotations
#' are to be extended for quantification.
#' @param verbose A logical variable indicating whether processing messages will
#' be printed.
#' @details
#' @return A list of two SummarizedExperiment object for transcript expression
#' and gene expression.
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
#'     genome = fa.file,  discovery = FALSE)
#' @export
bambu <- function(reads = NULL, rcFile = NULL, rcOutDir = NULL,
    annotations = NULL, genome = NULL, stranded = FALSE, ncore = 1,
    yieldSize = NULL, opt.discovery = NULL, opt.em = NULL,
    discovery = TRUE, verbose = FALSE) {
    annotations <- checkInputs(annotations, reads, readClass.file = rcFile, 
            readClass.outputDir = rcOutDir, genomeSequence = genome)
    isoreParameters <- setIsoreParameters(isoreParameters = opt.discovery)
    emParameters <- setEmParameters(emParameters = opt.em)
    bpParameters <- setBiocParallelParameters(reads, readClass.file = rcFile,
        ncore, verbose)
    if (bpParameters$workers > 1) ncore <- 1
    rm.readClassSe <- FALSE
    if (!is.null(reads)) {
        if (length(reads) > 10 & (is.null(rcOutDir))) {
            rcOutDir <- tempdir() #>=10 samples, save to temp folder
            message(paste0("There are more than 10 samples, read class files
                will be temporarily saved to ", rcOutDir,
                " for more efficient processing"))
            rm.readClassSe <- TRUE # remove temporary read class files 
        }
        readClassList <- bambu.processReads(reads, annotations, 
            genomeSequence = genome, 
            readClass.outputDir = rcOutDir, yieldSize, 
            bpParameters, stranded, verbose)
    } else {
        readClassList <- rcFile
    }
    if (discovery) {
        annotations <- bambu.extendAnnotations(readClassList, annotations,
            isoreParameters, verbose = verbose)
        if (!verbose) message("Finished extending annotations.")
    }
    if (!verbose) message("Start isoform quantification")
    countsSe <- bplapply(readClassList,
        bambu.quantify,annotations = annotations,
        min.exonDistance = isoreParameters[["min.exonDistance"]],
        min.primarySecondaryDist =
            isoreParameters[['min.primarySecondaryDist']], 
        min.primarySecondaryDistStartEnd =
            isoreParameters[['min.primarySecondaryDistStartEnd2']],
        emParameters = emParameters, ncore = ncore,
        verbose = verbose, BPPARAM = bpParameters)
    countsSe <- do.call(SummarizedExperiment::cbind, countsSe)
    rowRanges(countsSe) <- annotations
    if (!verbose) message("Finished isoform quantification.")
    if (rm.readClassSe) file.remove(unlist(readClassList))#Clean temp directory
    return(countsSe)
}
