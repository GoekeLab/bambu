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
#' @param opt.rcAssignment A list of controlling parameters for the read class
#' to transcript assignment process:
#' \describe{
#'     \item{min.exonDistance}{specifying minimum distance to known transcript
#'     to be considered a valid match, defaults to 35bp}
#'     \item{min.primarySecondaryDist}{specifying the minimum distance
#'     threshold between primary and secondary assignments, defaults to 5bp}
#'     \item{min.primarySecondaryDistStartEnd2}{specifying the minimum
#'     distance threshold for start/end positions used for read assignment,
#'     defaults to 5bp}
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
#' @param sampleData A character vector of paths to metadata CSV files (or \code{NA} if 
#' unavailable for specific samples); defaults to \code{NULL}. Files must contain a 
#' "sampleName" column for bulk data or a "barcode" column for single-cell/spatial data. 
#' For bulk data, one metadata CSV file for all samples is sufficient, whereas single-cell/spatial 
#' data requires one metadata CSV file per sample.
#' @param fusionMode A logical variable indicating whether run in fusion mode
#' @param verbose A logical variable indicating whether processing messages will
#' be printed.
#' @param opt.singlecell A list of single-cell specific parameters:
#' \describe{
#'     \item{extractBarcodeUMI}{Logical, whether to extract cell barcodes and
#'     UMIs from BAM tags or read names. Defaults to FALSE}
#'     \item{dedupUMI}{Logical, whether to perform UMI-based deduplication per
#'     barcode. Defaults to FALSE}
#'     \item{clusters}{A named list mapping cluster names to barcode vectors,
#'     used for cluster-level transcript discovery. Defaults to NULL}
#' }
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
    mode = NULL, opt.discovery = NULL, opt.rcAssignment = NULL, opt.em = NULL, opt.singlecell = NULL,
    rcOutDir = NULL, discovery = TRUE, assignDist = TRUE, quant = TRUE, stranded = FALSE,  
    ncore = 1, yieldSize = NULL, trackReads = FALSE, returnDistTable = FALSE, lowMemory = FALSE, 
    sampleData = NULL, fusionMode = FALSE, verbose = FALSE, quantData = NULL,
    processByChromosome = FALSE, processByBam = TRUE) {
    message(paste0("Running Bambu-v", "3.9.0"))
    if(!is.null(mode)){
        if(mode == "bulk"){
            processByChromosome <- FALSE
            processByBam <- TRUE
        }
        if(mode == "multiplexed"){
            opt.singlecell$extractBarcodeUMI <- TRUE
            opt.em <- list(degradationBias = FALSE)
            quant <- FALSE
            processByChromosome <- TRUE
        }
        if(mode == "fusion"){
            NDR <- 1
            fusionMode <- TRUE
            if(is.null(opt.discovery)) opt.discovery <- list()
            opt.discovery$remove.subsetTx <- FALSE
            opt.discovery$min.readCount <- 1
            opt.discovery$min.sampleNumber <- 0
        }
        if(mode == "debug"){
            verbose <- TRUE
            trackReads <- TRUE
            returnDistTable <- TRUE
        }
    }
    if(lowMemory)
        message("lowMemory has been deprecated and split into processByChromosome and processByBam. Please see Documentation")
    if("min.primarySecondaryDistStartEnd2" %in% names(opt.discovery))
        message("min.primarySecondaryDistStartEnd2 has been moved to opt.rcAssignment. Please pass this parameter via opt.rcAssignment instead.")
    if(is.null(annotations)){ 
        annotations <- GRangesList()
    } else {
        annotations <- checkInputs(annotations, reads,
            readClass.outputDir = rcOutDir,
            genomeSequence = genome, discovery = discovery,
            sampleData = sampleData, quantData = quantData)
    }
    opt.discovery <- setDiscoveryParameters(discoveryParameters = opt.discovery)
    #below line is to be compatible with earlier version of running bambu
    if(!is.null(isoreParameters$max.txNDR)) NDR = isoreParameters$max.txNDR
    extractBarcodeUMI <- isTRUE(opt.singlecell$extractBarcodeUMI)
    dedupUMI <- isTRUE(opt.singlecell$dedupUMI)
    clusters <- opt.singlecell$clusters
    opt.rcAssignment <- setRcAssignmentParameters(rcAssignmentParameters = opt.rcAssignment)
    opt.em <- setEmParameters(emParameters = opt.em)
    bpParameters <- setBiocParallelParameters(reads, ncore, verbose, extractBarcodeUMI)
	xgb.set.config(nthread = 1)
    # only when reads is not NULL, this proceed, otherwise, it will jump to quant step
    if(!is.null(reads)){ 
        rm.readClassSe <- FALSE
        readClassList <- reads
        isRDSs <- all(sapply(reads, class)=="RangedSummarizedExperiment")
        isBamFiles <- !isRDSs
        warnings <- NULL
        if(!isRDSs) 
            isBamFiles <- ifelse(!is(reads, "BamFileList"), 
                                 all(grepl(".bam$", reads)), FALSE)
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
                                                bpParameters = bpParameters, stranded = stranded, verbose = verbose,
                                                discoveryParameters = opt.discovery, trackReads = trackReads,
                                                fusionMode = fusionMode, 
                                                processByChromosome = processByChromosome, processByBam = processByBam, 
                                                extractBarcodeUMI = extractBarcodeUMI,
                                                dedupUMI = dedupUMI)
        }
        
        #warnings = handleWarnings(readClassList, verbose)
        if (!discovery & !assignDist & !quant) return(readClassList)
        if (discovery) {
            message("--- Start extending annotations ---")
            extendedAnnotations <- bambu.extendAnnotations(readClassList, annotations, NDR,
                                                           opt.discovery, stranded, bpParameters, fusionMode, verbose)
            metadata(extendedAnnotations)$warnings = warnings
            
            #### cluster based transcript discovery
            if(!is.null(clusters)){
                annotations.clusters <- isore.extendAnnotations.clusters(readClassList,
                                                                         annotations, clusters, NDR,
                                                                         opt.discovery, stranded, bpParameters, fusionMode, verbose = FALSE)
                metadata(extendedAnnotations)$clusters <- annotations.clusters    
            }
            annotations <- extendedAnnotations
            
            if (!quant & !assignDist) return(annotations)
        }
        if(assignDist){
            message("--- Start calculating equivilance classes ---")
            quantData <- bplapply(seq_along(readClassList), function(i){
              assignReadClasstoTranscripts(
                readClassList = readClassList[[i]],
                annotations = annotations,
                rcAssignmentParameters = opt.rcAssignment,
                verbose = verbose,
                # for bulk data, there is one sampleData (keep sampleData[1]), for single-cell, there is one per sample
                sampleMetadata = if(length(sampleData) == 1) sampleData[1] else sampleData[i],
                extractBarcodeUMI = extractBarcodeUMI,
                returnDistTable = returnDistTable,
                trackReads = trackReads
              )
            }, BPPARAM = bpParameters)
            if (!quant) return(quantData)
        }
    }
    
    if (quant) {
        message("--- Start isoform EM quantification ---")
        if(!is.null(NDR) & !discovery)# this step is used when reset NDR is needed 
            annotations <- setNDR(annotations, NDR,
                                  prefix = opt.discovery$prefix,
                baselineFDR = opt.discovery[["baselineFDR"]],
                defaultModels2 = opt.discovery[["defaultModels"]])
        if(length(annotations)==0) stop("No valid annotations, if running
                                    de novo please try less stringent parameters")
        if(is.null(quantData)) stop("quantData must be provided or assignDist = TRUE")
        GENEIDs.i <- as.numeric(factor(unique(mcols(annotations)$GENEID)))
        start.ptm <- proc.time()
        countsSeCompressed.all <- NULL
        ColNames <- c()
        colData.all <- list()
        for(i in seq_along(quantData)){
            quantData_i <- quantData[[i]]
            #load in the barcode clustering from file if provided
            iter <- seq_len(ncol(metadata(quantData_i)$countMatrix)) # iter is integer
            if(!is.null(clusters)){
              if(class(clusters[[i]])!="CompressedCharacterList"){ # !is.list(clusters) is FALSE for CompressedCharacterList
                clusterMaps <- NULL
                for(j in seq_along(metadata(quantData_i)$sampleNames)){ #load in a file per sample name provided
                  clusterMap <- fread(clusters[[j]], header = FALSE,
                                      data.table = FALSE)
                  # read.table(clusters[[j]],
                  #     sep = ifelse(grepl(".tsv$",clusters[[j]]), "\t", ","),
                  #     header = FALSE)
                  clusterMap[,1] <- paste0(metadata(quantData_i)$sampleNames[j],
                                           "_",clusterMap[,1])
                  clusterMaps <- rbind(clusterMaps, clusterMap)
                }
                clustering <- splitAsList(clusterMaps[,1], clusterMaps[,2])
                rm(clusterMaps)
                rm(clusterMap)
                iter <- clustering

              } else{ #if clusters is a list
                iter <- clusters[[i]]
              }
            }
            countsSeCompressed <- bplapply(iter, FUN = function(j){ # previous i changed to j to avoid duplicated assignment 
                #i = iter[i %in% colnames(metadata(quantData_i)$countMatrix)] #bug, after assignment, i become emptyprint(i)
                countMatrix <- unname(metadata(quantData_i)$countMatrix[,j]) # same here 
                incompatibleCountMatrix <- unname(metadata(quantData_i)$incompatibleCountMatrix[,j]) # same here
                if(!is.null(dim(countMatrix))){
                    countMatrix <- rowSums(countMatrix)
                    incompatibleCountMatrix <- rowSums(metadata(quantData_i)$incompatibleCountMatrix[,j]) # same here
                }
                return(bambu.quantify(readClassDt = metadata(quantData_i)$readClassDt, countMatrix = countMatrix, 
                                            incompatibleCountMatrix = data.table(GENEID.i = as.numeric(rownames(metadata(quantData_i)$incompatibleCountMatrix)), counts = incompatibleCountMatrix),
                                            txid.index = mcols(annotations)$txid, GENEIDs = GENEIDs.i,
                                            emParameters = opt.em, trackReads = trackReads, 
                                            verbose = verbose))}, 
                                            BPPARAM = bpParameters)
            end.ptm <- proc.time()
            message("Total Time ", round((end.ptm - start.ptm)[3] / 60, 3), " mins.")
            if(!is.null(clusters)){
                ColNames <- c(ColNames, names(iter))
                colData.all[[i]] <- data.frame(
                  id = names(countsSeCompressed), 
                  sampleName = names(countsSeCompressed),
                  row.names = names(countsSeCompressed)
                )
            } else{
                ColNames <- c(ColNames, colnames(quantData_i)) 
                colData.all[[i]] <- data.frame(colData(quantData_i))
            }
            countsSeCompressed.all <- c(countsSeCompressed.all, countsSeCompressed)
        }
        names(countsSeCompressed.all) <- ColNames   
        
        countsSe <- combineCountSes(countsSeCompressed.all, colData.all, annotations)
        if(returnDistTable){
            distTables = list()
            for(i in seq_along(quantData)){
                distTables[[i]] <- metadata(quantData[[i]])$distTable
            }
            metadata(countsSe)$distTables <- distTables
        }
        return(countsSe)
    }
  }

#' Single-cell isoform reconstruction and quantification
#' @title Single-cell isoform reconstruction and quantification with Bambu
#' @description Analyse single-cell long-read RNA-seq data with Bambu,
#' performing isoform discovery and quantification at single-cell resolution.
#' This function calls the main \code{\link{bambu}} function with cell
#' barcode/UMI extraction and UMI-based deduplication enabled by default,
#' and returns a \emph{SummarizedExperiment} object with per-cell transcript
#' expression estimates.
#'
#' We recommend processing single-cell data using the Bambu Nextflow pipeline,
#' which handles preprocessing, barcode demultiplexing, and alignment prior to
#' running this function. See \url{https://github.com/GoekeLab/bambu-singlecell-spatial}.
#' @param reads A string or vector of strings specifying paths to BAM files.
#' BAM files must contain cell barcode and UMI information, either as BAM tags
#' (\code{CB} for cell barcode, \code{UB} for UMI) or encoded in the read name
#' using the format \code{CB_UMI#READNAME} (note: \code{CB} and \code{UMI}
#' must not contain underscores).
#' @param annotations A path to a .gtf file or a \code{TxDb} object for
#' transcript annotations. Defaults to NULL.
#' @param genome A path to a fasta file or a \code{BSGenome} object.
#' Defaults to NULL.
#' @param NDR Numeric specifying the maximum NDR rate for novel transcript
#' discovery. Defaults to NULL.
#' @param discovery Logical, whether transcript discovery is performed.
#' Defaults to TRUE.
#' @param assignDist Logical, whether to assign reads to transcripts.
#' Defaults to TRUE.
#' @param quant Logical, whether quantification is performed. Defaults to TRUE.
#' @param clusters A named list mapping cluster names to barcode vectors,
#' used for cluster-level transcript discovery. Defaults to NULL.
#' @param stranded Logical, whether reads are stranded. Defaults to FALSE.
#' @param ncore Integer specifying the number of cores for parallel processing.
#' Defaults to 1.
#' @param ... Additional arguments passed to \code{\link{bambu}}, such as
#' \code{verbose}, \code{lowMemory}, \code{opt.discovery}, etc.
#' @examples
#' sc.bam <- system.file("extdata", "demultiplexed.bam", package = "bambu")
#' fa.file <- system.file("extdata",
#'     "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa",
#'     package = "bambu")
#' gr <- readRDS(system.file("extdata",
#'     "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds",
#'     package = "bambu"))
#' se <- bambu.singlecell(reads = sc.bam, annotations = gr, genome = fa.file)
#' @export
bambu.singlecell <- function(reads, annotations = NULL, genome = NULL, NDR = NULL,
    clusters = NULL, discovery = TRUE, assignDist = TRUE, quant = TRUE,
    stranded = FALSE, ncore = 1, ...) {
    bambu(reads = reads, annotations = annotations, genome = genome, NDR = NDR,
        opt.singlecell = list(extractBarcodeUMI = TRUE, dedupUMI = TRUE, clusters = clusters),
        discovery = discovery, assignDist = assignDist, quant = quant,
        stranded = stranded, ncore = ncore, ...)
}
