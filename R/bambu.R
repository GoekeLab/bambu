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
#' @param assignDist A logical variable indicating whether read-class-to-transcript
#' assignment will be performed, defaults to TRUE
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
#' @param opt.singlecell A list of parameters for bambu's single-cell module
#' (\code{extractBarcodeUMI}, \code{dedupUMI}, and \code{clusters}). See
#' \code{\link{bambu.singlecell}} for details.
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
    processByChromosome = FALSE) {
    message(paste0("Running Bambu-v", "3.9.0"))
    if(!is.null(mode)){
        if(mode == "bulk"){
            processByChromosome <- FALSE
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
        message("lowMemory has been replaced by processByChromosome. Please see Documentation")
    if("min.primarySecondaryDistStartEnd2" %in% names(opt.discovery))
        message("min.primarySecondaryDistStartEnd2 has been moved to opt.rcAssignment. Please pass this parameter via opt.rcAssignment instead.")
    if(is.null(annotations)){ 
        annotations <- GRangesList()
    } else {
        annotations <- checkInputs(annotations, reads,
            readClass.outputDir = rcOutDir,
            genomeSequence = genome, discovery = discovery,
            sampleData = sampleData, quantData = quantData,
            clusters = opt.singlecell$clusters)
    }
    opt.discovery <- setDiscoveryParameters(discoveryParameters = opt.discovery)
    #below line is to be compatible with earlier version of running bambu
    if(!is.null(opt.discovery$max.txNDR)) NDR = opt.discovery$max.txNDR
    opt.rcAssignment <- setRcAssignmentParameters(rcAssignmentParameters = opt.rcAssignment)
    opt.em <- setEmParameters(emParameters = opt.em)
    extractBarcodeUMI <- isTRUE(opt.singlecell$extractBarcodeUMI)
    dedupUMI <- isTRUE(opt.singlecell$dedupUMI)
    clusters <- opt.singlecell$clusters
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
                                                processByChromosome = processByChromosome,
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
            names(quantData) <- names(readClassList)
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
        start.ptm <- proc.time()
        countsSeCompressed.all <- NULL
        ColNames <- c()
        colData.all <- list()
        if (!is.null(clusters)) {
            clusterDf <- if (is.atomic(clusters) && !is.null(names(clusters))) {
                data.frame(id = names(clusters), cluster = as.character(clusters))
            } else if (is.data.frame(clusters)) {
                clusters
            } else {
                fread(clusters)
            }
        }
        for(i in seq_along(quantData)){
            quantData_i <- quantData[[i]]
            if(!is.null(clusters)){
              # cluster mode: one column per cluster, pooling that cluster's cells/spots
              sampleName <- names(quantData)[i]
              # find each clustered cell/spot's column in this sample; those from other
              # samples don't match here and drop out as NA
              idx <- match(clusterDf$id, getSampleData(quantData_i)$id)
              inSample <- !is.na(idx)                 # cells/spots belonging to sample i
              # columnGroups is a named list, one entry per cluster, each holding
              # that cluster's columns, e.g.
              #   cluster0: 1 4 5 6 ...
              #   cluster1: 13 19 20 ...
              columnGroups <- split(idx[inSample], clusterDf$cluster[inSample])
              # each cluster's barcodes, in the same index order as in columnGroups
              clusterBarcodes <- split(getSampleData(quantData_i)$barcode[idx[inSample]],
                                       clusterDf$cluster[inSample])
              clusterLabels <- names(columnGroups)
              # two samples can both have a cluster, so prefix the sample name back
              # so the columns don't clash when samples are combined into a single SE
              # later
              names(columnGroups) <- paste0(sampleName, "_", clusterLabels)
            } else {
              # no clusters: columnGroups is an integer vector with one column per
              # input column, i.e. one per cell/spot for single-cell/spatial or one
              # per sample for bulk (1, 2, 3, ...)
              columnGroups <- seq_len(nrow(getSampleData(quantData_i)))
            }
            countsSeCompressed <- bplapply(columnGroups, FUN = function(columnIdx){ # previous i changed to j to avoid duplicated assignment

                incompatibleCounts_i <- getIncompatibleCounts(quantData_i)[, columnIdx, drop = FALSE]

                return(bambu.quantify(readClassDt = getReadClassDt(quantData_i), columnIdx = columnIdx,
                                            incompatibleCounts = data.table(GENEID.i = rownames(incompatibleCounts_i), counts = rowSums(incompatibleCounts_i)),
                                            txid.index = mcols(annotations)$txid, GENEIDs = rownames(incompatibleCounts_i),
                                            emParameters = opt.em, trackReads = trackReads,
                                            verbose = verbose))}, 
                                            BPPARAM = bpParameters)
            end.ptm <- proc.time()
            message("Total Time ", round((end.ptm - start.ptm)[3] / 60, 3), " mins.")
            if(!is.null(clusters)){
                ColNames <- c(ColNames, names(columnGroups))
                colData.all[[i]] <- data.frame(
                  id = names(columnGroups),
                  sampleName = sampleName,
                  cluster = clusterLabels,
                  row.names = names(columnGroups)
                )
                colData.all[[i]]$barcodes <- unname(clusterBarcodes)
            } else{
                ColNames <- c(ColNames, rownames(getSampleData(quantData_i)))
                colData.all[[i]] <- data.frame(getSampleData(quantData_i))
            }
            countsSeCompressed.all <- c(countsSeCompressed.all, countsSeCompressed)
        }
        names(countsSeCompressed.all) <- ColNames   
        
        countsSe <- combineCountSes(countsSeCompressed.all, colData.all, annotations)
        if(trackReads){
            readToTranscriptMaps = list()
            for(i in seq_along(quantData)){
                readToTranscriptMaps[[i]] <- getReadToTranscriptMap(quantData[[i]])
            }
            names(readToTranscriptMaps) <- names(quantData)
            metadata(countsSe)$readToTranscriptMaps <- readToTranscriptMaps
        }
        if(returnDistTable){
            distTables = list()
            for(i in seq_along(quantData)){
                distTables[[i]] <- getDistTable(quantData[[i]])
            }
            names(distTables) <- names(quantData)
            metadata(countsSe)$distTables <- distTables
        }
        return(countsSe)
    }
  }

#' @title Single-cell / Spatial transcript discovery and quantification with Bambu
#' @description A function wrapper for \code{\link{bambu}} that runs on long-read RNA-seq data to
#' perform transcript discovery and quantification at the single-cell/spatial level.
#' The function calls \code{\link{bambu}} with cell
#' barcode/UMI extraction and UMI-based deduplication enabled by default, with an optional
#' parameter for clustered-level transcript quantification.
#'
#' For general users, we recommend to use bambu-pipe, an end-to-end single-cell/spatial pipeline that
#' handles preprocessing, barcode demultiplexing, and alignment prior to
#' transcript discovery and quantification. See \url{https://github.com/GoekeLab/bambu-pipe}.
#' @inheritParams bambu
#' @param reads For \code{output} in \code{"readClasses"},
#' \code{"extendedAnnotations"}, or \code{"quantData"}, the input reads as in
#' \code{\link{bambu}}: paths to BAM files, a \code{BamFileList}, or the read-class
#' objects returned by the \code{"readClasses"} stage. For \code{output} in
#' \code{"uniqueCounts"}, \code{"EM"}, or \code{"clusteredEM"}, the per-sample
#' \code{quantData} list returned by the \code{"quantData"} stage.
#' @param output Bambu performs transcript discovery and quantification on single-cell
#' data in distinct steps, that are each specified using the \code{output} option. Must
#' be one of \code{"readClasses"}, \code{"extendedAnnotations"}, \code{"quantData"},
#' \code{"uniqueCounts"}, \code{"EM"}, or \code{"clusteredEM"}. See the \strong{Details}
#' section for what each stage does and the \strong{Value} section for what it returns.
#' @param clusters Assignment of cells or spatial spots to clusters, supplied as one of:
#' a named vector whose names are the cell/spot identifiers and whose values are the
#' cluster labels; a \code{data.frame} object with an \code{id} and a \code{cluster}
#' column; or a path to a \code{.csv}/\code{.tsv}/\code{.txt} file holding
#' that \code{data.frame}. In every case the identifier follows the
#' \code{sampleName_barcode} format and the cluster label is the cluster that cell/spot
#' belongs to. Must be used with \code{output = "clusteredEM"}, where expression is
#' aggregated to the cluster level instead of per individual cell/spot.
#' @param ... Additional arguments passed to \code{\link{bambu}}, such as
#' \code{sampleData}, \code{opt.em}, \code{trackReads}, \code{returnDistTable},
#' \code{yieldSize}, and \code{verbose}. See \code{\link{bambu}} for the full set.
#' @details
#' Single-cell and spatial long read RNA-Seq protocols tag every read with a barcode
#' and a unique molecular identifier (UMI). The barcode identifies which cell
#' (single-cell data) or spatial spot (spatial data) a read came from; below, "cell"
#' and "single-cell" refer to either. \code{\link{bambu.singlecell}} reads the barcode
#' and UMI for each read and collapses duplicate molecules, so that transcript
#' quantification can be performed per barcode.
#'
#' \code{\link{bambu.singlecell}} runs one stage at a time, each specified with the
#' \code{output} argument for the returned object. The step-by-step workflow is
#' described below and in the examples:
#' \enumerate{
#'     \item \strong{Read Class Construction} (\code{output = "readClasses"}): build
#'     per-sample read classes from the demultiplexed reads. Requires \code{reads}
#'     (paths to BAM files), \code{genome} (reference genome sequence), and
#'     \code{annotations} (optional, but recommended).
#'     \item \strong{Transcript Discovery} (\code{output = "extendedAnnotations"}):
#'     discover novel transcripts and extend the reference annotations. Requires
#'     \code{reads} (the read classes from step 1) and \code{annotations} (optional,
#'     but recommended; use the same annotations as in step 1, or omit them in both).
#'     \item \strong{Read to Transcript Assignment} (\code{output = "quantData"}):
#'     assign reads to the transcripts, giving the per-sample read-to-transcript
#'     assignments. Requires \code{reads} (the read classes from step 1) and
#'     \code{annotations} (the extended annotations from step 2 or the original
#'     annotations, which cannot be NULL if transcript discovery was skipped).
#'     \item \strong{Transcript Quantification} with the Expectation-Maximization
#'     (EM) algorithm, or without it (EM-free):
#'         \itemize{
#'            \item \strong{EM-free unique counts} (\code{output = "uniqueCounts"}):
#'            transcript-level unique counts, ready for downstream analysis. Apply
#'            \code{\link{transcriptToGeneExpression}} for gene counts, or cluster cells
#'            by gene expression to run the clustered EM. Requires \code{reads} (the
#'            quantData list from step 3) and the same \code{annotations} used in step 3.
#'            \item \strong{Single-Cell EM} (\code{output = "EM"}): estimates
#'            single-cell transcript expression with the EM. Requires \code{reads} (the
#'            quantData list from step 3) and the same \code{annotations} used in step 3.
#'            Single-cell EM is resource intensive; for higher efficiency use \code{uniqueCounts},
#'            or use \code{clusteredEM} to run the EM directly on cell clusters for more
#'            accurate, stable estimates.
#'            \item \strong{Clustered EM} (\code{output = "clusteredEM"}): runs the EM
#'            directly on user-defined cell clusters, giving more stable estimates than
#'            quantifying single cells and clustering afterwards. Requires \code{reads}
#'            (the quantData list from step 3), the same \code{annotations} used in
#'            step 3, and \code{clusters}.
#'         }
#' }
#' This makes it practical to run different stages in different R sessions to minimise
#' memory usage, and also lets users reuse intermediate results and explore different
#' downstream settings without reprocessing large datasets from the raw alignments. See
#' the \strong{Value} section for the returned output from each stage and \strong{Examples}
#' for a detailed step-by-step run.
#' @return The object returned depends on the specified \code{output} argument:
#' \describe{
#'     \item{\code{output = "readClasses"}}{returns a list of read-class objects for
#'     each sample that represent a summarised representation of similar reads with
#'     cell-barcode-to-read mapping information. ReadClass objects can be used as input
#'     for transcript discovery and read to transcript assignment.}
#'     \item{\code{output = "extendedAnnotations"}}{returns the extended annotations
#'     after transcript discovery as a \code{GRangesList} object. The extended
#'     annotation object can be used as input for transcript to read assignment and
#'     quantification.}
#'     \item{\code{output = "quantData"}}{returns a list of \code{quantData} objects
#'     that contain the read-to-transcript assignments for each sample. The quantData
#'     object is used to calculate unique counts per transcript and gene-level counts,
#'     and also can be used as input for EM quantification.}
#'     \item{\code{output = "uniqueCounts"}}{returns a \code{SummarizedExperiment} of
#'     transcript-level unique counts (reads compatible with a single transcript), with
#'     one column per cell. Unique counts are ready for downstream analysis as is.
#'     Collapse them to gene-level counts with \code{\link{transcriptToGeneExpression}},
#'     or use them to cluster cells by gene expression profile (see \code{output = "clusteredEM"}).}
#'     \item{\code{output = "EM"}}{returns a \code{SummarizedExperiment} with transcript
#'     counts estimated using Bambu's EM for each barcode/single cell. This step is
#'     resource intensive, consider using unique counts or the clustered EM.}
#'     \item{\code{output = "clusteredEM"}}{returns a \code{SummarizedExperiment} with
#'     one column per cluster as specified in the cluster argument. \code{colData}
#'     describes the \code{cluster} label and the set of cell barcodes belonging to the
#'     cluster. Quantifying transcript expression for cell clusters provides more stable
#'     EM estimates.}
#' }
#' @seealso \code{\link{bambu}} for the underlying function and the full
#' parameter set; \code{\link{transcriptToGeneExpression}} for preparing gene-level counts
#' for cell clustering and clustered quantification, or for direct use in downstream analysis.
#' @examples
#' ## single-cell example data: two demultiplexed samples (PacBio and ONT)
#' sce.dir <- system.file("extdata", "single_cell", package = "bambu")
#' reads <- file.path(sce.dir, c(
#'     "GIS_cellMix_HepG2-A549-H9-HEYA8_5primeSingleCellcDNAPacBio_Rep1_Run2_demultiplexed_chr9_1_1000000.bam",
#'     "GIS_cellMix_HepG2-A549-H9-HEYA8_5primeSingleCellcDNA_Rep1_Run1_demultiplexed_chr9_1_1000000.bam"))
#' annotations <- system.file("extdata",
#'     "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
#' genome <- system.file("extdata",
#'     "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa",
#'     package = "bambu")
#'
#' ## Each stage is selected with an 'output' preset; its output feeds the next.
#' ## 1. read-class construction from the demultiplexed BAMs
#' readClasses <- bambu.singlecell(reads = reads, output = "readClasses",
#'     annotations = annotations, genome = genome)
#'
#' ## 2. extend the annotations with novel transcripts discovered by Bambu
#' extendedAnnotations <- bambu.singlecell(reads = readClasses,
#'     output = "extendedAnnotations", annotations = annotations)
#'
#' ## 3. assign reads to transcripts (per-sample quantData list)
#' quantData <- bambu.singlecell(reads = readClasses,
#'     output = "quantData", annotations = extendedAnnotations)
#'
#' ## 4. EM-free quantification. Steps 4a-4b give unique counts and gene counts ready
#' ## for downstream analysis; stop here, or continue to step 5 for EM-based
#' ## quantification.
#'
#' ## 4a. transcript-level unique counts (reads compatible with a single transcript).
#' ## These EM-free counts are ready to use downstream; collapse them to gene counts,
#' ## or cluster cells by gene expression to run the clusteredEM step.
#' uniqueCountsSe <- bambu.singlecell(reads = quantData,
#'     output = "uniqueCounts", annotations = extendedAnnotations)
#'
#' ## 4b. gene-level counts: transcriptToGeneExpression() sums the transcript counts
#' ## within each gene to give a gene-level SummarizedExperiment.
#' geneCountsSe <- transcriptToGeneExpression(uniqueCountsSe)
#'
#' ## 5. (optional) EM-based quantification. Steps 5a-5b estimate transcript expression
#' ## with the EM, either at single-cell resolution (5a) or per user-defined cluster
#' ## for more stable estimates (5b).
#'
#' ## 5a. single-cell transcript counts estimated with the EM.
#' EMCountsSe <- bambu.singlecell(reads = quantData,
#'     output = "EM", annotations = extendedAnnotations)
#'
#' ## 5b. per-cluster (pseudobulk) transcript counts estimated with the EM.
#' ## How the clusters are defined is up to the analysis. One common approach is to
#' ## cluster by cell-type using the gene-level counts from step 4 with Seurat (not a
#' ## bambu dependency, so this block runs only when Seurat is available):
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'     seurat <- Seurat::CreateSeuratObject(
#'         counts = assays(geneCountsSe)$counts, min.cells = 1)
#'     seurat <- Seurat::NormalizeData(seurat)
#'     ## this example is small, so cluster on all genes with exact PCA and a small
#'     ## number of components (npcs must stay below the gene count)
#'     seurat <- Seurat::ScaleData(seurat, features = rownames(seurat))
#'     seurat <- Seurat::RunPCA(seurat, features = rownames(seurat),
#'         npcs = 10, approx = FALSE)
#'     seurat <- Seurat::FindNeighbors(seurat, dims = 1:10)
#'     seurat <- Seurat::FindClusters(seurat, resolution = 0.05)
#'
#'     ## Seurat's active.ident is a named vector (names = ids in sampleName_barcode
#'     ## format, values = cluster labels), one of the accepted input data type for
#'     ## `clusters` argument in `bambu.singlecell()`
#'     clusters <- seurat@active.ident
#' }
#'
#' ## In general, clusteredEM accepts the clusters argument as a named vector (e.g.
#' ## seurat@active.ident), a data.frame, or a path to a .csv/.tsv/.txt file. Here we
#' ## illustrate the input data type with a stored .csv:
#' clusters <- file.path(sce.dir,
#'     "clusters_GIS_cellMix_HepG2-A549-H9-HEYA8_5primeSingleCell_multisample_chr9_1_1000000.csv")
#' clusteredEMCountsSe <- bambu.singlecell(reads = quantData,
#'     output = "clusteredEM", annotations = extendedAnnotations,
#'     clusters = clusters)
#' @export
bambu.singlecell <- function(reads, output, annotations = NULL, genome = NULL,
    clusters = NULL, NDR = NULL, stranded = FALSE, ncore = 1, ...) {
    # 'output' is required: it names the single stage to run and return
    if (missing(output) || is.null(output))
        stop("'output' must be one of 'readClasses', 'extendedAnnotations', ",
            "'quantData', 'uniqueCounts', 'EM', or 'clusteredEM'.")
    output <- match.arg(output, c("readClasses", "extendedAnnotations",
        "quantData", "uniqueCounts", "EM", "clusteredEM"))
    # For the EM/clusteredEM stages the quantData list is passed in through 'reads',
    # but bambu() currently takes quantData as a separate argument, so it is rerouted
    # (reads -> quantData) below and defaults to NULL otherwise.
    # TODO: once bambu() accepts a quantData list directly via 'reads', this
    # initializer and the reroute can be removed and 'reads' passed straight through.
    quantData <- NULL
    # validate the inputs each preset needs before setting the stage flags
    # (a genome for BAM input is checked downstream in checkInputs())
    if (output %in% c("readClasses", "extendedAnnotations", "quantData") &&
        is.null(reads))
        stop("output = '", output, "' requires 'reads' to be provided.")
    if (output %in% c("uniqueCounts", "EM", "clusteredEM") && is.null(reads))
        stop("output = '", output, "' requires the per-sample 'quantData' list ",
            "(from output = 'quantData') to be supplied as 'reads'.")
    # annotations is optional for extendedAnnotations (de novo discovery), but
    # required for the assignment and quantification stages
    if (output %in% c("quantData", "uniqueCounts", "EM",
        "clusteredEM") && (is.null(annotations) || length(annotations) == 0))
        stop("output = '", output, "' requires 'annotations' to be provided.")
    if (output %in% c("uniqueCounts", "EM") && !is.null(clusters))
        stop("output = '", output, "' ignores 'clusters'; use output = 'clusteredEM'.")
    if (output == "clusteredEM" && is.null(clusters))
        stop("output = 'clusteredEM' requires 'clusters' to be provided.")
    # each (named) preset returns exactly one object
    if (output == "readClasses") {
        discovery <- FALSE; assignDist <- FALSE; quant <- FALSE
    } else if (output == "extendedAnnotations") {
        discovery <- TRUE;  assignDist <- FALSE; quant <- FALSE
    } else if (output == "quantData") {
        discovery <- FALSE; assignDist <- TRUE;  quant <- FALSE
    } else if (output == "uniqueCounts") {
        return(getUniqueCountsSe(reads, annotations))
    } else {
        discovery <- FALSE; assignDist <- FALSE; quant <- TRUE
        # the quantification stages take the quantData list through 'reads'
        quantData <- reads
        reads <- NULL
    }
    bambu(reads = reads, annotations = annotations, genome = genome, NDR = NDR,
        quantData = quantData,
        opt.singlecell = list(extractBarcodeUMI = TRUE, dedupUMI = TRUE,
            clusters = clusters),
        discovery = discovery, assignDist = assignDist, quant = quant,
        stranded = stranded, ncore = ncore, ...)
}
