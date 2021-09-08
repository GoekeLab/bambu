
#' process reads
#' @param reads path to BAM file(s)
#' @param annotations path to GTF file or TxDb object
#' @param genomeSequence path to FA file or BSgenome object
#' @param readClass.outputDir path to readClass output directory
#' @param yieldSize yieldSize
#' @param bpParameters BioParallel parameter
#' @param stranded stranded
#' @param verbose verbose
#' @importFrom Rsamtools yieldSize BamFileList yieldSize<-
#' @importFrom methods is 
#' @importFrom BiocParallel bplapply
#' @importFrom BiocGenerics basename
#' @noRd
bambu.processReads <- function(reads, annotations, genomeSequence,
    readClass.outputDir=NULL, yieldSize=1000000, bpParameters, 
    stranded=FALSE, verbose=FALSE, min.readCount = 2, fitReadClassModel = T,
    trackReads = FALSE) {
    # ===# create BamFileList object from character #===#
    if (is(reads, "BamFile")) {
        if (!is.null(yieldSize)) {
            yieldSize(reads) <- yieldSize
        } else {
            yieldSize <- yieldSize(reads)
        }
        reads <- BamFileList(reads)
        names(reads) <- tools::file_path_sans_ext(BiocGenerics::basename(reads))
    } else if (is(reads, "BamFileList")) {
        if (!is.null(yieldSize)) {
            yieldSize(reads) <- yieldSize
        } else {
            yieldSize <- min(yieldSize(reads))
        }
    } else if (any(!grepl("\\.bam$", reads))) {
        stop("Bam file is missing from arguments.")
    } else {
        if (is.null(yieldSize)) yieldSize <- NA
        reads <- BamFileList(reads, yieldSize = yieldSize)
        names(reads) <- tools::file_path_sans_ext(BiocGenerics::basename(reads))
    }
    if (!verbose) message("Start generating read class files")
    readClassList <- bplapply(names(reads), function(bamFileName) {
        bambu.processReadsByFile(bam.file = reads[bamFileName],
        genomeSequence = genomeSequence,annotations = annotations,
        readClass.outputDir = readClass.outputDir,
        stranded = stranded, min.readCount = min.readCount, 
        fitReadClassModel = fitReadClassModel, verbose = verbose), 
        trackReads = trackReads)},
        BPPARAM = bpParameters)
    if (!verbose)
        message("Finished generating read classes from genomic alignments.")
    return(readClassList)
}

#' Preprocess bam files and save read class files
#' @inheritParams bambu
#' @importFrom GenomeInfoDb seqlevels seqlevels<- keepSeqlevels
#' @noRd
bambu.processReadsByFile <- function(bam.file, genomeSequence, annotations,
    readClass.outputDir = NULL, stranded = FALSE, min.readCount = 2, 
    fitReadClassModel = TRUE,  verbose = FALSE) {
    readGrgList <- prepareDataFromBam(bam.file[[1]], verbose = verbose, use.names = trackReads)
    seqlevelCheckReadsAnnotation(readGrgList, annotations)
    #check seqlevels for consistency, drop ranges not present in genomeSequence
    refSeqLevels <- seqlevels(genomeSequence)
    if(trackReads) readNames = names(readGrgList)
    unname(readGrgList)
    if (!all(seqlevels(readGrgList) %in% refSeqLevels)) {
        message("not all chromosomes from reads present in reference genome 
            sequence, reads without reference chromosome sequence are dropped")
        refSeqLevels <- intersect(refSeqLevels, seqlevels(readGrgList))
        readGrgList <- keepSeqlevels(readGrgList,
            value =  refSeqLevels,
            pruning.mode = "coarse")
        # reassign Ids after seqlevels are dropped
        mcols(readGrgList)$id <- seq_along(readGrgList) 
    }
    if (!all(seqlevels(annotations) %in% refSeqLevels)) {
    message("not all chromosomes from annotations present in reference genome 
    sequence, annotations without reference chrosomomse sequence are dropped")
    annotations <- keepSeqlevels(annotations,
        value = refSeqLevels,pruning.mode = "coarse")
    }
    # create error and strand corrected junction tables
    unlisted_junctions <- unlistIntrons(readGrgList, use.ids = TRUE)
    uniqueJunctions <- isore.constructJunctionTables(unlisted_junctions, 
        annotations,genomeSequence, stranded = stranded, verbose = verbose)
    # create SE object with reconstructed readClasses
    se <- isore.constructReadClasses(readGrgList, unlisted_junctions, 
        uniqueJunctions, runName = names(bam.file)[1],
        annotations, stranded, verbose)
    if(trackReads) metadata(se)$readNames = readNames
    GenomeInfoDb::seqlevels(se) <- refSeqLevels
    se <- scoreReadClasses(se,genomeSequence, 
                             annotations, 
                             defaultModels = defaultModels,
                             fit = fitReadClassModel,
                             min.readCount = min.readCount,
                             verbose = verbose)
    if (!is.null(readClass.outputDir)) {
        readClassFile <- paste0(readClass.outputDir,names(bam.file),
            "_readClassSe.rds")
        if (file.exists(readClassFile)) {
            show(paste(readClassFile, "exists, will be overwritten"))
            # warning is not printed, use show in addition
            warning(paste(readClassFile, "exists, will be overwritten"))
        } else {
            readClassFile <- BiocFileCache::bfcnew(BiocFileCache::BiocFileCache(
                readClass.outputDir, ask = FALSE),
                paste0(names(bam.file),"_readClassSe"), ext = ".rds")
        }
        saveRDS(se, file = readClassFile)
        se <- readClassFile
    }
    return(se)
}

#' Check seqlevels for reads and annotations
#' @importFrom GenomeInfoDb seqlevels
#' @noRd
seqlevelCheckReadsAnnotation <- function(reads, annotations){
    if (length(intersect(seqlevels(reads),
        seqlevels(annotations))) == 0)
        warning("Warning: no annotations with matching seqlevel styles, 
        all missing chromosomes will use de-novo annotations")
    if (!all(seqlevels(reads) %in% 
        seqlevels(annotations))) 
        message("not all chromosomes present in reference annotations,
            annotations might be incomplete. Please compare objects
            on the same reference")
}