
#' process reads
#' @param reads path to BAM file(s)
#' @param annotations path to GTF file or TxDb object
#' @param genomeSequence path to FA file or BSgenome object
#' @param readClass.outputDir path to readClass output directory
#' @param yieldSize yieldSize
#' @param bpParameters BioParallel parameter
#' @param stranded stranded
#' @param verbose verbose
#' @noRd
bambu.processReads <- function(reads, annotations, genomeSequence,
    readClass.outputDir=NULL, yieldSize=1000000, bpParameters, 
    stranded=FALSE, verbose=FALSE) {
    # ===# create BamFileList object from character #===#
    if (methods::is(reads, "BamFile")) {
        if (!is.null(yieldSize)) {
            Rsamtools::yieldSize(reads) <- yieldSize
        } else {
            yieldSize <- Rsamtools::yieldSize(reads)
        }
        reads <- Rsamtools::BamFileList(reads)
        names(reads) <- tools::file_path_sans_ext(BiocGenerics::basename(reads))
    } else if (methods::is(reads, "BamFileList")) {
        if (!is.null(yieldSize)) {
            Rsamtools::yieldSize(reads) <- yieldSize
        } else {
            yieldSize <- min(Rsamtools::yieldSize(reads))
        }
    } else if (any(!grepl("\\.bam$", reads))) {
        stop("Bam file is missing from arguments.")
    } else {
        if (is.null(yieldSize)) yieldSize <- NA
        reads <- Rsamtools::BamFileList(reads, yieldSize = yieldSize)
        names(reads) <- tools::file_path_sans_ext(BiocGenerics::basename(reads))
    }
    genomeSequence <- checkInputSequence(genomeSequence)
    if (!verbose) message("Start generating read class files")
    readClassList <- 
        BiocParallel::bplapply(names(reads), function(bamFileName) {
        bambu.processReadsByFile(bam.file = reads[bamFileName],
        readClass.outputDir = readClass.outputDir,
        genomeSequence = genomeSequence,annotations = annotations,
        stranded = stranded,verbose = verbose)},
        BPPARAM = bpParameters)
    if (!verbose)
        message("Finished generating read classes from genomic alignments.")
    return(readClassList)
}

#' Preprocess bam files and save read class files
#' @inheritParams bambu
#' @noRd
bambu.processReadsByFile <- function(bam.file, genomeSequence, annotations,
    readClass.outputDir = NULL, stranded = FALSE, verbose = FALSE) {
    readGrgList <- prepareDataFromBam(bam.file[[1]], verbose = verbose)
    seqlevelCheckReadsAnnotation(readGrgList, annotations)
    #check seqlevels for consistency, drop ranges not present in genomeSequence
    refSeqLevels <-  GenomeInfoDb::seqlevels(genomeSequence)
    if (!all(GenomeInfoDb::seqlevels(readGrgList) %in% refSeqLevels)) {
        message("not all chromosomes from reads present in reference genome 
            sequence, reads without reference chromosome sequence are dropped")
        readGrgList <- GenomeInfoDb::keepSeqlevels(readGrgList,
            value =  refSeqLevels,
            pruning.mode = "coarse")
        # reassign Ids after seqlevels are dropped
        mcols(readGrgList)$id <- seq_along(readGrgList) 
    }
    if (!all(GenomeInfoDb::seqlevels(annotations) %in% refSeqLevels)) {
    message("not all chromosomes from annotations present in reference genome 
    sequence, annotations without reference chrosomomse sequence are dropped")
    annotations <- GenomeInfoDb::keepSeqlevels(annotations,
        value = refSeqLevels,pruning.mode = "coarse")
    }
    # create error and strand corrected junction tables
    unlisted_junctions <- unlistIntrons(readGrgList, use.ids = TRUE)
    uniqueJunctions <- isore.constructJunctionTables(unlisted_junctions, 
        annotations,genomeSequence, stranded = stranded, verbose = verbose)
    # create SE object with reconstructed readClasses
    isore.constructReadClassesOutput <- 
        isore.constructReadClasses(readGrgList, unlisted_junctions, 
        uniqueJunctions, runName = names(bam.file)[1],
        annotations, stranded, verbose)
    se = isore.constructReadClassesOutput$se
    readGrgList = isore.constructReadClassesOutput$readGrgList
    GenomeInfoDb::seqlevels(se) <- refSeqLevels
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

    #txRange starts here!
    se = txrange.filterReadClasses(se, readGrgList, genomeSequence, annotations)

    return(se)

}

#' Check seqlevels for reads and annotations
#' @noRd
seqlevelCheckReadsAnnotation <- function(reads, annotations){
    if (length(intersect(GenomeInfoDb::seqlevels(reads),
        GenomeInfoDb::seqlevels(annotations))) == 0)
        stop("Error: please provide annotation with matched seqlevel styles.")
    if (!all(GenomeInfoDb::seqlevels(reads) %in% 
        GenomeInfoDb::seqlevels(annotations))) 
        message("not all chromosomes present in reference annotations,
            annotations might be incomplete. Please compare objects
            on the same reference")
}

  
#' Function to create a object that can be queried by getSeq
#' Either from fa file, or BSGenome object
#' @importFrom BiocParallel bppram bpvec
#' @noRd
checkInputSequence <- function(genomeSequence) {
    if (is.null(genomeSequence)) stop("Reference genome sequence is missing,
        please provide fasta file or BSgenome name, see available.genomes()")
    if (methods::is(genomeSequence, "character")) {
        if (grepl(".fa", genomeSequence)) {
            if (.Platform$OS.type == "windows") {
                genomeSequence <- Biostrings::readDNAStringSet(genomeSequence)
                newlevels <- unlist(lapply(strsplit(names(genomeSequence)," "),
                    "[[", 1))
                names(genomeSequence) <- newlevels
            } else {
                genomeSequence <- Rsamtools::FaFile(genomeSequence)
            }
        } else {
            genomeSequence <- BSgenome::getBSgenome(genomeSequence)
        }
    }
    return(genomeSequence)
}
