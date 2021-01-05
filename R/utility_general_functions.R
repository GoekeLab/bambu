#' @useDynLib bambu, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @noRd
.onUnload <- function(libpath) {
    library.dynam.unload("bambu", libpath)
}


#' process reads
#' @param reads path to BAM file(s)
#' @param annotations path to GTF file or TxDb object
#' @param genomeSequence path to FA file or BSgenome object
#' @param readClass.outputDir path to readClass output directory
#' @param yieldSize yieldSize
#' @param bpParameters BioParallel parameter
#' @param stranded stranded
#' @param ncore ncore
#' @param verbose verbose
#' @noRd
processReads <- function(reads, readClass.file, annotations, genomeSequence,
    readClass.outputDir, yieldSize, bpParameters, stranded, ncore, verbose) {
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

        if (!verbose) message("Start generating read class files")
        readClassList <- BiocParallel::bplapply(names(reads),
            function(bamFileName) {
            bambu.constructReadClass(bam.file = reads[bamFileName],
                readClass.outputDir = readClass.outputDir,
                genomeSequence = genomeSequence,annotations = annotations,
                stranded = stranded,ncore = ncore,verbose = verbose)},
        BPPARAM = bpParameters)
        if (!verbose)
            message("Finished generating read classes from genomic alignments.")

    return(readClassList)
}

#' @noRd
helpFun <- function(chr, chrRanges, bamFile) {
    return(GenomicAlignments::grglist(GenomicAlignments::readGAlignments(
        file = bamFile,
        param = Rsamtools::ScanBamParam(
            flag = Rsamtools::scanBamFlag(isSecondaryAlignment = FALSE),
            which = chrRanges[chr]),
        use.names = FALSE)))
}


