
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
bambu.processReads <- function(reads, readClass.file, annotations, genomeSequence,
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
                                            bambu.processReadsByFile(bam.file = reads[bamFileName],
                                                                     readClass.outputDir = readClass.outputDir,
                                                                     genomeSequence = genomeSequence,annotations = annotations,
                                                                     stranded = stranded,ncore = ncore,verbose = verbose)},
                                          BPPARAM = bpParameters)
  if (!verbose)
    message("Finished generating read classes from genomic alignments.")
  
  return(readClassList)
}

#' Preprocess bam files and save read class files
#' @inheritParams bambu
#' @noRd
bambu.processReadsByFile <- function(bam.file, genomeSequence, annotations,
                                     readClass.outputDir = NULL, stranded = FALSE, ncore = 1, verbose = FALSE) {
  readGrgList <- prepareDataFromBam(bam.file[[1]], ncore = ncore,
                                    verbose = verbose)
  if (length(intersect(GenomeInfoDb::seqlevels(readGrgList),
                       GenomeInfoDb::seqlevels(annotations))) == 0)
    stop("Error: please provide annotation with matched seqlevel styles.")
  se <- isore.constructReadClasses(
    readGrgList = readGrgList,
    runName = names(bam.file)[1],
    annotationGrangesList = annotations,
    genomeSequence = genomeSequence,
    stranded = stranded,
    ncore = ncore,
    verbose = verbose)
  GenomeInfoDb::seqlevels(se) <- unique(c(GenomeInfoDb::seqlevels(se),
                                          GenomeInfoDb::seqlevels(annotations)))
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
