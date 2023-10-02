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
    stranded=FALSE, verbose=FALSE, isoreParameters = setIsoreParameters(NULL),
    lowMemory=FALSE, trackReads = trackReads, fusionMode = fusionMode, 
    demultiplexed = demultiplexed, readGrgListFile = readGrgListFile) {
    genomeSequence <- checkInputSequence(genomeSequence)
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
    min.readCount = isoreParameters[["min.readCount"]]
    fitReadClassModel = isoreParameters[["fitReadClassModel"]]
    defaultModels = isoreParameters[["defaultModels"]]
    returnModel = isoreParameters[["returnModel"]]
    min.exonOverlap = isoreParameters[["min.exonOverlap"]]
    
    ### add ### 
    if (!is.null(readGrgListFile)){
        reads <- as.list(readGrgListFile) 
        names(reads) <- tools::file_path_sans_ext(BiocGenerics::basename(readGrgListFile))
    }
    ### add ### 
    readClassList <- bplapply(names(reads), function(bamFileName){
        bambu.processReadsByFile(bam.file = reads[bamFileName],
        genomeSequence = genomeSequence, annotations = annotations,
        readClass.outputDir = readClass.outputDir, yieldSize = yieldSize, 
        stranded = stranded, min.readCount = min.readCount, 
        fitReadClassModel = fitReadClassModel, min.exonOverlap = min.exonOverlap, 
        defaultModels = defaultModels, returnModel = returnModel, verbose = verbose, 
        lowMemory = lowMemory, trackReads = trackReads, fusionMode = fusionMode, 
        demultiplexed = demultiplexed, readGrgListFile = readGrgListFile)},
        BPPARAM = bpParameters)

    return(readClassList)
}

#' Preprocess bam files and save read class files
#' @inheritParams bambu
#' @importFrom GenomeInfoDb seqlevels seqlevels<- keepSeqlevels
#' @noRd
bambu.processReadsByFile <- function(bam.file, genomeSequence, annotations,
    readClass.outputDir = NULL, yieldSize = NULL, stranded = FALSE, min.readCount = 2, 
    fitReadClassModel = TRUE, min.exonOverlap = 10, defaultModels = NULL, returnModel = FALSE, 
    verbose = FALSE, lowMemory = FALSE, trackReads = FALSE, fusionMode = FALSE, 
    demultiplexed = FALSE, readGrgListFile = NULL) {
    ### add ### 
    if (is.null(readGrgListFile)){
    ### add ###   
      if(verbose) message(names(bam.file)[1])
      readGrgList <- prepareDataFromBam(bam.file[[1]],  yieldSize = yieldSize, verbose = verbose, use.names = trackReads, 
                                        demultiplexed = demultiplexed)
      warnings = c()
      warnings = seqlevelCheckReadsAnnotation(readGrgList, annotations)
      if(verbose & length(warnings) > 0) warning(paste(warnings,collapse = "\n"))
      #check seqlevels for consistency, drop ranges not present in genomeSequence
      refSeqLevels <- seqlevels(genomeSequence)
      if (!all(seqlevels(readGrgList) %in% refSeqLevels)) {
        refSeqLevels <- intersect(refSeqLevels, seqlevels(readGrgList))
        if (!all(seqlevels(annotations) %in% refSeqLevels)&(!(length(annotations)==0))) {
          refSeqLevels <- intersect(refSeqLevels, seqlevels(annotations))
          warningText = paste0("not all chromosomes from annotations present in ", 
                               "reference genome sequence, annotations without reference genomic sequence ",
                               "are dropped")
          warnings = c(warnings, warningText)
          if(verbose) warning(warningText)
          annotations <- keepSeqlevels(annotations, value = refSeqLevels,
                                       pruning.mode = "coarse")
        }
        warningText = paste0("not all chromosomes from reads present in reference ",
                             "genome sequence, reads without reference chromosome sequence are dropped")
        warnings = c(warnings, warningText)
        if(verbose) warning(warningText)
        readGrgList <- keepSeqlevels(readGrgList, value =  refSeqLevels,
                                     pruning.mode = "coarse")
      }
      #removes reads that are outside genome coordinates
      badReads = which(max(end(ranges(readGrgList)))>=
                         seqlengths(genomeSequence)[as.character(getChrFromGrList(readGrgList))])
      if(length(badReads) > 0 ){
        readGrgList = readGrgList[-badReads]
        warningText = paste0(length(badReads), " reads are mapped outside the provided ",
                             "genomic regions. These reads will be dropped. Check you are using the ",
                             "same genome used for the alignment")
        warnings = c(warnings, warningText)
        if(verbose) warning(warningText)
      }
      
      ### add ### 
      # reassign Ids after seqlevels are dropped
      mcols(readGrgList)$id <- seq_along(readGrgList) 
      ### add ###
      
      if(length(readGrgList) == 0) {
        stop("No reads left after filtering.")
      }
      
      ## add ###
      if (isTRUE(demultiplexed)){
        cellBarcodeAssign <- tibble(index = mcols(readGrgList)$id, CB = mcols(readGrgList)$CB) %>% nest(.by = "CB")

        if (!dir.exists("CB")){
          dir.create("CB")
        } else{
          unlink(paste("CB", "*", sep = "/"))
        }
        
        invisible(lapply(seq(nrow(cellBarcodeAssign)),
                  function(x){saveRDS(readGrgList[pull(cellBarcodeAssign$data[[x]])], paste0("CB/", cellBarcodeAssign$CB[[x]],".rds"))}))
      }
      ## add ###
      
    ### add ###    
    } else {
        
        readGrgList <- readRDS(bam.file[[1]])
    
    }
    ### add ###   
  
    # construct read classes for each chromosome seperately 
    if(lowMemory) {se <- lowMemoryConstructReadClasses(readGrgList, genomeSequence, 
                                                      annotations, stranded, verbose,bam.file)
    } else { 
        unlisted_junctions <- unlistIntrons(readGrgList, use.ids = TRUE)
        if(length(unlisted_junctions)==0){
            warningText = paste0("No aligned spliced reads detected!", 
                "Bambu expects spliced reads. If this is intended, ",
                "see Documentation on how to handle single-exon ",
                "transcripts")
            warnings = c(warnings, warningText)
            if (verbose) warning(warningText)
        }
        uniqueJunctions <- isore.constructJunctionTables(unlisted_junctions, 
                                                         annotations,genomeSequence, stranded = stranded, verbose = verbose)
        # create SE object with reconstructed readClasses
        se <- isore.constructReadClasses(readGrgList, unlisted_junctions, 
                                         uniqueJunctions, runName = names(bam.file)[1],
                                         annotations, stranded, verbose)
    }
    metadata(se)$warnings = warnings
    if(trackReads){
        metadata(se)$readNames = names(readGrgList)
        metadata(se)$readId = mcols(readGrgList)$id
    }
    rm(readGrgList)
    
    ### add ### 
    if (is.null(readGrgListFile)){
        GenomeInfoDb::seqlevels(se) <- refSeqLevels   
    }
    ### add ###
    
    # create SE object with reconstructed readClasses
    se <- scoreReadClasses(se, genomeSequence, annotations, 
                             defaultModels = defaultModels,
                             fit = fitReadClassModel,
                             returnModel = returnModel,
                             min.readCount = min.readCount,
                             min.exonOverlap = min.exonOverlap,
                             fusionMode = fusionMode,
                             verbose = verbose)
    
    if (!is.null(readClass.outputDir)) {
        readClassFile <- paste0(readClass.outputDir,names(bam.file),
                                "_readClassSe.rds")
        if (file.exists(readClassFile)) {
            show(paste(readClassFile, "exists, will be overwritten"))
            warning(readClassFile, "exists, will be overwritten")
        } else {
            readClassFile <- BiocFileCache::bfcnew(BiocFileCache::BiocFileCache(
                readClass.outputDir, ask = FALSE),
                paste0(names(bam.file),"_readClassSe"), ext = ".rds")
        }
        saveRDS(se, file = readClassFile)
        se <- readClassFile
    }

    if (!is.null(readGrgListFile)){
      se <- names(se)
    }
    
    return(se)
}

lowMemoryConstructReadClasses <- function(readGrgList, genomeSequence, 
                                          annotations, stranded, verbose,bam.file){
    readGrgList = split(readGrgList, getChrFromGrList(readGrgList))
    se = lapply(names(readGrgList),FUN = function(i){
        if(length(readGrgList[[i]]) == 0) return(NULL)
        # create error and strand corrected junction tables
        unlisted_junctions <- unlistIntrons(readGrgList[[i]], use.ids = TRUE)
        uniqueJunctions <- isore.constructJunctionTables(unlisted_junctions, 
                                                         annotations,genomeSequence, stranded = stranded, verbose = verbose)
        se.temp <- isore.constructReadClasses(readGrgList[[i]], 
                                              unlisted_junctions, uniqueJunctions, runName = names(bam.file)[1],
                                              annotations, stranded, verbose)
        return(se.temp)
    })
    se = se[!sapply(se, FUN = is.null)]
    se = do.call("rbind",se)
    rownames(se) = paste("rc", seq_len(nrow(se)), sep = ".")
    return(se)
}

#' Check seqlevels for reads and annotations
#' @importFrom GenomeInfoDb seqlevels
#' @noRd
seqlevelCheckReadsAnnotation <- function(reads, annotations){
    warnings = c()
    if (length(intersect(seqlevels(reads),
                         seqlevels(annotations))) == 0)
        warnings = c(warnings, paste0("no annotations with matching seqlevel styles, ",
        "all missing chromosomes will use de-novo annotations"))
    if (!all(seqlevels(reads) %in% 
             seqlevels(annotations))) 
        warnings = c(warnings, paste0("not all chromosomes present in reference annotations, ",
            "annotations might be incomplete. Please compare objects ",
            "on the same reference"))
    return(warnings)
}
