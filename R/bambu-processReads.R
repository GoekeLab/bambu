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
    stranded=FALSE, verbose=FALSE, discoveryParameters = setDiscoveryParameters(NULL),
    processByChromosome = FALSE, processByBam = TRUE, trackReads = trackReads, fusionMode = fusionMode, 
    extractBarcodeUMI = FALSE, dedupUMI = FALSE) {
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
    min.readCount <- discoveryParameters[["min.readCount"]]
    fitReadClassModel <- discoveryParameters[["fitReadClassModel"]]
    defaultModels <- discoveryParameters[["defaultModels"]]
    returnModel <- discoveryParameters[["returnModel"]]
    min.exonOverlap <- discoveryParameters[["min.exonOverlap"]]

    if(processByBam){
        readClassList <- bplapply(seq_along(reads), function(i) {
            bambu.processReadsByFile(bam.file = reads[i],
            genomeSequence = genomeSequence,annotations = annotations,
            stranded = stranded, min.readCount = min.readCount, 
            fitReadClassModel = fitReadClassModel, min.exonOverlap = min.exonOverlap, 
            defaultModels = defaultModels, returnModel = returnModel, verbose = verbose, 
            processByChromosome = processByChromosome, trackReads = trackReads, fusionMode = fusionMode, 
            extractBarcodeUMI = extractBarcodeUMI, dedupUMI = dedupUMI, index = 1)},
            BPPARAM = bpParameters)
        names(readClassList) <- names(reads)
    } else {
        readGrgList <- bplapply(seq_along(reads), function(i) {
            bambu.readsByFile(bam.file = reads[i],
            genomeSequence = genomeSequence,annotations = annotations,
            stranded = stranded, min.readCount = min.readCount,
            fitReadClassModel = fitReadClassModel, min.exonOverlap = min.exonOverlap,
            defaultModels = defaultModels, returnModel = returnModel, verbose = verbose,
            trackReads = trackReads, fusionMode = fusionMode,
            extractBarcodeUMI = extractBarcodeUMI, dedupUMI = dedupUMI, index = i)},
            BPPARAM = bpParameters)
        for(i in seq_along(readGrgList)){
            if(extractBarcodeUMI){
                mcols(readGrgList[[i]])$CB <- paste0(names(reads)[i], '_', mcols(readGrgList[[i]])$CB)
            } else{
                mcols(readGrgList[[i]])$CB <- names(reads)[i]
            }
            
            mcols(readGrgList[[i]])$CB <- as.factor(mcols(readGrgList[[i]])$CB)
            
        }
        readGrgList <- do.call(c, readGrgList)    
        mcols(readGrgList)$id <- seq_along(readGrgList) 
        if(extractBarcodeUMI){
          mcols(readGrgList)$sampleID <- as.numeric(mcols(readGrgList)$CB)
        } else {
          mcols(readGrgList)$sampleID <- i
        }
        readClassList <- constructReadClasses(readGrgList, genomeSequence = genomeSequence,annotations = annotations,
            stranded = stranded, min.readCount = min.readCount,
            fitReadClassModel = fitReadClassModel, min.exonOverlap = min.exonOverlap,
            defaultModels = defaultModels, returnModel = returnModel, verbose = verbose,
            processByChromosome = processByChromosome, trackReads = trackReads, fusionMode = fusionMode)
        metadata(readClassList)$samples <- names(reads)
        metadata(readClassList)$sampleNames <- names(reads)
        if(extractBarcodeUMI) metadata(readClassList)$samples <- levels(mcols(readGrgList)$CB)
        readClassList <- list(readClassList)
    }
        
    if (!is.null(readClass.outputDir)) {
        for(i in seq_along(readClassList)){
            readClassFile <- "combinedSamples"
            readClassFile <- BiocFileCache::bfcnew(BiocFileCache::BiocFileCache(
                readClass.outputDir, ask = FALSE),
                paste0(readClassFile,"_readClassSe"), ext = ".rds")
            saveRDS(readClassList[[i]], file = readClassFile)
            readClassList[[i]] <- readClassFile
        }
    }
    #TODO don't output list, current there because discovery needs it
    return(readClassList)
}

#' Preprocess bam files and save read class files
#' @inheritParams bambu
#' @importFrom GenomeInfoDb seqlevels seqlevels<- keepSeqlevels
#' @noRd
bambu.processReadsByFile <- function(bam.file, genomeSequence, annotations,
    yieldSize = NULL, stranded = FALSE, min.readCount = 2, 
    fitReadClassModel = TRUE, min.exonOverlap = 10, defaultModels = NULL, returnModel = FALSE, 
    verbose = FALSE, processByChromosome = FALSE, trackReads = FALSE, fusionMode = FALSE,
    extractBarcodeUMI = FALSE, dedupUMI = FALSE, index = 0) {
    if(verbose) message(names(bam.file)[1])
    readGrgList <- prepareDataFromBam(bam.file[[1]], verbose = verbose, yieldSize = yieldSize, use.names = trackReads, extractBarcodeUMI = extractBarcodeUMI, dedupUMI = dedupUMI)
    if(verbose) message(paste0("Number of alignments/reads: ",length(readGrgList)))
    warnings <- c()
    warnings <- seqlevelCheckReadsAnnotation(readGrgList, annotations)
    if(verbose & length(warnings) > 0) warning(paste(warnings,collapse = "\n"))
    #check seqlevels for consistency, drop ranges not present in genomeSequence
    refSeqLevels <- seqlevels(genomeSequence)
    if (!all(seqlevels(readGrgList) %in% refSeqLevels)) {
        refSeqLevels <- intersect(refSeqLevels, seqlevels(readGrgList))
        if (!all(seqlevels(annotations) %in% refSeqLevels)&(!(length(annotations)==0))) {
            refSeqLevels <- intersect(refSeqLevels, seqlevels(annotations))
            warningText <- paste0("not all chromosomes from annotations present in ", 
            "reference genome sequence, annotations without reference genomic sequence ",
            "are dropped")
            warnings <- c(warnings, warningText)
            if(verbose) warning(warningText)
            annotations <- keepSeqlevels(annotations, value = refSeqLevels,
                                         pruning.mode = "coarse")
        }
        warningText <- paste0("not all chromosomes from reads present in reference ",
        "genome sequence, reads without reference chromosome sequence are dropped")
        warnings <- c(warnings, warningText)
        if(verbose) warning(warningText)
        readGrgList <- keepSeqlevels(readGrgList, value =  refSeqLevels,
                                     pruning.mode = "coarse")
        # reassign Ids after seqlevels are dropped
        mcols(readGrgList)$id <- seq_along(readGrgList) 
    }
    #removes reads that are outside genome coordinates
    badReads <- which(max(end(ranges(readGrgList)))>
                         seqlengths(genomeSequence)[as.character(getChrFromGrList(readGrgList))])
    if(length(badReads) > 0 ){
        readGrgList <- readGrgList[-badReads]
        warningText <- paste0(length(badReads), " reads are mapped outside the provided ",
                       "genomic regions. These reads will be dropped. Check you are using the ",
                       "same genome used for the alignment")
        warnings <- c(warnings, warningText)
        if(verbose) warning(warningText)
    }
    if(length(readGrgList) == 0)
        stop("No reads left after filtering.")

    mcols(readGrgList)$id <- seq_along(readGrgList) 

    if(extractBarcodeUMI){ 
        mcols(readGrgList)$columnID <- as.numeric(mcols(readGrgList)$CB)
    } else {
        mcols(readGrgList)$columnID <- index
    }
        
    runName <- names(bam.file)[1]
    # construct read classes for each chromosome seperately 
    if(processByChromosome){
        se <- lowMemoryConstructReadClasses(readGrgList, genomeSequence, 
                                                      annotations, stranded, verbose, bam.file)
    } else{
        unlisted_junctions <- unlistIntrons(readGrgList, use.ids = TRUE)
        uniqueJunctions <- isore.constructJunctionTables(unlisted_junctions, 
                                                         annotations,genomeSequence, stranded = stranded, verbose = verbose)
        se <- isore.constructReadClasses(readGrgList, 
                                              unlisted_junctions, uniqueJunctions, runName = runName,
                                              annotations, stranded, verbose)

    }

    metadata(se)$warnings <- warnings
    if(trackReads){
        metadata(se)$readNames <- names(readGrgList)
        metadata(se)$readId <- mcols(readGrgList)$id
    }
    refSeqLevels <- seqlevels(genomeSequence)
    GenomeInfoDb::seqlevels(se) <- refSeqLevels
    # create SE object with reconstructed readClasses
    se <- scoreReadClasses(se, genomeSequence, annotations, 
                             defaultModels = defaultModels,
                             fit = fitReadClassModel,
                             returnModel = returnModel,
                             min.readCount = min.readCount,
                             min.exonOverlap = min.exonOverlap,
                             fusionMode = fusionMode,
                             verbose = verbose)

    if (extractBarcodeUMI) {
        barcodes <- levels(mcols(readGrgList)$CB)
        metadata(se)$sampleData <- tibble(
          id = paste(names(bam.file)[1], barcodes, sep = '_'),
          sampleName = names(bam.file)[1],
          barcode = barcodes
        )
    } else{
        metadata(se)$sampleData <- tibble(
          id = names(bam.file)[1],
          sampleName = names(bam.file)[1]
        )
    }

    return(se)
}

#' Preprocess bam files and save read class files
#' @inheritParams bambu
#' @importFrom GenomeInfoDb seqlevels seqlevels<- keepSeqlevels
#' @noRd
bambu.readsByFile <- function(bam.file, genomeSequence, annotations,
    yieldSize = NULL, stranded = FALSE, min.readCount = 2, 
    fitReadClassModel = TRUE, min.exonOverlap = 10, defaultModels = NULL, returnModel = FALSE, 
    verbose = FALSE, trackReads = FALSE, fusionMode = FALSE,
    extractBarcodeUMI = FALSE, dedupUMI = FALSE, index = 0) {
    readGrgList <- prepareDataFromBam(bam.file[[1]], verbose = verbose, yieldSize = yieldSize, use.names = trackReads, extractBarcodeUMI = extractBarcodeUMI, dedupUMI = dedupUMI)

    if(verbose) message("Number of alignments/reads: ",length(readGrgList))
    
    warnings <- c()
    warnings <- seqlevelCheckReadsAnnotation(readGrgList, annotations)
    
    if(verbose & length(warnings) > 0) warning(paste(warnings,collapse = "\n"))
    #check seqlevels for consistency, drop ranges not present in genomeSequence
    refSeqLevels <- seqlevels(genomeSequence)
    if (!all(seqlevels(readGrgList) %in% refSeqLevels)) {
        refSeqLevels <- intersect(refSeqLevels, seqlevels(readGrgList))
        if (!all(seqlevels(annotations) %in% refSeqLevels)&(!(length(annotations)==0))) {
          refSeqLevels <- intersect(refSeqLevels, seqlevels(annotations))
          warningText <- paste0("not all chromosomes from annotations present in ", 
                               "reference genome sequence, annotations without reference genomic sequence ",
                               "are dropped")
          warnings <- c(warnings, warningText)
          if(verbose) warning(warningText)
          annotations <- keepSeqlevels(annotations, value = refSeqLevels,
                                       pruning.mode = "coarse")
        }
        warningText <- paste0("not all chromosomes from reads present in reference ",
                             "genome sequence, reads without reference chromosome sequence are dropped")
        warnings <- c(warnings, warningText)
        if(verbose) warning(warningText)
        readGrgList <- keepSeqlevels(readGrgList, value =  refSeqLevels,
                                     pruning.mode = "coarse")
        # reassign Ids after seqlevels are dropped
        mcols(readGrgList)$id <- seq_along(readGrgList) 
    }
    #removes reads that are outside genome coordinates
    badReads <- which(max(end(ranges(readGrgList)))>=
                         seqlengths(genomeSequence)[as.character(getChrFromGrList(readGrgList))])
      if(length(badReads) > 0 ){
        readGrgList <- readGrgList[-badReads]
        warningText <- paste0(length(badReads), " reads are mapped outside the provided ",
                             "genomic regions. These reads will be dropped. Check you are using the ",
                             "same genome used for the alignment")
        warnings <- c(warnings, warningText)
        if(verbose) warning(warningText)
      }
      
      ### add ### 
      # reassign Ids after seqlevels are dropped
      mcols(readGrgList)$id <- seq_along(readGrgList) 
      ### add ###
      if(verbose) message("Number of post-filter alignments/reads: ",length(readGrgList))
      if(length(readGrgList) == 0)
        stop("No reads left after filtering.")
      
      ## add ###
      #if (extractBarcodeUMI){
      #  cellBarcodeAssign <- tibble(index = mcols(readGrgList)$id, CB = mcols(readGrgList)$CB) %>% nest(.by = "CB")

        # if (!dir.exists("CB")){
        #   dir.create("CB")
        # } else{
        #   unlink(paste("CB", "*", sep = "/"))
        # }
        
        # invisible(lapply(seq(nrow(cellBarcodeAssign)),
        #           function(x){saveRDS(readGrgList[pull(cellBarcodeAssign$data[[x]])], paste0("CB/", cellBarcodeAssign$CB[[x]],".rds"))}))
      #} 
    return(readGrgList)
}

#' Construct read classes
#' @noRd
constructReadClasses <- function(readGrgList, genomeSequence, annotations,
    stranded = FALSE, min.readCount = 2, 
    fitReadClassModel = TRUE, min.exonOverlap = 10, defaultModels = NULL, returnModel = FALSE, 
    verbose = FALSE, processByChromosome = FALSE, trackReads = FALSE, fusionMode = FALSE, runName = "sample"){
    
    if(processByChromosome){
        # construct read classes for each chromosome seperately 
        se <- lowMemoryConstructReadClasses(readGrgList, genomeSequence, 
                                            annotations, stranded, verbose, runName, fusionMode)
    } else{
        unlisted_junctions <- unlistIntrons(readGrgList, use.ids = TRUE)
        uniqueJunctions <- isore.constructJunctionTables(unlisted_junctions, 
                                                         annotations,genomeSequence, stranded = stranded, verbose = verbose)
        se <- isore.constructReadClasses(readGrgList, 
                                              unlisted_junctions, uniqueJunctions, runName = runName,
                                              annotations, stranded, verbose)

    }
    metadata(se)$warnings <- warnings
    if(trackReads){
        metadata(se)$readNames <- names(readGrgList)
        metadata(se)$readId <- mcols(readGrgList)$id
    }
    rm(readGrgList)
    refSeqLevels <- seqlevels(genomeSequence)
    GenomeInfoDb::seqlevels(se) <- refSeqLevels
    # create SE object with reconstructed readClasses
    se <- scoreReadClasses(se, genomeSequence, annotations, 
                             defaultModels = defaultModels,
                             fit = fitReadClassModel,
                             returnModel = returnModel,
                             min.readCount = min.readCount,
                             min.exonOverlap = min.exonOverlap,
                             fusionMode = fusionMode,
                             verbose = verbose)
    return(se)
}


#' Low memory mode for construct read classes (processByChromosome)
#' @noRd
lowMemoryConstructReadClasses <- function(readGrgList, genomeSequence, 
                                          annotations, stranded, verbose, bam.file, fusionMode = FALSE){
    if(fusionMode){
        readGrgList <- list(readGrgList)
        names(readGrgList) <- c("fusion")
    } else{
        readGrgList <- split(readGrgList, getChrFromGrList(readGrgList))
    }
    runName <- names(bam.file)[1]
    se <- lapply(names(readGrgList),FUN = function(i){
        if(length(readGrgList[[i]]) == 0) return(NULL)
        # create error and strand corrected junction tables
        unlisted_junctions <- unlistIntrons(readGrgList[[i]], use.ids = TRUE)
        uniqueJunctions <- isore.constructJunctionTables(unlisted_junctions, 
                                                         annotations,genomeSequence, stranded = stranded, verbose = verbose)
        se.temp <- isore.constructReadClasses(readGrgList[[i]], 
                                              unlisted_junctions, uniqueJunctions, runName = runName,
                                              annotations, stranded, verbose)
        return(se.temp)
    })
    se <- se[!sapply(se, FUN = is.null)]
    se <- do.call("rbind",se)
    rownames(se) <- paste("rc", seq_len(nrow(se)), sep = ".")
    return(se)
}

#' Check seqlevels for reads and annotations
#' @importFrom GenomeInfoDb seqlevels
#' @noRd
seqlevelCheckReadsAnnotation <- function(reads, annotations){
    warnings <- c()
    if (length(intersect(seqlevels(reads),
                         seqlevels(annotations))) == 0)
        warnings <- c(warnings, paste0("no annotations with matching seqlevel styles, ",
        "all missing chromosomes will use de-novo annotations"))
    if (!all(seqlevels(reads) %in% 
             seqlevels(annotations))) 
        warnings <- c(warnings, paste0("not all chromosomes present in reference annotations, ",
            "annotations might be incomplete. Please compare objects ",
            "on the same reference"))
    return(warnings)
}


#' Split read class files
#' @importFrom dplyr Matrix
#' @noRd
splitReadClassFiles = function(readClassFile){
    distTable <- metadata(metadata(readClassFile)$readClassDist)$distTable  
    eqClasses <- distTable %>% group_by(eqClassById) %>% 
        distinct(eqClassById, readCount,GENEID, totalWidth, firstExonWidth, .keep_all = TRUE)
    eqClasses$columnIds <- rowData(readClassFile)$columnIds[match(eqClasses$readClassId, rownames(readClassFile))]
    eqClasses <- eqClasses %>% summarise(nobs = sum(readCount),
                                                columnIds = list(unlist(columnIds)))
    counts.table <- tableFunction(eqClasses$columnIds)
    metadata(readClassFile)$columnIds <- lapply(counts.table, function(x) as.numeric(names(x)))
    metadata(readClassFile)$columnCounts <- lapply(counts.table, function(x) as.numeric(x))
    counts <- sparseMatrix(
        i = rep(seq_along(counts.table), lengths(counts.table)),
        j = as.numeric(names(unlist(counts.table))),
        x = unlist(counts.table),
        dims = c(nrow(eqClasses), nrow(metadata(readClassFile)$sampleData)))
    #incompatible counts
    distTable <- metadata(metadata(readClassFile)$readClassDist)$distTable.incompatible
    if(nrow(distTable)==0) {
        counts.incompatible <- sparseMatrix(i= integer(0), j = integer(0), x = numeric(0),
        dims = c(0, length(metadata(readClassFile)$sampleData$id)))
        rownames(counts.incompatible) <- character(0)
    } else{
        distTable$columnIds <- rowData(readClassFile)$columnIds[match(distTable$readClassId, rownames(readClassFile))]
        distTable <- distTable %>% group_by(GENEID.i) %>% summarise(counts = sum(readCount),
                    columnIds = list(unlist(columnIds)))
        counts.table <- lapply(distTable$columnIds, FUN = function(x){table(x)})
        counts.incompatible <- sparseMatrix(
            i = rep(seq_along(counts.table), lengths(counts.table)),
            j = as.numeric(names(unlist(counts.table))),
            x = unlist(counts.table),
            dims = c(nrow(distTable), length(metadata(readClassFile)$sampleData$id)))
        colnames(counts.incompatible) <- metadata(readClassFile)$sampleData$id
        rownames(counts.incompatible) <- distTable$GENEID.i 
    }
    colnames(counts) <- metadata(readClassFile)$sampleData$id
    metadata(readClassFile)$eqClassById <- eqClasses$eqClassById
    #rownames(counts) = eqClasses$eqClassById
    metadata(readClassFile)$incompatibleCountMatrix <- counts.incompatible  
    return(readClassFile)
}


#' Split read class files by RC
#' @importFrom Matrix
#' @noRd
splitReadClassFilesByRC <- function(readClassFile){
    counts.table <- tableFunction(rowData(readClassFile)$columnIds)
    counts <- sparseMatrix(
        i = rep(seq_along(counts.table), lengths(counts.table)),
        j = as.numeric(names(unlist(counts.table))),
        x = unlist(counts.table),
        dims = c(nrow(readClassFile), length(metadata(readClassFile)$samples)))
    return(counts)
}

#' table sample IDs list column
#' @noRd
tableFunction <- function(xList){
    return(lapply(xList, function(x) table(x)))
}
