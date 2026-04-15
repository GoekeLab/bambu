
## Functions to set basic parameters and check inputs
#' setBiocParallelParameters
#' @importFrom BiocParallel bpparam
#' @noRd
setBiocParallelParameters <- function(reads, ncore, verbose, demultiplexed){
    bpParameters <- bpparam()
    #===# set parallel options: otherwise use parallel to distribute samples
    # when demultiplexed is FALSE, isFALSE(demultiplexed) is TRUE
    bpParameters$workers <- ifelse(length(reads) == 1 & isFALSE(demultiplexed), 1, ncore)
    bpParameters$progressbar <- ifelse(length(reads) > 1 & !verbose, TRUE, FALSE)
    return(bpParameters)
}


#' setIsoreparameters
#' @noRd
setIsoreParameters <- function(isoreParameters){
    # ===# set default controlling parameters for isoform reconstruction  #===#
    isoreParameters.default <- list(
        remove.subsetTx = TRUE, 
        min.readCount = 2,
        min.readFractionByGene = 0.05,
        min.sampleNumber = 1,
        min.exonOverlap = 10,
        min.txScore.multiExon = 0,
        min.txScore.singleExon = 1,
        fitReadClassModel = TRUE,
        defaultModels = defaultModels,
        returnModel = FALSE,
        baselineFDR = 0.1,
        min.readFractionByEqClass = 0,
        prefix = "Bambu") 
    isoreParameters <- 
        updateParameters(isoreParameters, isoreParameters.default)
    return(isoreParameters)
}


#' setRcAssignmentParameters
#' @noRd
setRcAssignmentParameters <- function(rcAssignmentParameters){
    rcAssignmentParameters.default <- list(
        min.exonDistance = 35,
        min.primarySecondaryDist = 5,
        min.primarySecondaryDistStartEnd1 = 5, # for extending annotations
        min.primarySecondaryDistStartEnd2 = 5) # for read class assignment
    rcAssignmentParameters <-
        updateParameters(rcAssignmentParameters, rcAssignmentParameters.default)
    return(rcAssignmentParameters)
}


#' setEmParameters
#' @noRd
setEmParameters <- function(emParameters){
    emParameters.default <- list(degradationBias = TRUE, maxiter = 10000, 
        conv = 10^(-2), minvalue = 10^(-8), sig.digit = 5)
    emParameters <- updateParameters(emParameters, emParameters.default)
    return(emParameters)
}

#' check parameters for isore and em
#' @param Parameters parameters inputted by user
#' @param Parameters.default default parameters
#' @noRd
updateParameters <- function(Parameters, Parameters.default) {
    if (!is.null(Parameters)) {
        for (i in names(Parameters)) {
            if(!(i %in% names(Parameters.default))) message("Setting parameter that does not exist. Check the spelling - ", i)
            Parameters.default[[i]] <- Parameters[[i]]
        }
    }
    Parameters <- Parameters.default
    return(Parameters)
}

#' check valid inputs
#' @param annotations path to GTF file or TxDb object
#' @param reads path to BAM file(s)
#' @param readClass.file path to readClass file(s)
#' @param readClass.outputDir path to readClass output directory
#' @importFrom methods is
#' @noRd
checkInputs <- function(annotations, reads, readClass.outputDir, genomeSequence, 
                        discovery, sampleNames, sampleData, quantData){
    # ===# Check annotation inputs #===#
    if (!is.null(annotations)) {
        if (is(annotations, "CompressedGRangesList")) {
            ## check if annotations is as expected
            if (!all(c("TXNAME", "GENEID", "txid","eqClassById") %in% 
                     colnames(mcols(annotations)))) 
                stop("The annotations is not properly prepared.\nPlease 
                    see ?prepareAnnotations for help")
            if(anyDuplicated(mcols(annotations)$TXNAME)) {
                warning('Annotations contain duplicated transcript/gene names
                        Please re-create your annotation object')
            }
        } 
        else if (is(annotations, "TxDb") | grepl(".gtf$", annotations)) {
            if (grepl(".gtf$", annotations)) 
                message("If you are running bambu multiple times we recommend ",
                "processing your annotation file first with ",
                "annotations = prepareAnnotations(gtf.file)")
            annotations <- prepareAnnotations(annotations)
        } else {
            stop("The annotations is not a GRangesList object a TxDb or a path to a .gtf.")
        }
        if(discovery & (any(grepl("^BambuGene", names(annotations))) | 
            any(grepl("^BambuTx", mcols(annotations)$TXNAME)))){
                message("Detected Bambu derived annotations in the annotations. ", 
                "Set a new prefix with opt.discovery(list(prefix='newPrefix')) ",
                "to prevent ambigious id assignment.")
        }
    } else {
        stop("Annotations is missing.")
    }
    # ===# Check whether provided readClass.outputDir exists  #===#
    if (!is.null(readClass.outputDir)) {
        if (!dir.exists(readClass.outputDir)) 
            stop("output folder does not exist")
    }
    if(!is.null(reads)){
        if (is(reads, "BamFileList")){
            if(is.null(genomeSequence)){
                stop("A genome must be provided when running bambu from bam files")
            }
        } else{
            # ===# Check whether provided read files are all in the same format (.bam or .rds) #===#
            isRDSs <- all(sapply(reads, class)=="RangedSummarizedExperiment")
            # there is a bug here, when reads is NULL, isRDSs == TRUE
            if(!isRDSs){
                if (!all(grepl(".bam$", reads)) & !all(grepl(".rds$", reads)))
                    stop("Reads should either be: a vector of paths to .bam files, ", 
                         "a vector of paths to Bambu RCfile .rds files, ",
                         "or a list of loaded Bambu RCfiles")
                # if bam files are loaded in check that a genome is provided
                if (all(grepl(".bam$", reads)) & is.null(genomeSequence)){
                    stop("A genome must be provided when running bambu from bam files")
                }
            }
        }
    }else if(is.null(quantData)){
        stop("Please provide either reads or quantData!",
             "Reads should either be: a vector of paths to .bam files, ", 
             "a vector of paths to Bambu RCfile .rds files, ",
             "or a list of loaded Bambu RCfiles. ",
             "quantData should be output from bambu with ",
             "assignDist = TRUE and quant = FALSE")
    }
    
    ## check genomeSequence can't be FaFile in Windows as faFile will be dealt
    ## strangely in windows system
    if (.Platform$OS.type == "windows") {
        if (is(genomeSequence, "FaFile")) 
            warning("Note that use of FaFile using Rsamtools in Windows is a bit
            fuzzy, recommend to provide the path as a string variable to avoid
            use of Rsamtools for opening.")
    }

    #check single-cell and spatial inputs match
    if(!is.null(sampleNames)){
        if(length(reads)!=length(sampleNames)){
            stop("There are not the same number of sampleNames as input files to reads. ",
            "Make sure these two arguments are vectors of the same length")
        }
    }

    if(!is.null(sampleData)){
        if (!all(grepl("\\.(csv|tsv|txt)$", na.omit(sampleData), ignore.case = TRUE))){
            stop("Not all paths for sample metadata files are .csv/.tsv/.txt files")
        }
        if(length(sampleData)==1 & length(reads)>1){ # one sample metadata for all samples
            message("Using the same sample metadata file for all input samples")
        } else if(length(reads)!=length(sampleData)){ # multiple sample metadatas for multiple samples
            stop(
                "The number of sample metadata files does not match the number of input read files. ",
                "These two arguments (sampleData & reads) must be vectors of the same length. ",
                "If a specific sample has no metadata, please use 'NA' as a placeholder in the sampleData vector."
            )
        }
    }
    return(annotations)
}


#' Function to create a object that can be queried by getSeq
#' Either from fa file, or BSGenome object
#' @importFrom methods is
#' @importFrom Rsamtools FaFile
#' @noRd
checkInputSequence <- function(genomeSequence) {
    if (is.null(genomeSequence)) stop("Reference genome sequence is missing,
        please provide fasta file or BSgenome name, see available.genomes()")
    if(is.character(genomeSequence)){
    if (genomeSequence %in% BSgenome::available.genomes()) {
        genomeSequence <- BSgenome::getBSgenome(genomeSequence)
        return(genomeSequence)
    } 
    tryCatch(
    {
        if (.Platform$OS.type == "windows") {
        genomeSequence <- Biostrings::readDNAStringSet(genomeSequence)
        newlevels <- unlist(lapply(strsplit(names(genomeSequence)," "),
                                    "[[", 1))
        names(genomeSequence) <- newlevels
        } else {
            indexFileExists <- file.exists(paste0(genomeSequence,".fai"))
            if (!indexFileExists) indexFa(genomeSequence)
            genomeSequence <- FaFile(genomeSequence)
        }
    },
    error=function(cond) {
        stop("Input genome file not readable.",
            " Requires a FASTA or BSgenome name")
    }
    )}
    return(genomeSequence)
}


#' Function that gathers warnings from several read class lists and outputs the counts
#' @noRd
handleWarnings <- function(readClassList, verbose){
    warnings <- list()
    sampleNames <- c()
    for(i in seq_along(readClassList)){
        readClassSe <- readClassList[[i]]
        if (is.character(readClassSe))
            readClassSe <- readRDS(file = readClassSe)
        
        warnings[[i]] <- metadata(readClassSe)$warnings
        
        if(is.null(metadata(readClassSe)$warnings)) 
            warnings[[i]] <- NA
        sampleNames <- c(sampleNames, colnames(readClassSe))
    }
    names(warnings) <- sampleNames
    if(verbose & any(lengths(warnings)>0)){
        message("--- per sample warnings during read class construction ---")
        for(i in seq_along(warnings)){
            if(lengths(warnings)[i]>0){
                message("Warnings for: ", sampleNames[i])
                sapply(warnings[[i]], message)
            }
        }
    } else {
        message("Detected ", sum(lengths(warnings)), " warnings across the samples during ",
        "read class construction. Access warnings with metadata(bambuOutput)$warnings")
    }
    return(warnings)
}

#' Calculate the dist table used for Bambu Quantification
calculateDistTable <- function(readClassList, annotations, isoreParameters, rcAssignmentParameters, verbose, returnDistTable){
    readClassDist <- isore.estimateDistanceToAnnotations(readClassList, annotations,
                                                            min.exonDistance = rcAssignmentParameters[["min.exonDistance"]],
                                                            min.primarySecondaryDist = rcAssignmentParameters[['min.primarySecondaryDist']],
                                                            min.primarySecondaryDistStartEnd = rcAssignmentParameters[['min.primarySecondaryDistStartEnd2']],
                                                            verbose = verbose)
        metadata(readClassDist)$distTable <- modifyIncompatibleAssignment(metadata(readClassDist)$distTable)
        if(returnDistTable) metadata(readClassDist)$distTableOld <- metadata(readClassDist)$distTable
                #convert string gene ids into index to save memory
        GENEIDs <- factor(unique(mcols(annotations)$GENEID))
        GENEID.i <- as.numeric(GENEIDs)
        metadata(readClassDist)$distTable$GENEID.i <- GENEID.i[match(metadata(readClassDist)$distTable$GENEID, GENEIDs)]
        metadata(readClassDist)$distTable.incompatible <- data.table(as.data.frame(metadata(readClassDist)$distTable)) %>% 
            filter(grepl("unidentified", annotationTxId)) %>% distinct(readClassId, .keep_all = TRUE)
        metadata(readClassDist)$distTable <- genEquiRCsBasedOnObservedReads(readClassDist)     
        return(readClassDist)
}

#' Combine combined count se object from multiple samples, cells or spatial locations
#' @noRd
combineCountSes <- function(countsSe, colDataList, annotations){
    countsData <- c("counts", "CPM", "fullLengthCounts", "uniqueCounts", "incompatibleCounts")
    sampleNames <- names(countsSe)
    countsDataMat <- lapply(countsData, FUN = function(k){
        countsVecList <- lapply(countsSe, function(j){j[[k]]})
        countsMat <- sparseMatrix(i = unlist(lapply(countsVecList, function(j) j@i)),
                                j = unlist(lapply(seq_along(countsVecList), function(j) rep(j, length(countsVecList[[j]]@i)))),
                                x = unlist(lapply(countsVecList, function(j) j@x)),
                                dims = c(length(countsVecList[[1]]), length(countsVecList)))
        if(all(is.na(countsMat)))
            countsMat <- sparseMatrix(i=NULL, j = NULL, dims = c(length(countsVecList[[1]]), length(countsVecList)))
        
        colnames(countsMat) <- sampleNames
        
        if (k == "incompatibleCounts")
            rownames(countsMat) <- unique(mcols(annotations)$GENEID)
        
        return(countsMat)
    })
    names(countsDataMat) <- countsData
    combinedCountsSe <- SummarizedExperiment(assays = SimpleList(counts = countsDataMat$counts, 
                                                        CPM = countsDataMat$CPM, 
                                                        fullLengthCounts = countsDataMat$fullLengthCounts, 
                                                        uniqueCounts = countsDataMat$uniqueCounts))
    metadata(combinedCountsSe)$incompatibleCounts <- countsDataMat$incompatibleCounts
    rowRanges(combinedCountsSe) <- annotations

    colData(combinedCountsSe) <- DataFrame(bind_rows(colDataList))
    
    return(combinedCountsSe)
}

#' Generate the colData using the external sampleMetadata.csv provided by the user in the sampleMetadata argument
#' @param readClassList A list object containingmetadata about read classes.
#' @param sampleMetadata A path to a CSV file or NULL/NA if there is no metadata for the sample.
#' @param demultiplexed Logical; indicates if data is demultiplexed.
#'
#' @return A DataFrame containing colData for the sample.
#' @export
generateColData <- function(readClassList, sampleMetadata, demultiplexed) {
  sampleMetadataDf <- if (is.null(sampleMetadata) || is.na(sampleMetadata)) {
    if (demultiplexed) tibble(barcode = character()) else tibble(sampleName = character())
  } else {
    fread(sampleMetadata)
  }

  joinKey <- if (demultiplexed) "barcode" else "sampleName"

  colData <- tibble(
      id = metadata(readClassList)$sampleData$id, 
      sampleName = metadata(readClassList)$sampleData$sampleName
  ) 

  if (demultiplexed) {
      colData <- colData %>%
        mutate(barcode = metadata(readClassList)$sampleData$barcode)
  }
  
  colData <- colData %>%
    left_join(sampleMetadataDf, by = joinKey) %>%
    as.data.frame()
  
  rownames(colData) <- colData$id
  
  colData
}

# Quick wrapper function (https://stackoverflow.com/questions/13273833/merging-multiple-data-tables)
#' @noRd 
merge_wrapper <- function(x,y){
    merge.data.table(x,y,by = "GENEID",all=TRUE)
}

