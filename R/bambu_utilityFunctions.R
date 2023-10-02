
## Functions to set basic parameters and check inputs
#' setBiocParallelParameters
#' @importFrom BiocParallel bpparam
#' @noRd
setBiocParallelParameters <- function(reads, ncore, verbose){
    if(ncore >= 2) message("WARNING - If you change the number of cores (ncore) ",
    "between Bambu runs and there is no progress please restart your R session ",
    "to resolve the issue that originates from the XGboost package.")
    bpParameters <- bpparam()
    #===# set parallel options: otherwise use parallel to distribute samples
    bpParameters$workers <- ifelse(length(reads) == 1, 1, ncore)
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
        min.exonDistance = 35,
        min.exonOverlap = 10, 
        min.primarySecondaryDist = 5,
        min.primarySecondaryDistStartEnd1 = 5, # for creating new annotations
        min.primarySecondaryDistStartEnd2 = 5, # for read assignment
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
checkInputs <- function(annotations, reads, readClass.outputDir, genomeSequence){
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
        if(any(grepl("^BambuGene", names(annotations))) | 
            any(grepl("^BambuTx", mcols(annotations)$TXNAME))){
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

    if (is(reads, "BamFileList")){
        if(is.null(genomeSequence)){
            stop("A genome must be provided when running bambu from bam files")
        }
    } else{
    # ===# Check whether provided read files are all in the same format (.bam or .rds) #===#
        isRDSs = all(sapply(reads, class)=="RangedSummarizedExperiment")
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
    ## check genomeSequence can't be FaFile in Windows as faFile will be dealt
    ## strangely in windows system
    if (.Platform$OS.type == "windows") {
        if (is(genomeSequence, "FaFile")) 
            warning("Note that use of FaFile using Rsamtools in Windows is a bit
            fuzzy, recommend to provide the path as a string variable to avoid
            use of Rsamtools for opening.")
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
            "Requires a FASTA or BSgenome name")
    }
    )}
    return(genomeSequence)
}


#' Function that gathers warnings from several read class lists and outputs the counts
#' @noRd
handleWarnings <- function(readClassList, verbose){
    warnings = list()
    sampleNames = c()
    for(i in seq_along(readClassList)){
        readClassSe = readClassList[[i]]
        if (is.character(readClassSe)) 
            readClassSe <- readRDS(file = readClassSe)
        warnings[[i]] = metadata(readClassSe)$warnings
        sampleNames = c(sampleNames, colnames(readClassList[[i]]))
    }
    names(warnings) = sampleNames

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


#' Combine count se object while preserving the metadata objects
#' @noRd
combineCountSes <- function(countsSe, trackReads = FALSE, returnDistTable = FALSE){
    sampleNames = sapply(countsSe, FUN = function(x){colnames(x)})
    if(trackReads){
        readToTranscriptMaps = lapply(countsSe, FUN = function(se){metadata(se)$readToTranscriptMap})
        names(readToTranscriptMaps) = sampleNames
        countsSe = lapply(countsSe, FUN = function(se){
            metadata(se)$readToTranscriptMap=NULL
            return(se)})
    }
    if(returnDistTable){
        distTables = lapply(countsSe, FUN = function(se){metadata(se)$distTable})
        names(distTables) = sampleNames
        countsSe = lapply(countsSe, FUN = function(se){
            metadata(se)$distTable=NULL
            return(se)})
    }
    # combine incompatible counts
    incompatibleCounts = Reduce(merge_wrapper, lapply(countsSe, FUN = function(se){metadata(se)$incompatibleCounts}))
    countsSe = lapply(countsSe, FUN = function(se){
        metadata(se)$incompatibleCounts=NULL
        return(se)})
    countsSe <- do.call(SummarizedExperiment::cbind, countsSe)
    if(trackReads) metadata(countsSe)$readToTranscriptMaps = readToTranscriptMaps
    if(returnDistTable) metadata(countsSe)$distTables = distTables
    metadata(countsSe)$incompatibleCounts = incompatibleCounts
    return(countsSe)
}

# Quick wrapper function (https://stackoverflow.com/questions/13273833/merging-multiple-data-tables)
#' @noRd 
merge_wrapper <- function(x,y){
    merge.data.table(x,y,by = "GENEID",all=TRUE)
}

