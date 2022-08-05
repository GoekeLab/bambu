
## Functions to set basic parameters and check inputs
#' setBiocParallelParameters
#' @importFrom BiocParallel bpparam
#' @noRd
setBiocParallelParameters <- function(reads, ncore, verbose){
    bpParameters <- bpparam()
    #===# set parallel options: otherwise use parallel to distribute samples
    bpParameters$workers <- ifelse(length(reads) == 1, 1, ncore)
    bpParameters$progressbar <- (!verbose)
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
        min.exonOverlap = 10, #
        min.primarySecondaryDist = 5,
        min.primarySecondaryDistStartEnd1 = 5, # for creating new annotations
        min.primarySecondaryDistStartEnd2 = 5, # for read assignment
        min.exonOverlap = 10,
        min.txScore.multiExon = 0,
        min.txScore.singleExon = 1,
        fitReadClassModel = TRUE,
        min.readFractionByEqClass = 0,
        prefix = "") 
    isoreParameters <- 
        updateParameters(isoreParameters, isoreParameters.default)
    return(isoreParameters)
}


#' setEmParameters
#' @noRd
setEmParameters <- function(emParameters){
    emParameters.default <- list(degradationBias = TRUE, maxiter = 10000, 
        conv = 10^(-2), minvalue = 10^(-8))
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
            if (!all(c("TXNAME", "GENEID", "eqClass") %in% 
                     colnames(mcols(annotations)))) 
                stop("The annotations is not properly prepared.\nPlease 
                    prepareAnnnotations using prepareAnnotations function.")
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
    } else {
        stop("Annotations is missing.")
    }
    # ===# Check whether provided readClass.outputDir exists  #===#
    if (!is.null(readClass.outputDir)) {
        if (!dir.exists(readClass.outputDir)) 
            stop("output folder does not exist")
    }
    # ===# Check whether provided read files are all in the same format (.bam or .rds) #===#
    if (!all(sapply(reads, class)=="RangedSummarizedExperiment") 
        & !all(grepl(".bam$", reads)) & !all(grepl(".rds$", reads)))
            stop("Reads should either be: a vector of paths to .bam files, ", 
            "a vector of paths to Bambu RCfile .rds files, ",
            "or a list of loaded Bambu RCfiles")
    # if bam files are loaded in check that a genome is provided
    if (all(grepl(".bam$", reads)) & is.null(genomeSequence)){
        stop("A genome must be provided when running bambu from bam files")
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

#' Function to load in a vector of read class files paths
#' @noRd
loadReadClassFiles <- function(rcFiles){
    #If vector of paths load in the files
    if (all(grepl(".rds", rcFiles))) 
        rcFiles = lapply(rcFiles, FUN = function(x){readRDS(x)})
    return(rcFiles)
}