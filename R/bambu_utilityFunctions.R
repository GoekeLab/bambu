
## Functions to set basic parameters and check inputs
#' setBiocParallelParameters
#' @importFrom BiocParallel bpparam
#' @noRd
setBiocParallelParameters <- function(reads, readClass.file, ncore, verbose){
    bpParameters <- bpparam()
    #===# set parallel options: otherwise use parallel to distribute samples
    bpParameters$workers <- ifelse(max(length(reads),
        length(readClass.file)) == 1, 1, ncore)
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
        max.txNDR = 0.1,
        fitReadClassModel = TRUE,
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
checkInputs <- function(annotations, reads, readClass.file,
                        readClass.outputDir, genomeSequence){
    # ===# Check annotation inputs #===#
    if (!is.null(annotations)) {
        if (is(annotations, "TxDb")) {
            annotations <- prepareAnnotations(annotations)
        } else if (is(annotations, "CompressedGRangesList")) {
            ## check if annotations is as expected
            if (!all(c("TXNAME", "GENEID", "eqClass") %in% 
                colnames(mcols(annotations)))) 
                stop("The annotations is not properly prepared.\nPlease 
                    prepareAnnnotations using prepareAnnotations function.")
        } else {
            stop("The annotations is not a GRangesList object.")
        }
    } else {
        stop("Annotations is missing.")
    }
    ## When SE object from bambu.quantISORE is provided ##
    if (!is.null(reads) & (!is.null(readClass.file))) stop("At least bam file or
        path to readClass file needs to be provided.")
    # ===# Check whether provided readClass.outputDir exists  #===#
    if (!is.null(readClass.outputDir)) {
        if (!dir.exists(readClass.outputDir)) 
            stop("output folder does not exist")
    }
    # ===# Check whether provided readclass files are all in rds format #===#
    if (!is.null(readClass.file)) {
        if (!all(grepl(".rds", readClass.file))) 
            stop("Read class files should be provided in rds format.")
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
    if (is(genomeSequence, "character")) {
        if (grepl(".fa", genomeSequence)) {
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
        } else {
            genomeSequence <- BSgenome::getBSgenome(genomeSequence)
        }
    }
    return(genomeSequence)
}
