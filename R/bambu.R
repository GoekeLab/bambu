#' Main function
#' @title long read isoform reconstruction and quantification
#' @description This function takes bam file of genomic alignments and performs
#' isoform recontruction and gene and transcript expression quantification.
#' It also allows saving of read class files of alignments, extending provided
#' annotations, and quantification based on extended annotations. When multiple
#' samples are provided, extended annotations will be combined across samples to
#' allow comparison.
#' @param reads A string or a vector of strings specifying the paths of bam
#' files for genomic alignments, or a \code{\link{BamFile}} object or a
#' \code{\link{BamFileList}}  object (see \code{\link{Rsamtools}}).
#' @param readClass.file A string or a vector of strings specifying the read
#' class files that are saved during previous run of \code{\link{bambu}}.
#' @param readClass.outputDir A string variable specifying the path to where
#' read class files will be saved.
#' @param annotations A \code{\link{TxDb}} object or A GRangesList object
#' obtained by \code{\link{prepareAnnotations}}.
#' @param genomeSequence A fasta file or a BSGenome object.
#' @param stranded A boolean for strandedness, defaults to FALSE.
#' @param ncore specifying number of cores used when parallel processing 
#' is used, defaults to 1.
#' @param yieldSize see \code{\link{Rsamtools}}.
#' @param isoreParameters A list of controlling parameters for isoform
#' reconstruction process:
#' \itemize{
#'     \item prefix specifying prefix for new gene Ids (genePrefix.number),
#'     defaults to empty
#'     \item remove.subsetTx indicating whether filter to remove read classes
#'     which are a subset of known transcripts(), defaults to TRUE
#'     \item min.readCount specifying minimun read count to consider a read
#'     class valid in a sample, defaults to 2
#'     \item min.readFractionByGene specifying minimum relative read count per
#'     gene, highly expressed genes will have many high read count low relative
#'     abundance transcripts that can be filtered, defaults to 0.05
#'     \item min.sampleNumber specifying minimum sample number with minimum read
#'     count, defaults to 1
#'     \item min.exonDistance specifying minum distance to known transcript 
#'     to be considered valid as new, defaults to 35
#'     \item min.exonOverlap specifying minimum number of bases shared with
#'     annotation to be assigned to the same gene id, defaults 10 base pairs
#' }
#' @param emParameters A list of controlling parameters for quantification
#' algorithm estimation process:
#' \itemize{
#'     \item maxiter specifying maximum number of run interations,
#'     defaults to 10000.
#'     \item bias specifying whether to correct for bias, defaults to FALSE.
#'     \item conv specifying the covergence trheshold control,
#'     defaults to 0.0001.
#' }
#' @param extendAnnotations A logical variable indicating whether annotations
#' are to be extended for quantification.
#' @param verbose A logical variable indicating whether processing messages will
#' be printed.
#' @details
#' @return A list of two SummarizedExperiment object for transcript expression
#' and gene expression.
#' @examples
#'
#' ## =====================
#' test.bam <- system.file("extdata",
#'     "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam",
#'     package = "bambu"
#' )
#' gr <- readRDS(system.file("extdata",
#'     "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds",
#'     package = "bambu"
#' ))
#' fa.file <- system.file("extdata",
#'     "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa",
#'     package = "bambu"
#' )
#' se <- bambu(
#'     reads = test.bam, annotations = gr,
#'     genomeSequence = fa.file, extendAnnotations = FALSE
#' )
#' @export
bambu <- function(reads = NULL, readClass.file = NULL,
    readClass.outputDir = NULL, annotations = NULL, genomeSequence = NULL,
    stranded = FALSE, ncore = 1, yieldSize = NULL, isoreParameters = NULL,
    emParameters = NULL, extendAnnotations = TRUE, verbose = FALSE) {
    annotations <-
        checkInputs(annotations, reads, readClass.file, readClass.outputDir)
    # ===# set default controlling parameters for isoform reconstruction  #===#
    isoreParameters.default <- list(
        remove.subsetTx = TRUE, min.readCount = 2,
        min.readFractionByGene = 0.05, min.sampleNumber = 1,
        min.exonDistance = 35, min.exonOverlap = 10, prefix = "") 
    isoreParameters <- checkParameters(isoreParameters, isoreParameters.default)
    emParameters.default <- list(bias = TRUE, maxiter = 10000, conv = 10^(-4))
    emParameters <- checkParameters(emParameters, emParameters.default)
    rm.readClassSe <- FALSE # indicator to remove temporary read class files
    bpParameters <- BiocParallel::bpparam()
    #===# set parallel options: otherwise use parallel to distribute samples
    bpParameters$workers <- ifelse(length(reads) == 1, 1, ncore)
    bpParameters$progressbar <- (!verbose)
    if (bpParameters$workers > 1) ncore <- 1
    readClassList <- processReads( reads, annotations, genomeSequence,
        readClass.outputDir,yieldSize, bpParameters, stranded,
        ncore, verbose)
    if (extendAnnotations) {
        annotations <- bambu.extendAnnotations(readClassList, annotations,
            isoreParameters, verbose = verbose)
        if (!verbose) message("Finished extending annotations.")
        gc(verbose = FALSE)
    }
    if (!verbose) message("Start isoform quantification")
    countsSe <- BiocParallel::bplapply(readClassList,
        bambu.quantify,annotations = annotations,
        min.exonDistance = isoreParameters[["min.exonDistance"]],
        emParameters = emParameters, ncore = ncore,
        verbose = verbose, BPPARAM = bpParameters
    )
    countsSe <- do.call(SummarizedExperiment::cbind, countsSe)
    rowRanges(countsSe) <- annotations
    if (!verbose) message("Finished isoform quantification.")
    # ===# Clean up temp directory
    if (rm.readClassSe) {
        file.remove(unlist(readClassList))
    }
    return(countsSe)
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
processReads <- function(reads, annotations, genomeSequence,
    readClass.outputDir, yieldSize, bpParameters, stranded, ncore, verbose) {
    if (!is.null(reads)) { # calculate readClass objects
        # ===# create BamFileList object from character #===#
        if (is(reads, "BamFile")) {
            if (!is.null(yieldSize)) {
                Rsamtools::yieldSize(reads) <- yieldSize
            } else {
                yieldSize <- Rsamtools::yieldSize(reads)
            }
        reads <- Rsamtools::BamFileList(reads)
        names(reads) <- tools::file_path_sans_ext(BiocGenerics::basename(reads))
        } else if (is(reads, "BamFileList")) {
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
        #===# When more than 10 samples, files saved to temporary directory
        if (length(reads) > 10 & (is.null(readClass.outputDir))) {
            readClass.outputDir <- tempdir()
            message(paste0("There are more than 10 samples, read class files
                will be temporarily saved to ", readClass.outputDir,
                " for more efficient processing"))
            rm.readClassSe <- TRUE # remove temporary read class files 
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
    } else {
        readClassList <- reads
    }
    return(readClassList)
}

#' check parameters for isore and em
#' @param Parameters parameters inputted by user
#' @param Parameters.default default parameters
#' @noRd
checkParameters <- function(Parameters, Parameters.default) {
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
#' @noRd
checkInputs <- function(annotations, reads, readClass.file,
    readClass.outputDir) {
    # ===# Check annotation inputs #===#
    if (!is.null(annotations)) {
        if (is(annotations, "TxDb")) {
            annotations <- prepareAnnotations(annotations)
        } else if (is(annotations, "CompressedGRangesList")) {
            ## check if annotations is as expected
            if (!all(c("TXNAME", "GENEID", "eqClass") %in% 
                colnames(mcols(annotations)))) 
                stop("The annotations is not properly prepared.\nPlease 
                prepareAnnnotations using prepareAnnotations or 
                prepareAnnotationsFromGTF functions.")
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
    return(annotations)
}

#' Extend annotations
#' @inheritParams bambu
#' @noRd
bambu.extendAnnotations <- function(readClassList, annotations,
                                isoreParameters, verbose = FALSE) {
    combinedTxCandidates <- NULL
    start.ptm <- proc.time()
    for (readClassIndex in seq_along(readClassList)) {
        readClass <- readClassList[[readClassIndex]]
        if (is.character(readClass)) {
            readClass <- readRDS(file = readClass)
            seqlevelsStyle(readClass) <- seqlevelsStyle(annotations)[1]
        }
        combinedTxCandidates <- isore.combineTranscriptCandidates(readClass,
            readClassSeRef = combinedTxCandidates, verbose = verbose
        )
        rm(readClass)
        gc(verbose = FALSE)
    }
    end.ptm <- proc.time()
    if (verbose) message("combining transcripts in ",
        round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
    annotations <- isore.extendAnnotations(
        se = combinedTxCandidates,
        annotationGrangesList = annotations,
        remove.subsetTx = isoreParameters[["remove.subsetTx"]],
        min.readCount = isoreParameters[["min.readCount"]],
        min.readFractionByGene = isoreParameters[["min.readFractionByGene"]],
        min.sampleNumber = isoreParameters[["min.sampleNumber"]],
        min.exonDistance = isoreParameters[["min.exonDistance"]],
        min.exonOverlap = isoreParameters[["min.exonOverlap"]],
        prefix = isoreParameters[["prefix"]],
        verbose = verbose
    )
    rm(combinedTxCandidates)
    gc(verbose = FALSE)
    return(annotations)
}

#' Perform quantification
#' @inheritParams bambu
#' @noRd
bambu.quantify <- function(readClass, annotations, emParameters,
    min.exonDistance = 35, ncore = 1, verbose = FALSE) {
    if (is.character(readClass)) {
        readClass <- readRDS(file = readClass)
        seqlevelsStyle(readClass) <- seqlevelsStyle(annotations)[1]
    }

    readClass <- isore.estimateDistanceToAnnotations(readClass, annotations,
        min.exonDistance = min.exonDistance, verbose = verbose)
    gc(verbose = FALSE)
    readClassDt <- getEmptyClassFromSE(readClass, annotations)
    colNameRC <- colnames(readClass)
    colDataRC <- colData(readClass)
    rm(readClass)
    gc(verbose = FALSE)
    counts <- bambu.quantDT(readClassDt,
        emParameters = emParameters,
        ncore = ncore, verbose = verbose)
    rm(readClassDt)
    gc(verbose = FALSE)
    if (length(setdiff(counts$tx_name, names(annotations))) > 0) 
        stop("The provided annotation is incomplete")
    counts <- counts[data.table(tx_name = names(annotations)), on = "tx_name"]
    rm(annotations)
    gc(verbose = FALSE)
    counts[is.na(estimates), `:=`(estimates = 0, CPM = 0)]

    seOutput <- SummarizedExperiment::SummarizedExperiment(
        assays =
            SimpleList(counts = matrix(counts$estimates, ncol = 1,
                dimnames = list(NULL, colNameRC)),CPM = matrix(counts$CPM,
                ncol = 1, dimnames = list(NULL, colNameRC))),
        colData = colDataRC)
    return(seOutput)
}

#' Preprocess bam files and save read class files
#' @inheritParams bambu
#' @noRd
bambu.constructReadClass <- function(bam.file, genomeSequence, annotations,
    readClass.outputDir = NULL, stranded = FALSE, ncore = 1, verbose = FALSE) {
    readGrgList <- prepareDataFromBam(bam.file[[1]], ncore = ncore,
        verbose = verbose)
    
    if (length(intersect(seqlevels(readGrgList),seqlevels(annotations))) == 0)
        stop("Error: please provide annotation with matched seqlevel styles.")
    
    se <- isore.constructReadClasses(
        readGrgList = readGrgList,
        runName = names(bam.file)[1],
        annotationGrangesList = annotations,
        genomeSequence = genomeSequence,
        stranded = stranded,
        ncore = ncore,
        verbose = verbose)
    seqlevels(se) <- unique(c(seqlevels(se), seqlevels(annotations)))
    if (!is.null(readClass.outputDir)) {
        readClassFile <- fs::path(readClass.outputDir, paste0(
            names(bam.file),"_readClassSe"), ext = "rds")
        if (file.exists(readClassFile)) {
            show(paste(readClassFile, "exists, will be overwritten"))
            # warning is not printed, use show in addition
            warning(paste(readClassFile, "exists, will be overwritten"))
        }
        saveRDS(se, file = readClassFile)
        se <- readClassFile
    }
    return(se)
}



#' Process data.table object
#' @param readClassDt A data.table object
#' @inheritParams bambu
#' @noRd
bambu.quantDT <- function(readClassDt = readClassDt, emParameters = NULL,
    ncore = 1, verbose = FALSE) {
    if (is.null(readClassDt)) {
        stop("Input object is missing.")
    } else if (any(!(c("gene_id", "tx_id", "read_class_id","nobs") %in% 
        colnames(readClassDt)))) {
        stop("Columns gene_id, tx_id, read_class_id, nobs,
            are missing from object.")
    }
    ## check quantification parameters
    emParameters.default <- list(bias = FALSE,maxiter = 10000,conv = 10^(-4))
    if (!is.null(emParameters)) {
        for (i in names(emParameters)) {
            emParameters.default[[i]] <- emParameters[[i]]
        }
    }
    emParameters <- emParameters.default
    ## ----step2: match to simple numbers to increase claculation efficiency
    geneVec <- unique(readClassDt$gene_id)
    txVec <- unique(readClassDt$tx_id)
    readclassVec <- unique(readClassDt$read_class_id)
    readClassDt <- as.data.table(readClassDt)
    readClassDt[, gene_sid := match(gene_id, geneVec)]
    readClassDt[, tx_sid := match(tx_id, txVec)]
    readClassDt[, read_class_sid := match(read_class_id, readclassVec)]
    readClassDt[, `:=`(tx_id = NULL, gene_id = NULL, read_class_id = NULL)]

    ## ----step3: aggregate read class
    temp <- aggReadClass(readClassDt)
    readClassDt <- temp[[1]]
    eqClassVec <- temp[[2]]

    ## ----step4: quantification
    start.time <- proc.time()
    outList <- abundance_quantification(readClassDt,ncore = ncore,
        bias = emParameters[["bias"]], maxiter = emParameters[["maxiter"]],
        conv = emParameters[["conv"]])
    end.time <- proc.time()
    if (verbose) message("Finished EM estimation in ",
        round((end.time - start.time)[3] / 60, 1), " mins.")

    theta_est <- outList[[1]]
    theta_est[, `:=`(tx_name = txVec[as.numeric(tx_sid)],
        gene_name = geneVec[gene_sid])]
    theta_est[, `:=`(tx_sid = NULL, gene_sid = NULL)]
    theta_est <- theta_est[, .(tx_name, estimates)]
    theta_est[, `:=`(CPM = estimates / sum(estimates) * (10^6))]
    return(theta_est)
}
