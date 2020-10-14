## Internal functions for bambu =========================
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
        }
        combinedTxCandidates <- isore.combineTranscriptCandidates(readClass,
            readClassSeRef = combinedTxCandidates, verbose = verbose
        )
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
    return(annotations)
}

#' Perform quantification
#' @inheritParams bambu
#' @noRd
bambu.quantify <- function(readClass, annotations, emParameters,
    min.exonDistance = 35, ncore = 1, verbose = FALSE) {
    if (is.character(readClass)) {
        readClass <- readRDS(file = readClass)
    }
    readClass <- isore.estimateDistanceToAnnotations(readClass, annotations,
        min.exonDistance = min.exonDistance, verbose = verbose)
    readClassDt <- getEmptyClassFromSE(readClass, annotations)
    colNameRC <- colnames(readClass)
    colDataRC <- colData(readClass)
    counts <- bambu.quantDT(readClassDt,
        emParameters = emParameters,
        ncore = ncore, verbose = verbose)
    if (length(setdiff(counts$tx_name, names(annotations)))) 
        stop("The provided annotation is incomplete")
    counts <- counts[data.table(tx_name = names(annotations)), on = "tx_name"]
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