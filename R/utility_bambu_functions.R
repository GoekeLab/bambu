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
        min.primarySecondaryDist = 
            isoreParameters[['min.primarySecondaryDist']], 
        min.primarySecondaryDistStartEnd = 
            isoreParameters[['min.primarySecondaryDistStartEnd1']],
        prefix = isoreParameters[["prefix"]],
        verbose = verbose
    )
    return(annotations)
}

#' Perform quantification
#' @inheritParams bambu
#' @noRd
bambu.quantify <- function(readClass, annotations, emParameters,ncore = 1,
    verbose = FALSE, min.exonDistance = 35, min.primarySecondaryDist = 5, 
    min.primarySecondaryDistStartEnd = 5) {
    if (is.character(readClass)) readClass <- readRDS(file = readClass)
    readClass <- isore.estimateDistanceToAnnotations(readClass, annotations,
        min.exonDistance = min.exonDistance,
        min.primarySecondaryDist = min.primarySecondaryDist,
        min.primarySecondaryDistStartEnd = min.primarySecondaryDistStartEnd,
        verbose = verbose)
    txLength <- data.table(annotationTxId = names(annotations),
        txLength = sum(width(annotations)))
    ## For unique part of rc, we use first exon as the proxy
    readClassDt <- getEmptyClassFromSE(readClass, annotations)
    readClassDt <- 
        modifyReadClassWtFullLengthTranscript(readClassDt, annotations)
    readClassDt[, tx_len := txLength[match(gsub("Start","",tx_id),
        annotationTxId)]$txLength]
    d_rate <- 
        calculateExpectedCoverageRatio(readClass, annotations, txLength)
    colNameRC <- colnames(readClass)
    colDataRC <- colData(readClass)
    counts <- bambu.quantDT(readClassDt, emParameters = emParameters,
        ncore = ncore, verbose = verbose, d_rate)
    counts <- bambu.formatOutput(counts, annotations)
    seOutput <- SummarizedExperiment::SummarizedExperiment(
        assays = SimpleList(counts = matrix(counts$estimates, ncol = 1,
            dimnames = list(NULL, colNameRC)), CPM = matrix(counts$CPM,
            ncol =  1, dimnames = list(NULL, colNameRC)),
            fullLengthCounts = matrix(counts$FullLengthCounts, ncol = 1,
            dimnames = list(NULL, colNameRC)),
            partialLengthCounts = matrix(counts$PartialLengthCounts, 
            ncol = 1, dimnames = list(NULL, colNameRC))), colData = colDataRC)
    return(seOutput)
}

#' @noRd
bambu.formatOutput <- function(counts,annotations){
    Est <- c("FullLength", "PartialLength")
    counts[, estimate_type := ifelse(grepl("Start",tx_name),Est[1],Est[2])]
    counts[, tx_name := gsub("Start","",tx_name)]
    countsEstimates <- dcast(counts, tx_name ~ estimate_type,
        fun.aggregate = max, fill = 0, na.rm = TRUE, value.var = "estimates")
    setnames(countsEstimates, Est, paste0(Est,"Counts"))
    countsCPM <- dcast(counts, tx_name ~ estimate_type, 
        fun.aggregate = max,fill = 0, na.rm = TRUE, value.var = "CPM")
    setnames(countsCPM, Est, paste0(Est,"CPM"))
    counts <- merge(countsEstimates, countsCPM, by = "tx_name")
    counts[is.na(FullLengthCounts), `:=`(FullLengthCounts = 0,
        FullLengthCPM = 0) ]
    counts[is.na(PartialLengthCounts),`:=`(PartialLengthCounts = 0,
        PartialLengthCPM = 0) ]
    counts[, `:=`(estimates = FullLengthCounts + PartialLengthCounts,
                CPM = FullLengthCPM + PartialLengthCPM)]
    if (length(setdiff(counts$tx_name, names(annotations)))) 
        stop("The provided annotation is incomplete")
    counts <- counts[data.table(tx_name = names(annotations)), on = "tx_name"]
    counts[is.na(estimates),`:=`(estimates = 0, CPM = 0, FullLengthCounts = 0,
        PartialLengthCounts = 0, FullLengthCPM = 0, PartialLengthCPM = 0)]
    return(counts)
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
    ncore = 1, verbose = FALSE, d_rate = NULL) {
    if (is.null(readClassDt)) {
        stop("Input object is missing.")
    } else if (any(!(c("gene_id", "tx_id", "read_class_id","nobs") %in% 
        colnames(readClassDt)))) {
        stop("Columns gene_id, tx_id, read_class_id, nobs,
            are missing from object.")
    }
    ## ----step2: match to simple numbers to increase claculation efficiency
    geneVec <- unique(readClassDt$gene_id)
    ori_txvec <- unique(gsub("Start","",readClassDt$tx_id))
    txVec <- unique(readClassDt$tx_id)
    readclassVec <- unique(readClassDt$read_class_id)
    readClassDt <- as.data.table(readClassDt)
    readClassDt[, gene_sid := match(gene_id, geneVec)]
    readClassDt[, tx_sid := match(tx_id, txVec)]
    readClassDt[, tx_ori := match(gsub("Start","",tx_id),ori_txvec)]
    readClassDt[, read_class_sid := match(read_class_id, readclassVec)]
    readClassDt[, fullTx := grepl("Start",tx_id)]
    readClassDt[, `:=`(tx_id = NULL, gene_id = NULL, read_class_id = NULL)]

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