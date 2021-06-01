#' Combine transcript candidates across samples
#' @param readClassList readClassList
#' @param stranded stranded
#' @param min.readCount minimum number of read counts per transcript
#' @param min.readFractionByGene minimum fraction of transcript usage
#' @param bpParameters biocParallel parameters
#' @param verbose verbose
#' @importFrom GenomicRanges GRanges
#' @importFrom SummarizedExperiment rbind
#' @importFrom dplyr as_tibble mutate_if select mutate %>% filter
#' @importFrom tidyr separate
#' @noRd
isore.combineTranscriptCandidates <- function(readClassList,
    stranded, ## stranded used for unspliced reduce  
    min.readCount = isoreParameters[["min.readCount"]], 
    min.readFractionByGene = isoreParameters[["min.readFractionByGene"]],
    bpParameters,
    verbose){
    combinedSplicedTranscripts <- 
        combineSplicedTranscriptModels(readClassList, bpParameters, 
        min.readCount, min.readFractionByGene, verbose)
    combinedUnsplicedTranscripts <- 
        combineUnsplicedTranscriptModels(readClassList,  bpParameters, 
        stranded, min.readCount, min.readFractionByGene, verbose)
    combinedTranscripts <- 
        bind_rows(combinedSplicedTranscripts, combinedUnsplicedTranscripts)
    return(combinedTranscripts)
}


#' combine spliced transcript models
#' @noRd
combineSplicedTranscriptModels <- function(readClassList, bpParameters, 
    min.readCount = isoreParameters[["min.readCount"]], 
    min.readFractionByGene = isoreParameters[["min.readFractionByGene"]],
    min.geneFDR = isoreParameters[["min.geneFDR"]],
    min.txFDR = isoreParameters[["min.txFDR"]],
    verbose){
    options(scipen = 999) #maintain numeric basepair locations not sci.notfi.
    start.time <- proc.time()
    featureTibbleList <- 
        bplapply(seq_along(readClassList), function(sample_id){
        extractFeaturesFromReadClassSE(readClassSe = readClassList[[sample_id]],
        sample_id = sample_id)}, BPPARAM = bpParameters)
    end.time <- proc.time()
    if (verbose) message("creating spliced feature tibble objects for all 
        samples in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
    ## update combinedFeatureTibble by sample, as at each step 
    ## start and end are updated iteratively as the fead count weighted average
    start.time <- proc.time()
    for ( s in seq_along(featureTibbleList)[-1]){
        combinedFeatureTibble <- 
        if_else(s==1, featureTibbleList[[s]],
        bind_rows(list(combinedFeatureTibble,featureTibbleList[[s]]))) %>% 
        group_by(intronStarts, intronEnds, chr, strand) %>%
        mutate(start = median(rep(start, times = readCount)),
            end = median(rep(end, times = readCount)), #weighted median is used 
        NSampleReadCount = sum(readCount >= min.readCount), 
        # number of samples passed read count criteria
        NSampleReadProp = sum(geneReadProp >= min.readFractionByGene),
        # number of samples passed gene read prop criteria (transcript usage)
        NSampleGeneFDR = sum(geneFDR <= min.geneFDR),
        NSampleTxFDR = sum(txFDR <= min.txFDR)) %>% 
        # number of samples passed gene and tx score criteria
        ungroup() 
        ## remember to ungroup to avoid unnecessary wrong selection later
    }
    combinedFeatureTibble <- combinedFeatureTibble %>% 
        separate(row_id, c("sample","rcName"), sep = "\\-") %>%
        select(-sample)
    end.time <- proc.time()
    if (verbose) message("combing spliced feature tibble objects across all
        samples in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
    return(combinedFeatureTibble)
}

#' extract important features from readClassSe object for each sample
#' @noRd
extractFeaturesFromReadClassSE <- function(readClassSe, sample_id){
    if (is.character(readClassSe)) 
        readClassSe <- readRDS(file = readClassSe)
    dimNames <- list(rownames(readClassSe), colnames(readClassSe))
    rowRangesSe <- rowRanges(readClassSe)
    start <- matrix(min(start(rowRangesSe)),
        dimnames = dimNames)
    end <- matrix(max(end(rowRangesSe)),
        dimnames = dimNames) 
    rowData <- as_tibble(rowData(readClassSe))
    rowData$start <- rowMins(start)
    rowData$end <- rowMaxs(end)
    rowData$row_id <- paste("s",sample_id,"-",dimNames[[1]],sep="")
    rowData$sample_name <- dimNames[[2]]
    featureTibble <- rowData %>% dplyr::select(chr = chr.rc, start, end,
        strand = strand.rc, intronStarts, intronEnds, confidenceType,
        readCount, geneReadProp, txScore, txFDR, geneScore, geneFDR,
        row_id, sample_name) %>%
        filter(readCount > 1, # only use readCount>1 and highconfidence reads
            confidenceType == "highConfidenceJunctionReads") %>%
        mutate(sample_id = sample_id)
    return(featureTibble)
}


#' combine unspliced transcript models
#' @importFrom tidyr separate
#' @importFrom dplyr %>% select mutate 
#' @importFrom summarizedExperiment rowRanges
#' @importFrom biocParallel bplapply
#' @noRd
combineUnsplicedTranscriptModels <- 
    function(readClassList,  bpParameters, unstranded,
    min.readCount = isoreParameters[["min.readCount"]], 
    min.readFractionByGene = isoreParameters[["min.readFractionByGene"]],
    verbose){
    start.time <- proc.time()
    newUnsplicedSeList <- bplapply(seq_along(readClassList), function(sample_id)
        extractNewUnsplicedRanges(readClassSe = readClassList[[sample_id]],
        sample_id = sample_id), BPPARAM = bpParameters)
    end.time <- proc.time()
    if (verbose) message("extract new unspliced ranges object for all
        samples in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
    rangesList <- bplapply(newUnsplicedSeList, function(newUnsplicedSe){
        rr <- unlist(rowRanges(newUnsplicedSe))
        rr$row_id <- names(rr)
        return(rr)
    }, BPPARAM = bpParameters)
    colDataNames <- unlist(lapply(newUnsplicedSeList, function(x) colnames(x)))
    start.time <- proc.time()
    combinedNewUnsplicedSe <- reduceUnsplicedRanges(rangesList)
    end.time <- proc.time()
    if (verbose) message("reduce new unspliced ranges object across all
        samples in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
    start.time <- proc.time()
    combinedUnsplicedTibble <- 
        makeUnsplicedTibble(combinedNewUnsplicedSe,newUnsplicedSeList,
        min.readCount, min.readFractionByGene, bpParameters)
    combinedUnsplicedTibble <- combinedUnsplicedTibble %>% 
        separate(row_id, c("sample","rcName"), sep = "\\-") %>%
        mutate(sample_id = as.integer(gsub("s","",sample))) %>%
        mutate(sample_name = colDataNames[sample_id]) %>%
        select(-sample)
    end.time <- proc.time()
    if (verbose) message("combine new unspliced tibble object across all
        samples in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
    return(combinedUnsplicedTibble)
}


#' extract new unspliced ranges from readClassSe object for each sample
#' @importFrom dplyr as_tibble
#' @importFrom summarizedExperiment rowData rownames
#' @noRd
extractNewUnsplicedRanges <- function(readClassSe, sample_id){
    if (is.character(readClassSe)) 
        readClassSe <- readRDS(file = readClassSe)
    rowData <- as_tibble(rowData(readClassSe))
    pre_names <- rownames(readClassSe)
    rownames(readClassSe) <- paste0("s",sample_id,"-",pre_names)
    newUnsplicedSe <- 
        readClassSe[which(rowData$confidenceType == "unsplicedNew" & 
        rowData$readCount > 1)]
    return(newUnsplicedSe)
}


#' reduce unspliced ranges
#' @importFrom dplyr as_tibble %>% mutate group_by summarise ungroup
#' @noRd
reduceUnsplicedRanges <- function(rangesList){
    combinedNewUnsplicedSe <- rangesList[[1]]
    for ( s in seq_along(rangesList)[-1]){
        new_range <- rangesList[[s]]
        new_range$exon_rank <- new_range$exon_endRank <- 
            names(new_range) <- NULL
        unlistedSe <- c(combinedNewUnsplicedSe, new_range)
        combinedNewUnsplicedSe <- 
            reduce(unlistedSe, with.revmap=TRUE,ignore.strand = !stranded)
        ## map it back to find the corresponding sample and id, update name
        rcNames <- as_tibble(as.data.frame(combinedNewUnsplicedSe$revmap)) %>%
            mutate(row_id = unlistedSe$row_id[value]) %>%
            group_by(group) %>%
            summarise(combinedName = paste(row_id, collapse = "+")) %>%
            ungroup()
        ## at each iteration update the new combined name of unspliced se 
        combinedNewUnsplicedSe$row_id <- rcNames$combinedName
        combinedNewUnsplicedSe$revmap <- NULL
    }
    return(combinedNewUnsplicedSe)
}

#' make unspliced tibble 
#' @importFrom tidyr separate_rows
#' @importFrom dplyr as_tibble rename mutate select %>% group_by left_join
#'              ungroup
makeUnsplicedTibble <- function(combinedNewUnsplicedSe,newUnsplicedSeList,
                                min.readCount, min.readFractionByGene,
                                bpParameters){
    newUnsplicedTibble <- as_tibble(combinedNewUnsplicedSe) %>%
        rename(chr = seqnames) %>%
        select(chr, start, end, strand, row_id) %>%
        separate_rows(row_id, sep = "\\+") 
    rowDataCombined <-
        do.call("rbind",bplapply(newUnsplicedSeList, function(newUnsplicedSe) {
        rr <- rowData(newUnsplicedSe[intersect(rownames(newUnsplicedSe), 
            newUnsplicedTibble$row_id)])
        rr <- as_tibble(rr) %>% select(intronStarts, intronEnds, confidenceType,
            readCount, geneReadProp, txScore, txFDR, geneScore, geneFDR) %>%
            mutate(row_id = rownames(rr))
        return(rr)
    } , BPPARAM = bpParameters))
    newUnsplicedTibble <- newUnsplicedTibble %>% 
        left_join(rowDataCombined, by =  "row_id") %>%
        group_by(chr, strand, start, end) %>% 
        mutate(NSampleReadCount = sum(sum(readCount) >= min.readCount), 
            NSampleReadProp = sum(sum(geneReadProp) >= 
            min.readFractionByGene)) %>%
        ungroup()
    return(newUnsplicedTibble)
}
