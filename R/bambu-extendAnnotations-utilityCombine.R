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
    min.readCount , min.readFractionByGene , min.geneFDR,
    min.txFDR, bpParameters,verbose){
    combinedSplicedTranscripts <- 
        combineSplicedTranscriptModels(readClassList, bpParameters, 
        min.readCount, min.readFractionByGene, min.geneFDR, min.txFDR, verbose)
    combinedSplicedTranscripts$confidenceType <-  "highConfidenceJunctionReads"
    combinedUnsplicedTranscripts <- 
        combineUnsplicedTranscriptModels(readClassList, bpParameters, 
        stranded, min.readCount, min.readFractionByGene, min.geneFDR, 
        min.txFDR, verbose)
    combinedUnsplicedTranscripts$confidenceType <- "unsplicedNew"
    combinedTranscripts <- 
        bind_rows(combinedSplicedTranscripts, combinedUnsplicedTranscripts)
    return(combinedTranscripts)
}


#' combine spliced transcript models
#' @noRd
combineSplicedTranscriptModels <- function(readClassList, bpParameters, 
    min.readCount, min.readFractionByGene , min.geneFDR, min.txFDR,
    verbose){
    options(scipen = 999) #maintain numeric basepair locations not sci.notfi.
    start.ptm <- proc.time()
    featureTibbleList <- 
        bplapply(seq_along(readClassList), function(sample_id){
        extractFeaturesFromReadClassSE(readClassSe = readClassList[[sample_id]],
        sample_id = sample_id)}, BPPARAM = bpParameters)
    end.ptm <- proc.time()
    if (verbose) message("creating spliced feature tibble objects for all 
        samples in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
    ## update combinedFeatureTibble by sample, as at each step 
    ## start and end are updated iteratively as the fead count weighted average
    start.ptm <- proc.time()
    combinedFeatureTibble <- NULL
    for (s in seq_along(featureTibbleList)){
        combinedFeatureTibble <- combineFeatureTibble(combinedFeatureTibble,
            featureTibbleList[[s]])
    }
    end.ptm <- proc.time()
    if (verbose) message("combing spliced feature tibble objects across all
        samples in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
    return(combinedFeatureTibble)
}

#' Function to combine featureTibble and create the NSample variables 
#' @noRd
combineFeatureTibble <- function(combinedFeatureTibble, newSampleTibble){
    # step 1: summarise within each sample
    sampleName <- unique(newSampleTibble$sample_name)
    group_var <- c("intronStarts", "intronEnds", "chr", "strand")
    sum_var <- c("readCount","start","end","NSampleReadCount",
                 "NSampleReadProp","NSampleGeneFDR","NSampleTxFDR",
                 sampleName)
    featureTibbleSummarised <- newSampleTibble %>% 
        mutate(NSampleReadCount = (readCount >= min.readCount), 
            # number of samples passed read count criteria
            NSampleReadProp = (geneReadProp >= min.readFractionByGene),
            # number of samples passed gene read prop criteria
            NSampleGeneFDR = (geneFDR <= min.geneFDR),
            NSampleTxFDR = (txFDR <= min.txFDR),
            {{sampleName}} := readCount) %>%
        select(all_of(c(group_var, sum_var)))
    # step2: combine with previous results 
    combinedFeatureTibble <- bind_rows(list(combinedFeatureTibble,
                                           featureTibbleSummarised))
    countVars <- setdiff(colnames(combinedFeatureTibble), c(group_var,
            "start","end"))
    combinedFeatureTibble <- combinedFeatureTibble %>% 
        mutate(readCount_tmp = readCount) %>%
        group_by(intronStarts, intronEnds, chr, strand) %>%
        mutate(start = median(rep(start, times = readCount_tmp)),
            end = median(rep(end, times = readCount_tmp))) %>% 
        group_by(intronStarts, intronEnds, chr, strand, start, end) %>%
        summarise_at(countVars, sum, na.rm = TRUE) %>%
        ungroup()
    ## remember to ungroup to avoid unnecessary wrong selection later
    return(combinedFeatureTibble)
}



#' Left join with each sample to get the readCount

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
    rowData$sample_name <- dimNames[[2]]
    featureTibble <- rowData %>% dplyr::select(chr = chr.rc, start, end,
        strand = strand.rc, intronStarts, intronEnds, confidenceType,
        readCount, geneReadProp, txFDR,geneFDR,
        sample_name) %>%
        filter(readCount > 1, # only use readCount>1 and highconfidence reads
            confidenceType == "highConfidenceJunctionReads") 
    return(featureTibble)
}


#' combine unspliced transcript models
#' @importFrom tidyr separate
#' @importFrom dplyr %>% select mutate 
#' @importFrom summarizedExperiment rowRanges
#' @importFrom biocParallel bplapply
#' @noRd
combineUnsplicedTranscriptModels <- 
    function(readClassList,  bpParameters, stranded, min.readCount, 
             min.readFractionByGene, min.geneFDR, min.txFDR,
             verbose){
        start.ptm <- proc.time()
        newUnsplicedSeList <- 
            bplapply(seq_along(readClassList), function(sample_id)
            extractNewUnsplicedRanges(readClassSe = readClassList[[sample_id]],
            sample_id = sample_id), BPPARAM = bpParameters)
        end.ptm <- proc.time()
        if (verbose) message("extract new unspliced ranges object for all
        samples in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
        rangesList <- bplapply(newUnsplicedSeList, function(newUnsplicedSe){
            rr <- unlist(rowRanges(newUnsplicedSe))
            rr$row_id <- names(rr)
            return(rr)
        }, BPPARAM = bpParameters)
        colDataNames <- unlist(lapply(newUnsplicedSeList, function(x) colnames(x)))
        start.ptm <- proc.time()
        combinedNewUnsplicedSe <- reduceUnsplicedRanges(rangesList, stranded)
        end.ptm <- proc.time()
        if (verbose) message("reduce new unspliced ranges object across all
        samples in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
        start.ptm <- proc.time()
        combinedUnsplicedTibble <- 
            makeUnsplicedTibble(combinedNewUnsplicedSe,newUnsplicedSeList, 
            colDataNames, min.readCount, min.readFractionByGene, min.geneFDR,
            min.txFDR, bpParameters)
        end.ptm <- proc.time()
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
reduceUnsplicedRanges <- function(rangesList, stranded){
    unlistedSe <- do.call("c",rangesList)
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
    return(combinedNewUnsplicedSe)
}

#' make unspliced tibble 
#' @importFrom tidyr separate_rows pivot_wider
#' @importFrom dplyr as_tibble rename mutate select %>% group_by left_join
#'              ungroup
makeUnsplicedTibble <- function(combinedNewUnsplicedSe,newUnsplicedSeList,
        colDataNames,min.readCount, min.readFractionByGene,
        min.geneFDR, min.txFDR,bpParameters){
    newUnsplicedTibble <- as_tibble(combinedNewUnsplicedSe) %>%
        rename(chr = seqnames) %>%
        select(chr, start, end, strand, row_id) %>%
        separate_rows(row_id, sep = "\\+") 
    rowDataCombined <-
        do.call("rbind",bplapply(newUnsplicedSeList, function(newUnsplicedSe) {
            rr <- rowData(newUnsplicedSe[intersect(rownames(newUnsplicedSe), 
                                                   newUnsplicedTibble$row_id)])
            rr <- as_tibble(rr) %>% select(confidenceType,
                readCount, geneReadProp, txScore, txFDR, geneScore, geneFDR) %>%
                mutate(row_id = rownames(rr))
            return(rr)
        } , BPPARAM = bpParameters))
    newUnsplicedTibble <- newUnsplicedTibble %>% 
        left_join(rowDataCombined, by =  "row_id") %>%
        separate(row_id, c("sample","rcName"), sep = "\\-") %>%
        mutate(sample_id = as.integer(gsub("s","",sample))) %>%
        mutate(sample_name = colDataNames[sample_id]) %>%
        select(-sample, -sample_id) %>%
        mutate(readCount_tmp = readCount) %>%
        group_by(chr,strand, start, end, sample_name) %>%
        summarise(readCount = sum(readCount),
                  geneReadProp = sum(geneReadProp),
                  geneFDR = median(geneFDR, times = readCount_tmp),
                  txFDR = median(txFDR, times = readCount_tmp)) %>%
        group_by(chr, strand, start, end) %>% 
        mutate(NSampleReadCount = sum(readCount >= min.readCount), 
               NSampleReadProp = sum(geneReadProp >= 
                                         min.readFractionByGene),
               NSampleGeneFDR = sum(geneFDR<= 
                                        min.geneFDR),
               NSampleTxFDR = sum(txFDR <= 
                                      min.txFDR)) %>%
        select(-geneReadProp, -geneFDR, -txFDR)  %>%
        pivot_wider(names_from = sample_name, values_from = readCount)
    return(newUnsplicedTibble)
}
