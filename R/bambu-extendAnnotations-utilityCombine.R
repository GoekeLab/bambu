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
    min.readCount , min.readFractionByGene,
    min.txScore.multiExon, min.txScore.singleExon, bpParameters ,verbose){
    combinedSplicedTranscripts <- 
        combineSplicedTranscriptModels(readClassList, bpParameters, 
        min.readCount, min.readFractionByGene, 
        min.txScore.multiExon, min.txScore.singleExon, verbose)
    combinedSplicedTranscripts[,confidenceType := "highConfidenceJunctionReads"]
    combinedUnsplicedTranscripts <- 
        combineUnsplicedTranscriptModels(readClassList, bpParameters, 
        stranded, min.readCount, min.readFractionByGene, 
        min.txScore.multiExon, min.txScore.singleExon, verbose) %>% data.table()
    combinedUnsplicedTranscripts[, confidenceType := "unsplicedNew"]
    combinedTranscripts <- as_tibble(rbindlist(list(combinedSplicedTranscripts,
        combinedUnsplicedTranscripts), fill = TRUE))
    return(combinedTranscripts)
}


#' combine spliced transcript models
#' @noRd
combineSplicedTranscriptModels <- function(readClassList, bpParameters, 
        min.readCount, min.readFractionByGene, min.txScore.multiExon, 
        min.txScore.singleExon, verbose){
    bpParameters$progressbar = FALSE
    options(scipen = 999) #maintain numeric basepair locations not sci.notfi.
    start.ptm <- proc.time()
    n_sample <- length(readClassList)
    nGroups = max(ceiling(n_sample/10),min(bpworkers(bpParameters), 
                                            round(n_sample/2)))
    indexList <- sample(rep(seq_len(nGroups), length.out=n_sample))
    indexList <- splitAsList(seq_len(n_sample), indexList)
    combinedFeatureTibbleList <- bplapply(seq_along(indexList), function(g){
        indexVec <- indexList[[g]]
        return(sequentialCombineFeatureTibble(readClassList[indexVec],
            indexVec, intraGroup = TRUE, 
            min.readCount = min.readCount, 
            min.readFractionByGene = min.readFractionByGene, 
            min.txScore.multiExon = min.txScore.multiExon,
            min.txScore.singleExon = min.txScore.singleExon))
    }, BPPARAM = bpParameters)
    combinedFeatureTibble <- 
        sequentialCombineFeatureTibble(combinedFeatureTibbleList, 
            indexList = NULL, intraGroup = FALSE) 
    combinedFeatureTibble <- updateStartEndReadCount(combinedFeatureTibble)
    end.ptm <- proc.time()
    if (verbose) message("combing spliced feature tibble objects across all ",
        "samples in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
    return(combinedFeatureTibble)
}

#' Sequentially combine feature tibbles 
#' @noRd
sequentialCombineFeatureTibble <- function(readClassList,
        indexList,intraGroup,min.readCount,min.readFractionByGene,
        min.txScore.multiExon, min.txScore.singleExon){
    combinedFeatureTibble <- NULL
    for (s in seq_along(readClassList)){
        combinedListNew <- readClassList[[s]]
        if(intraGroup){
            combinedListNew <- 
                extractFeaturesFromReadClassSE(readClassSe = combinedListNew,
                    sample_id = indexList[s], min.readCount = min.readCount,
                    min.readFractionByGene = min.readFractionByGene,
                    min.txScore.multiExon = min.txScore.multiExon,
                    min.txScore.singleExon = min.txScore.singleExon)
        }
        combinedFeatureTibble <- combineFeatureTibble(combinedFeatureTibble,
            combinedListNew, intraGroup)
    }
    return(combinedFeatureTibble)
}



#' Collapse the long-form sample chain rows to one row per intron chain
#' @noRd
updateStartEndReadCount <- function(combinedFeatureTibble){
    sampleChainRows <- rbindlist(combinedFeatureTibble$sampleChainRows)
    chains <- as.data.table(combinedFeatureTibble$chains)
    readCount <- sampleChainRows[, .(readCount = sum(readCount)), keyby = chain]$readCount
    return(data.table(
        start = readCountWeightedMedian(sampleChainRows$chain, sampleChainRows$start, sampleChainRows$readCount),
        end = readCountWeightedMedian(sampleChainRows$chain, sampleChainRows$end, sampleChainRows$readCount),
        readCount = readCount,
        chains[, .(intronStarts, intronEnds, chr, strand, maxTxScore,
            maxTxScore.noFit, NSampleReadCount, NSampleReadProp,
            NSampleTxScore)]))
}

#' Get the median of each chain's values, weighted by read count, without
#' interpolation: element floor(N/2)+1 of the values repeated by read count.
#' Chains must be numbered 1 to k.
#' @noRd
readCountWeightedMedian <- function(chain, value, readCount){
    ord <- order(chain, value, method = "radix")
    cum <- cumsum(as.numeric(readCount[ord]))
    last <- cumsum(tabulate(chain))
    prior <- c(0, cum[last])[seq_along(last)]
    return(value[ord][findInterval(prior + floor((cum[last] - prior) / 2) + 0.5,
        cum) + 1L])
}


#' Function to combine featureTibble and create the NSample variables.
#' Returns a list: chains has one row per intron chain, and sampleChainRows
#' holds the start, end and readCount of each chain in each sample in long
#' form, where sampleChainRows$chain is a row number of chains.
#' @noRd
combineFeatureTibble <- function(combinedFeatureTibble,
        featureTibbleSummarised, intraGroup = TRUE){
    new <- featureTibbleSummarised
    if (intraGroup) {
        new <- list(chains = select(new, -start, -end, -readCount),
            sampleChainRows = list(data.table(chain = seq_len(nrow(new)),
                start = new$start, end = new$end, readCount = new$readCount)))
    }
    if (is.null(combinedFeatureTibble)) return(new)
    combinedTable <- full_join(combinedFeatureTibble$chains,
            mutate(new$chains, newRow = row_number()),
            by = c('intronStarts', 'intronEnds', 'chr', 'strand'),
            suffix = c('.combined', '.new')) %>%
        mutate(NSampleReadCount = pmax0NA(NSampleReadCount.combined) +
                    pmax0NA(NSampleReadCount.new),
                NSampleReadProp = pmax0NA(NSampleReadProp.combined) +
                    pmax0NA(NSampleReadProp.new),
                NSampleTxScore = pmax0NA(NSampleTxScore.combined) +
                    pmax0NA(NSampleTxScore.new),
                maxTxScore = pmax(maxTxScore.combined,
                    maxTxScore.new, na.rm = TRUE),
                maxTxScore.noFit = pmax(maxTxScore.noFit.combined,
                    maxTxScore.noFit.new, na.rm = TRUE)) %>%
        select(intronStarts, intronEnds, chr, strand, NSampleReadCount,
            NSampleReadProp, NSampleTxScore, maxTxScore, maxTxScore.noFit,
            newRow)
    rowOfNew <- integer(nrow(new$chains))
    inNew <- which(!is.na(combinedTable$newRow))
    rowOfNew[combinedTable$newRow[inNew]] <- inNew
    newSampleChainRows <- lapply(new$sampleChainRows,
        function(x) x[, chain := rowOfNew[chain]])
    return(list(chains = select(combinedTable, -newRow),
        sampleChainRows = c(combinedFeatureTibble$sampleChainRows,
            newSampleChainRows)))
}
#' pmax replace NAs with 0
#' @noRd
pmax0NA <- function(vec){
    vec[is.na(vec)] <- 0
    return(pmax(vec))
}
#' pmin replace NAs with 0
#' @noRd
pmin0NA <- function(vec){
    vec[is.na(vec)] <- 0
    return(pmin(vec))
}
#' extract important features from readClassSe object for each sample
#' @noRd
extractFeaturesFromReadClassSE <- function(readClassSe, sample_id,
        min.readCount, min.readFractionByGene, 
        min.txScore.multiExon, min.txScore.singleExon){
    if (is.character(readClassSe)) 
        readClassSe <- readRDS(file = readClassSe)
    dimNames <- list(rownames(readClassSe), colnames(readClassSe))
    rowRangesSe <- rowRanges(readClassSe)
    rowData <- as_tibble(rowData(readClassSe)) %>% 
        mutate(start = unname(min(start(rowRangesSe))), 
                end= unname(max(end(rowRangesSe))))
    group_var <- c("intronStarts", "intronEnds", "chr", "strand")
    sum_var <- c("start","end","NSampleReadCount", "maxTxScore", 
                "maxTxScore.noFit", "readCount","NSampleReadProp",
                "NSampleTxScore")
    featureTibble <- rowData %>% 
        dplyr::select(chr = chr.rc, start, end, strand = strand.rc, 
            intronStarts, intronEnds, confidenceType, readCount, geneReadProp, 
            txScore, txScore.noFit, numExons) %>%
        filter(readCount >= 1, # only use readCount>1 and highconfidence reads
            confidenceType == "highConfidenceJunctionReads") %>% 
        mutate(NSampleReadCount = (readCount >= min.readCount), 
            # number of samples passed read count criteria
            NSampleReadProp = (geneReadProp >= min.readFractionByGene),
            NSampleTxScore = ((txScore > min.txScore.multiExon & numExons >= 2) |
            (txScore > min.txScore.singleExon & numExons == 1)), 
            maxTxScore = txScore, maxTxScore.noFit = txScore.noFit) %>%
        select(all_of(c(group_var, sum_var))) 
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
            min.readFractionByGene, min.txScore.multiExon,
            min.txScore.singleExon, verbose){
        start.ptm <- proc.time()
        bpParameters$progressbar = FALSE
        newUnsplicedSeList <- 
            bplapply(seq_along(readClassList), function(sample_id)
                extractNewUnsplicedRanges(readClassSe = 
                readClassList[[sample_id]], sample_id = sample_id), 
                BPPARAM = bpParameters)
        end.ptm <- proc.time()
        if (verbose) message("extract new unspliced ranges object for all ",
        "samples in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
        rangesList <- bplapply(newUnsplicedSeList, function(newUnsplicedSe){
            rr <- unlist(rowRanges(newUnsplicedSe))
            rr$row_id <- names(rr)
            return(rr)
        }, BPPARAM = bpParameters)
        colDataNames <-unlist(lapply(newUnsplicedSeList, colnames))
        start.ptm <- proc.time()
        combinedNewUnsplicedSe <- reduceUnsplicedRanges(rangesList, stranded)
        end.ptm <- proc.time()
        if (verbose) message("reduce new unspliced ranges object across all ",
        "samples in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
        start.ptm <- proc.time()
        combinedUnsplicedTibble <- 
            makeUnsplicedTibble(combinedNewUnsplicedSe,newUnsplicedSeList, 
                colDataNames, min.readCount, min.readFractionByGene,
                min.txScore.multiExon, min.txScore.singleExon, bpParameters)
        end.ptm <- proc.time()
        if (verbose) message("combine new unspliced tibble object across all ",
        "samples in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
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
                                rowData$readCount > 1 &
                                (!rowData$equal))]
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
#' @noRd
makeUnsplicedTibble <- function(combinedNewUnsplicedSe,newUnsplicedSeList,
        colDataNames,min.readCount, min.readFractionByGene,
        min.txScore.multiExon, min.txScore.singleExon, bpParameters){
        bpParameters$progressbar = FALSE
    newUnsplicedTibble <- as_tibble(combinedNewUnsplicedSe) %>%
        rename(chr = seqnames) %>% select(chr, start, end, strand, row_id) %>%
        separate_rows(row_id, sep = "\\+") 
    rowDataCombined <-
        do.call("rbind",bplapply(newUnsplicedSeList, function(newUnsplicedSe) {
            rr <- rowData(newUnsplicedSe[intersect(rownames(newUnsplicedSe), 
                                        newUnsplicedTibble$row_id)])
            rr <- as_tibble(rr) %>% select(confidenceType,readCount, 
                    geneReadProp, txScore, txScore.noFit) %>%
                mutate(row_id = rownames(rr))
            return(rr)
        } , BPPARAM = bpParameters))
    newUnsplicedTibble <- newUnsplicedTibble %>% 
        left_join(rowDataCombined, by =  "row_id") %>%
        mutate(readCount_tmp = readCount) %>%
        group_by(chr,strand, start, end) %>%
        summarise(readCount = sum(readCount),
                  maxTxScore = weighted.mean(txScore, readCount_tmp),
                  maxTxScore.noFit = weighted.mean(txScore.noFit, readCount_tmp),
                  NSampleReadCount = sum(readCount_tmp >= min.readCount), 
                  NSampleReadProp = sum(geneReadProp >= 
                                          min.readFractionByGene),
                  NSampleTxScore = sum(txScore > min.txScore.singleExon))
    
    return(newUnsplicedTibble)
}
