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
#'
#' Dispatches to a sparse reduction, falling back to the original dense
#' implementation when the sparse path cannot reproduce it exactly (see
#' combineSplicedTranscriptModelsSparse). Both return the same 12 columns,
#' in the same order, with the same row order.
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
    combinedFeatureTibble <- NULL
    ## With fewer than 3 samples nGroups is 1, so the dense path performs no
    ## full_join at all and its NSample* columns stay logical while maxTxScore
    ## keeps NAs. Reproducing that column-type quirk is not worth it; defer.
    if (n_sample >= 3)
        combinedFeatureTibble <- combineSplicedTranscriptModelsSparse(
            readClassList, indexList, bpParameters, min.readCount,
            min.readFractionByGene, min.txScore.multiExon,
            min.txScore.singleExon)
    if (is.null(combinedFeatureTibble))
        combinedFeatureTibble <- combineSplicedTranscriptModelsDense(
            readClassList, indexList, bpParameters, min.readCount,
            min.readFractionByGene, min.txScore.multiExon,
            min.txScore.singleExon)
    end.ptm <- proc.time()
    if (verbose) message("combing spliced feature tibble objects across all ",
        "samples in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
    return(combinedFeatureTibble)
}

#' combine spliced transcript models by joining per-sample columns
#'
#' The original implementation: builds one wide table carrying a start, end and
#' readCount column per sample, then reduces it row-wise.
#' @noRd
combineSplicedTranscriptModelsDense <- function(readClassList, indexList,
        bpParameters, min.readCount, min.readFractionByGene,
        min.txScore.multiExon, min.txScore.singleExon){
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
    return(updateStartEndReadCount(combinedFeatureTibble))
}

#' combine spliced transcript models without materialising per-sample columns
#'
#' The dense path builds a table whose rows are distinct intron chains and whose
#' columns are 3n+9 for n samples. A read class is observed in only a few
#' samples, so that table is mostly NA, and the join fold holds a second copy of
#' it while the first is live.
#'
#' This stores only the populated (chain, sample) cells in long form and reduces
#' them by chain, so memory scales with observations rather than rows x samples.
#'
#' Returns NULL when it cannot guarantee an identical result, so the caller can
#' fall back to the dense implementation.
#' @noRd
combineSplicedTranscriptModelsSparse <- function(readClassList, indexList,
        bpParameters, min.readCount, min.readFractionByGene,
        min.txScore.multiExon, min.txScore.singleExon){
    res <- bplapply(seq_along(indexList), function(g){
        indexVec <- as.integer(indexList[[g]])
        return(sparseGroupFeatures(readClassList[indexVec], indexVec,
            min.readCount = min.readCount,
            min.readFractionByGene = min.readFractionByGene,
            min.txScore.multiExon = min.txScore.multiExon,
            min.txScore.singleExon = min.txScore.singleExon))
    }, BPPARAM = bpParameters)
    if (any(vapply(res, is.null, logical(1)))) return(NULL)
    ## Global first-appearance key ids, over the same group-major traversal the
    ## dense path uses, so row order is preserved.
    allk <- rbindlist(lapply(res, `[[`, "gk"), idcol = "g")
    for (i in seq_along(res)) res[[i]]$gk <- NULL
    allk[, gi := .GRP, by = c("intronStarts", "intronEnds", "chr", "strand")]
    nkey <- max(allk$gi)
    keyDT <- allk[!duplicated(gi), .(intronStarts, intronEnds, chr, strand)]
    if (nrow(keyDT) != nkey) return(NULL)
    agg <- allk[, .(nsrc = sum(nsrc), nsrp = sum(nsrp), nstx = sum(nstx),
        mts = max(mts), mtsnf = max(mtsnf)), by = gi]
    if (!identical(agg$gi, seq_len(nkey))) return(NULL)
    cellCount <- vapply(seq_along(res), function(i) nrow(res[[i]]$cells), 0L)
    gvec <- allk$g
    giVec <- allk$gi
    rm(allk)
    GI <- integer(sum(cellCount)); ST <- integer(sum(cellCount))
    EN <- integer(sum(cellCount)); RC <- integer(sum(cellCount))
    pos <- 1L
    for (i in seq_along(res)){
        cl <- res[[i]]$cells
        keyMap <- giVec[gvec == i] # group-local id -> global id
        if (nrow(cl)){
            sl <- pos:(pos + nrow(cl) - 1L)
            GI[sl] <- keyMap[cl$li]
            ST[sl] <- cl$start; EN[sl] <- cl$end; RC[sl] <- cl$readCount
            pos <- pos + nrow(cl)
        }
        res[[i]]$cells <- NULL
    }
    rm(res, gvec, giVec)
    if (anyNA(ST) || anyNA(EN) || anyNA(RC)) return(NULL)
    ord <- order(GI, method = "radix")
    GI <- GI[ord]; ST <- ST[ord]; EN <- EN[ord]; RC <- RC[ord]
    rm(ord)
    reduced <- reduceSparseCells(GI, ST, EN, RC, nkey)
    if (is.null(reduced)) return(NULL)
    ## A missing score is stored as -Inf while reducing so max() stays valid.
    mts <- agg$mts; mts[is.infinite(mts) & mts < 0] <- NA_real_
    mtsnf <- agg$mtsnf; mtsnf[is.infinite(mtsnf) & mtsnf < 0] <- NA_real_
    ## pmax() on a logical returns double, so the dense path's NSample* columns
    ## are double once any join has happened. Match that.
    return(data.table(start = reduced$start, end = reduced$end,
        readCount = reduced$readCount, intronStarts = keyDT$intronStarts,
        intronEnds = keyDT$intronEnds, chr = keyDT$chr,
        strand = keyDT$strand, maxTxScore = as.numeric(mts),
        maxTxScore.noFit = as.numeric(mtsnf),
        NSampleReadCount = as.numeric(agg$nsrc),
        NSampleReadProp = as.numeric(agg$nsrp),
        NSampleTxScore = as.numeric(agg$nstx)))
}

#' summarise one group of samples into a key table and its populated cells
#' @noRd
sparseGroupFeatures <- function(readClassList, indexVec, min.readCount,
        min.readFractionByGene, min.txScore.multiExon, min.txScore.singleExon){
    tabs <- vector("list", length(indexVec))
    for (s in seq_along(indexVec)){
        featureTibble <- extractFeaturesFromReadClassSE(
            readClassSe = readClassList[[s]], sample_id = indexVec[s],
            min.readCount = min.readCount,
            min.readFractionByGene = min.readFractionByGene,
            min.txScore.multiExon = min.txScore.multiExon,
            min.txScore.singleExon = min.txScore.singleExon)
        setDT(featureTibble)
        ## One row per key per sample is what makes per-key counters equal to
        ## per-row counters; a duplicate would fan out in the dense join.
        if (anyDuplicated(featureTibble,
            by = c("intronStarts", "intronEnds", "chr", "strand")))
            return(NULL)
        tabs[[s]] <- featureTibble
    }
    big <- rbindlist(tabs)
    rm(tabs)
    big[, li := .GRP, by = c("intronStarts", "intronEnds", "chr", "strand")]
    gk <- big[, .(intronStarts = intronStarts[1L],
        intronEnds = intronEnds[1L], chr = chr[1L], strand = strand[1L],
        nsrc = sum(NSampleReadCount), nsrp = sum(NSampleReadProp),
        nstx = sum(NSampleTxScore, na.rm = TRUE), # NA counts as 0, as pmax0NA
        mts = if (all(is.na(maxTxScore))) -Inf else
            max(maxTxScore, na.rm = TRUE),
        mtsnf = if (all(is.na(maxTxScore.noFit))) -Inf else
            max(maxTxScore.noFit, na.rm = TRUE)), by = li]
    if (!identical(gk$li, seq_len(nrow(gk)))) return(NULL)
    gk[, li := NULL]
    return(list(gk = gk, cells = big[, .(li, start, end, readCount)]))
}

#' reduce long-form cells to one start, end and readCount per key
#' @noRd
reduceSparseCells <- function(GI, ST, EN, RC, nkey, cellChunk = 2e7){
    readCountSum <- numeric(nkey)
    summed <- rowsum(as.numeric(RC), GI, reorder = FALSE)
    readCountSum[as.integer(rownames(summed))] <- summed[, 1L]
    rm(summed)
    if (any(readCountSum > .Machine$integer.max)) return(NULL)
    startOut <- rep(Inf, nkey)
    endOut <- rep(Inf, nkey)
    ## Chunk on key boundaries so no key is split across calls.
    brk <- c(0L, which(GI[-1L] != GI[-length(GI)]), length(GI))
    step <- max(1L, as.integer(cellChunk))
    b0 <- 1L
    while (b0 <= length(brk) - 1L){
        b1 <- b0
        while (b1 < length(brk) - 1L && (brk[b1 + 2L] - brk[b0]) <= step)
            b1 <- b1 + 1L
        idx <- (brk[b0] + 1L):brk[b1 + 1L]
        keys <- GI[idx]
        loK <- keys[1L]; hiK <- keys[length(keys)]
        startOut[loK:hiK] <- upperMedianByGroup(keys - loK + 1L, ST[idx],
            RC[idx], hiK - loK + 1L)
        endOut[loK:hiK] <- upperMedianByGroup(keys - loK + 1L, EN[idx],
            RC[idx], hiK - loK + 1L)
        b0 <- b1 + 1L
    }
    if (any(is.infinite(startOut)) || any(is.infinite(endOut))) return(NULL)
    return(list(start = as.integer(startOut), end = as.integer(endOut),
        readCount = as.integer(readCountSum)))
}

#' readCount-weighted median per group, without expanding the values
#'
#' readCountWeightedMedian() repeats each value by its read count, takes the
#' type 7 median and snaps up to the nearest observed value. Because the weights
#' are integers that is exactly the element at floor(N/2)+1 of the expanded
#' sorted vector, which can be located from the cumulative weights alone.
#' @noRd
upperMedianByGroup <- function(groupIndex, values, weights, ngroup){
    out <- rep(Inf, ngroup)
    if (!length(groupIndex)) return(out)
    ord <- order(groupIndex, values, method = "radix")
    groupIndex <- groupIndex[ord]
    values <- values[ord]
    cumWeight <- cumsum(as.numeric(weights[ord]))
    starts <- c(1L, which(groupIndex[-1L] !=
        groupIndex[-length(groupIndex)]) + 1L)
    priorWeight <- c(0, cumWeight)[starts]
    lastIdx <- c(starts[-1L] - 1L, length(cumWeight))
    groupTotal <- cumWeight[lastIdx] - priorWeight
    target <- priorWeight + floor(groupTotal / 2) + 1
    hit <- findInterval(target - 0.5, cumWeight) + 1L
    nonEmpty <- groupTotal > 0
    out[groupIndex[starts][nonEmpty]] <- values[hit[nonEmpty]]
    return(out)
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
            combinedListNew, index = indexList[s], intraGroup)
    }
    return(combinedFeatureTibble)
}



#' @noRd 
updateStartEndReadCount <- function(combinedFeatureTibble){
    setDT(combinedFeatureTibble)
    combinedFeatureTibble[, rowID := .I]
    
    colNames <- colnames(combinedFeatureTibble)
    readCountCols <- sort(colNames[grep("^readCount", colNames)]) # to make sure it's ordered by sample name
    startCols <- sort(colNames[grep("^start", colNames)])
    endCols <- sort(colNames[grep("^end", colNames)])
    
    startEndDt <- combinedFeatureTibble[, 
        .(start = readCountWeightedMedian(.SD,x,y),
        end = readCountWeightedMedian(.SD,z,y),
        readCount = sum(.SD[,y], na.rm = TRUE)),
        by = rowID,  env = I(list(x = startCols, y = readCountCols,z = endCols))]

    combinedFeatureTibble <- startEndDt[combinedFeatureTibble[,.(intronStarts, intronEnds, chr, strand, maxTxScore, 
                                                                 maxTxScore.noFit, NSampleReadCount, NSampleReadProp, 
                                                                 NSampleTxScore, rowID)], on = "rowID"]
    combinedFeatureTibble[, rowID := NULL]
    return(combinedFeatureTibble)
}

#' Function to get median value without interpolation using certain column names
#' @noRd
readCountWeightedMedian <- function(dt, valuevar, timesvar){
    sortVector <- rep(na.omit(unlist(dt[,..valuevar])), 
                times = as.integer(na.omit(unlist(dt[,..timesvar]))))
    return(min(sortVector[sortVector>=quantile(sortVector, probs = 0.5)]))
}


#' Function to combine featureTibble and create the NSample variables 
#' @noRd
combineFeatureTibble <- function(combinedFeatureTibble,
        featureTibbleSummarised, index=1, intraGroup = TRUE){ 
    if (is.null(combinedFeatureTibble)) { 
        combinedTable <- featureTibbleSummarised %>% 
            select(intronStarts, intronEnds, chr, strand, maxTxScore, 
            maxTxScore.noFit, NSampleReadCount, NSampleReadProp,NSampleTxScore, 
            starts_with('start'), starts_with('end'), starts_with('readCount'))
    } else { 
        combinedTable = full_join(combinedFeatureTibble, 
            featureTibbleSummarised, by = c('intronStarts', 'intronEnds', 'chr',
            'strand'), suffix=c('.combined','.new')) %>% 
            mutate(NSampleReadCount=pmax0NA(NSampleReadCount.combined) + 
                        pmax0NA(NSampleReadCount.new), 
                    NSampleReadProp = pmax0NA(NSampleReadProp.combined) + 
                        pmax0NA(NSampleReadProp.new), 
                    NSampleTxScore = pmax0NA(NSampleTxScore.combined) + 
                        pmax0NA(NSampleTxScore.new),
                    maxTxScore = pmax(maxTxScore.combined, 
                        maxTxScore.new, na.rm = TRUE),
                    maxTxScore.noFit = pmax(maxTxScore.noFit.combined, 
                        maxTxScore.noFit.new, na.rm = TRUE)) %>% 
            select(intronStarts, intronEnds, chr, strand,
            NSampleReadCount, NSampleReadProp, NSampleTxScore, maxTxScore, 
            maxTxScore.noFit, starts_with('start'), starts_with('end'), 
            starts_with('readCount')) 
    } 
    if(intraGroup) 
        combinedTable <- 
            rename_with(combinedTable, ~gsub('^(end|start|readCount)$',
                                            paste0('\\1\\.',index), .x)) 
    return(combinedTable) 
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
