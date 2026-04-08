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
    min.intronChainScore.multiExon, min.intronChainScore.singleExon, bpParameters ,verbose){
    combinedSplicedTranscripts <- 
        combineSplicedTranscriptModels(readClassList, bpParameters, 
        min.readCount, min.readFractionByGene, 
        min.intronChainScore.multiExon, min.intronChainScore.singleExon, verbose) %>% data.table()
    combinedSplicedTranscripts[,confidenceType := "highConfidenceJunctionReads"]
    # when single exon min score is greater than 1, skip unspliced transcripts combination
    # this is a very customized config, useful when data is very big 
    if (min.intronChainScore.singleExon > 1) 
        return(combinedSplicedTranscripts)
    combinedUnsplicedTranscripts <- 
        combineUnsplicedTranscriptModels(readClassList, bpParameters, 
        stranded, min.readCount, min.readFractionByGene, 
        min.intronChainScore.multiExon, min.intronChainScore.singleExon, verbose) %>% data.table()
    combinedUnsplicedTranscripts[, confidenceType := "unsplicedNew"]
    combinedTranscripts <- as_tibble(rbindlist(list(combinedSplicedTranscripts,
        combinedUnsplicedTranscripts), fill = TRUE))
    return(combinedTranscripts)
}


#' combine spliced transcript models
#' @noRd
combineSplicedTranscriptModels <- function(readClassList, bpParameters, 
        min.readCount, min.readFractionByGene, min.intronChainScore.multiExon, 
        min.intronChainScore.singleExon, verbose){
    bpParameters$progressbar <- FALSE
    options(scipen = 999) #maintain numeric basepair locations not sci.notfi.
    start.ptm <- proc.time()
    n_sample <- length(readClassList)
    nGroups <- max(ceiling(n_sample/10),min(bpworkers(bpParameters), 
                                            round(n_sample/2)))
    indexList <- sample(rep(seq_len(nGroups), length.out=n_sample))
    indexList <- splitAsList(seq_len(n_sample), indexList)
    combinedFeatureTibbleList <- bplapply(seq_along(indexList), function(g){
        indexVec <- indexList[[g]]
        return(sequentialCombineFeatureTibble(readClassList[indexVec],
            indexVec, intraGroup = TRUE, 
            min.readCount = min.readCount, 
            min.readFractionByGene = min.readFractionByGene, 
            min.intronChainScore.multiExon = min.intronChainScore.multiExon,
            min.intronChainScore.singleExon = min.intronChainScore.singleExon))
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
        indexList, intraGroup, min.readCount, min.readFractionByGene,
        min.intronChainScore.multiExon, min.intronChainScore.singleExon){
    combinedFeatureTibble <- NULL
    for (s in seq_along(readClassList)){
        combinedListNew <- readClassList[[s]]
        if(intraGroup){
            combinedListNew <- 
                extractFeaturesFromReadClassSE(readClassSe = combinedListNew,
                    sample_id = indexList[s], min.readCount = min.readCount,
                    min.readFractionByGene = min.readFractionByGene,
                    min.intronChainScore.multiExon = min.intronChainScore.multiExon,
                    min.intronChainScore.singleExon = min.intronChainScore.singleExon)
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
    startCols <- sort(colNames[grep("^start\\.[0-9]", colNames)])
    endCols <- sort(colNames[grep("^end\\.[0-9]", colNames)])

    startEndDt <- combinedFeatureTibble[, 
        .(start = readCountWeightedMedian(.SD,x,y),
        end = readCountWeightedMedian(.SD,z,y),
        readCount = sum(.SD[,y], na.rm = TRUE)),
        by = rowID,  env = I(list(x = startCols, y = readCountCols, z = endCols))]
    combinedFeatureTibble <- startEndDt[combinedFeatureTibble[,.(intronStarts, intronEnds, chr, strand, maxIntronChainScore, maxIntronChainScore.noFit, 
                                                                 firstExonGroup, lastExonGroup, sampleTesId, startRegionId, endRegionId, 
                                                                 compatible, equal,
                                                                 NSampleReadCount, NSampleReadProp, 
                                                                 NSampleIntronChainScore, rowID)], on = "rowID"]
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
            select(intronStarts, intronEnds, chr, strand, firstExonGroup, lastExonGroup, sampleTesId, startRegionId, endRegionId, compatible, equal,
            maxIntronChainScore, maxIntronChainScore.noFit, 
            NSampleReadCount, NSampleReadProp,NSampleIntronChainScore, 
            starts_with('start'), starts_with('end'), starts_with('readCount'))
    } else { 
        combinedTable <- full_join(combinedFeatureTibble, 
            featureTibbleSummarised, by = c('intronStarts', 'intronEnds', 'chr',
            'strand', 'firstExonGroup', 'lastExonGroup', "startRegionId", "sampleTesId","endRegionId", "compatible", "equal"), suffix=c('.combined','.new')) %>% 
            mutate(NSampleReadCount=pmax0NA(NSampleReadCount.combined) + 
                        pmax0NA(NSampleReadCount.new), 
                    NSampleReadProp = pmax0NA(NSampleReadProp.combined) + 
                        pmax0NA(NSampleReadProp.new), 
                    NSampleIntronChainScore = pmax0NA(NSampleIntronChainScore.combined) + 
                        pmax0NA(NSampleIntronChainScore.new),
                    maxIntronChainScore = pmax(maxIntronChainScore.combined,
                                                  maxIntronChainScore.new, na.rm = TRUE),
                    maxIntronChainScore.noFit = pmax(maxIntronChainScore.noFit.combined,
                                                        maxIntronChainScore.noFit.new, na.rm = TRUE)) %>% 
            select(intronStarts, intronEnds, chr, strand,
            NSampleReadCount, NSampleReadProp, NSampleIntronChainScore, 
            maxIntronChainScore, maxIntronChainScore.noFit, 
            starts_with('start'), starts_with('end'), 
            starts_with('readCount'), firstExonGroup, lastExonGroup, sampleTesId, startRegionId, endRegionId, compatible, equal) 
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
        min.intronChainScore.multiExon, min.intronChainScore.singleExon){
    if (is.character(readClassSe)) 
        readClassSe <- readRDS(file = readClassSe)
    dimNames <- list(rownames(readClassSe), colnames(readClassSe))
    rowRangesSe <- rowRanges(readClassSe)
    rowData <- as_tibble(rowData(readClassSe)) %>% 
        mutate(start = unname(min(start(rowRangesSe))), 
                end= unname(max(end(rowRangesSe))))

    group_var <- c("intronStarts", "intronEnds", "chr", "strand", "firstExonGroup", "lastExonGroup", "sampleTesId",
        "startRegionId", "endRegionId", "compatible", "equal")
    sum_var <- c("start","end","NSampleReadCount", 
                "maxIntronChainScore", "maxIntronChainScore.noFit",
                "readCount", "NSampleReadProp",
                "NSampleIntronChainScore")
    featureTibble <- rowData %>% 
        dplyr::select(chr = chr.rc, start, end, strand = strand.rc, firstExonGroup, lastExonGroup, sampleTesId,
            startRegionId, endRegionId, compatible, equal,
            intronStarts, intronEnds, confidenceType, readCount, geneReadProp, 
            intronChainScore, intronChainScore.noFit, numExons) %>%
        filter(readCount >= 1, # only use readCount>1 and highconfidence reads
            confidenceType == "highConfidenceJunctionReads") %>% 
        mutate(NSampleReadCount = (readCount >= min.readCount), 
            # number of samples passed read count criteria
            NSampleReadProp = (geneReadProp >= min.readFractionByGene),
            NSampleIntronChainScore = ((intronChainScore > min.intronChainScore.multiExon & numExons >= 2) |
            (intronChainScore > min.intronChainScore.singleExon & numExons == 1)), 
            maxIntronChainScore = intronChainScore, maxIntronChainScore.noFit = intronChainScore.noFit) %>%
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
            min.readFractionByGene, min.intronChainScore.multiExon,
            min.intronChainScore.singleExon, verbose){
        start.ptm <- proc.time()
        bpParameters$progressbar <- FALSE
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
                min.intronChainScore.multiExon, min.intronChainScore.singleExon, bpParameters)
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
        min.intronChainScore.multiExon, min.intronChainScore.singleExon, bpParameters){
        bpParameters$progressbar <- FALSE
    newUnsplicedTibble <- as_tibble(combinedNewUnsplicedSe) %>%
        rename(chr = seqnames) %>% select(chr, start, end, strand, row_id) %>%
        separate_rows(row_id, sep = "\\+") 
    rowDataCombined <-
        do.call("rbind",bplapply(newUnsplicedSeList, function(newUnsplicedSe) {
            rr <- rowData(newUnsplicedSe[intersect(rownames(newUnsplicedSe), 
                                        newUnsplicedTibble$row_id)])
            rr <- as_tibble(rr) %>% select(confidenceType,readCount, 
                    geneReadProp, intronChainScore, intronChainScore.noFit) %>%
                mutate(row_id = rownames(rr))
            return(rr)
        } , BPPARAM = bpParameters))
    newUnsplicedTibble <- newUnsplicedTibble %>% 
        left_join(rowDataCombined, by =  "row_id") %>%
        mutate(readCount_tmp = readCount) %>%
        group_by(chr,strand, start, end) %>%
        summarise(readCount = sum(readCount),
                  maxIntronChainScore = weighted.mean(intronChainScore, readCount_tmp),
                  maxIntronChainScore.noFit = weighted.mean(intronChainScore.noFit, readCount_tmp),
                  NSampleReadCount = sum(readCount_tmp >= min.readCount), 
                  NSampleReadProp = sum(geneReadProp >= 
                                          min.readFractionByGene),
                  NSampleIntronChainScore = sum(intronChainScore > min.intronChainScore.singleExon))
    return(newUnsplicedTibble)
}
