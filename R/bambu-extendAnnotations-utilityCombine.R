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
                                       min.txScore.multiExon, min.txScore.singleExon, verbose) %>% data.table()
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
                                           min.readCount, min.readFractionByGene, 
                                           min.txScore.multiExon, min.txScore.singleExon,
                                           verbose){
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
    if (verbose) message("combing spliced feature tibble objects across all
        samples in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
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
                                                      combinedListNew, index = indexList[s], intraGroup)
    }
    return(combinedFeatureTibble)
}



#' @noRd 
updateStartEndReadCount <- function(combinedFeatureTibble){
    combinedFeatureTibble <- combinedFeatureTibble %>% mutate(rowID = row_number())
    
    startEndCountTibble <- combinedFeatureTibble %>% 
        select(rowID, starts_with("start"),starts_with("end"),starts_with("readCount")) %>%
        tidyr::pivot_longer(c(starts_with("start"),starts_with("end"),starts_with("readCount")),
                            names_to = c(".value","set"),
                            names_pattern = "(.*)\\.(.)") %>%
        group_by(rowID) %>% 
        mutate(sumReadCount = sum(readCount,na.rm = TRUE))
    
    startTibble <- select(startEndCountTibble, rowID, start, readCount, sumReadCount) %>% 
        arrange(start) %>%
        filter(cumsum(readCount)/sumReadCount>=0.5) %>% 
        filter(row_number()==1)
    endTibble <- select(startEndCountTibble, rowID, end, readCount, sumReadCount) %>% 
        arrange(end) %>% 
        filter(cumsum(readCount)/sumReadCount>=0.5) %>% 
        filter(row_number()==1)
    
    combinedFeatureTibble <- combinedFeatureTibble %>% 
        dplyr::select(intronStarts, intronEnds, chr, strand, maxTxScore, NSampleReadCount, 
                      NSampleReadProp, NSampleTxScore, rowID) %>%
        full_join(select(startTibble, rowID, start), by = "rowID") %>% 
        full_join(select(endTibble, rowID, end, readCount=sumReadCount), by = "rowID") %>%
        select(-rowID)
    return(combinedFeatureTibble)
}



#' Function to combine featureTibble and create the NSample variables 
#' @noRd
combineFeatureTibble <- function(combinedFeatureTibble,
                                 featureTibbleSummarised, index=1, intraGroup = TRUE){ 
    if (is.null(combinedFeatureTibble)) { 
        combinedTable <- featureTibbleSummarised %>% 
            select(intronStarts, intronEnds, chr, strand, maxTxScore, NSampleReadCount,
                   NSampleReadProp,NSampleTxScore, starts_with('start'),
                   starts_with('end'), starts_with('readCount'))
    } else { 
        combinedTable = full_join(combinedFeatureTibble, 
                                  featureTibbleSummarised, by = c('intronStarts',
                                                                  'intronEnds', 'chr', 'strand'),
                                  suffix=c('.combined','.new')) %>% 
            mutate(NSampleReadCount=pmax0NA(NSampleReadCount.combined) + 
                       pmax0NA(NSampleReadCount.new), 
                   NSampleReadProp = pmax0NA(NSampleReadProp.combined) + 
                       pmax0NA(NSampleReadProp.new), 
                   NSampleTxScore = pmax0NA(NSampleTxScore.combined) + 
                       pmax0NA(NSampleTxScore.new),
                    maxTxScore = pmax(maxTxScore.combined, 
                        maxTxScore.new, na.rm = TRUE)) %>% 
            select(intronStarts, intronEnds, chr, strand,
            NSampleReadCount, NSampleReadProp, NSampleTxScore, maxTxScore, 
            starts_with('start'), starts_with('end'), starts_with('readCount')) 
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
                 "readCount","NSampleReadProp","NSampleTxScore")
    featureTibble <- rowData %>% 
        #filter(!equal) %>% # filter not compatible ones, i.e., overlapping with annotations?? if we are going to include subset tx, then can we still do the filtering?? maybe not equal but can be compatible 
        dplyr::select(chr = chr.rc, start, end, equal,
                      strand = strand.rc, intronStarts, intronEnds, confidenceType,
                      readCount, geneReadProp, txScore, numExons) %>%
        filter(readCount > 1, # only use readCount>1 and highconfidence reads
               confidenceType == "highConfidenceJunctionReads") %>% 
        mutate(NSampleReadCount = (readCount >= min.readCount), 
               # number of samples passed read count criteria
               NSampleReadProp = (geneReadProp >= min.readFractionByGene),
               NSampleTxScore = ((txScore >= min.txScore.multiExon & numExons >= 2) |
                (txScore >= min.txScore.singleExon & numExons == 1)), 
                maxTxScore = txScore) %>%
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
        colDataNames <-unlist(lapply(newUnsplicedSeList, colnames))
        start.ptm <- proc.time()
        combinedNewUnsplicedSe <- reduceUnsplicedRanges(rangesList, stranded)
        end.ptm <- proc.time()
        if (verbose) message("reduce new unspliced ranges object across all
        samples in ", round((end.ptm - start.ptm)[3] / 60, 1)," mins.")
        start.ptm <- proc.time()
        combinedUnsplicedTibble <- 
            makeUnsplicedTibble(combinedNewUnsplicedSe,newUnsplicedSeList, 
                                colDataNames, min.readCount, min.readFractionByGene,
                                min.txScore.multiExon, min.txScore.singleExon, bpParameters)
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
makeUnsplicedTibble <- function(combinedNewUnsplicedSe,newUnsplicedSeList,
                                colDataNames,min.readCount, min.readFractionByGene,
                                min.txScore.multiExon, min.txScore.singleExon, 
                                bpParameters){
    newUnsplicedTibble <- as_tibble(combinedNewUnsplicedSe) %>%
        rename(chr = seqnames) %>%
        select(chr, start, end, strand, row_id) %>%
        separate_rows(row_id, sep = "\\+") 
    rowDataCombined <-
        do.call("rbind",bplapply(newUnsplicedSeList, function(newUnsplicedSe) {
            rr <- rowData(newUnsplicedSe[intersect(rownames(newUnsplicedSe), 
                                                   newUnsplicedTibble$row_id)])
            rr <- as_tibble(rr) %>% select(confidenceType,
                                           readCount, geneReadProp, txScore) %>%
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
                  txScore = weighted.mean(txScore, readCount_tmp)) %>%
        group_by(chr, strand, start, end) %>% 
        summarise(readCount = sum(readCount),
                    NSampleReadCount = sum(readCount >= min.readCount), 
                    NSampleReadProp = sum(geneReadProp >= 
                                            min.readFractionByGene),
                    NSampleTxScore = sum(txScore >= min.txScore.singleExon))
    return(newUnsplicedTibble)
}
