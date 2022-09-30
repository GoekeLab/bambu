library(dplyr)
library(purrr)
library(ROCit)
library(Biostrings)

# **** This is the code script to predict transcription state site of novel transcript. Refer
# to the slides for some of the preliminary analysis **** 

# Load the necessary data 
annotations <- readRDS("./output/annotations.rds") 
fastaFile <- readDNAStringSet("./dataset/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa")                                                       
RCFile <- readRDS("./output/A549_directRNA_replicate4.rds") # read class file, run bambu to get it. 

# Helper function
onlyAnnotatedRCWithMoreThanOneExon <- function(rcFile){
    # Input: rcFile = read class file
    # Output: returns a Summarized Experiment object which includes only annotated read 
    # class and excludes single exon read class. 
  
    rcFile <- rcFile[!is.na(rowData(rcFile)$equal)]
    rcFile <- rcFile[rowData(rcFile)$equal == TRUE]
    rcFile <- rcFile[rowData(rcFile)$numExons > 1]
    
    return(rcFile)
}

getMetadata <- function(rcFile){
  # Input: rcFile = read class file 
  # Output: returns a dataframe that contains metadata information for each read:
  #         peaks - the start site of each read
  #         readIds - the read ID as mentioned in the rcFile 
  #         rc - the read class name that the read belongs to 
  #         GENEID - the gene name that the read belongs to 
  #         chr - the chromosome name that the read belongs to 
  #         readStrands - positive strand if TRUE, negative strand if FALSE 
  
  rcFileDf <- mcols(rcFile)[, c("starts", "ends", "readIds", "GENEID", "chr.rc", "readStrands", "strand.rc")]
  peaksMetadata <- data.frame(rcFileDf) %>% 
    tibble::rownames_to_column("rc") %>% 
    as.data.frame()  %>%
    tidyr::unnest_longer(c(starts, ends, readIds, readStrands)) %>% 
    arrange(readIds) %>% 
    mutate(peaks = ifelse(readStrands == TRUE, starts, ends)) %>% 
    rename("chr" = "chr.rc") %>% 
    select(peaks, readIds, rc, GENEID, chr, readStrands)
  
  return(peaksMetadata)
}

getAnnotatedStartSite <- function(rcFile, rc){
  # Input: rcFile = read class file 
  #        rc  = a selected vector of read class name
  # Output: returns a dataframe that gives the peaks of the annotated transcripts associated with
  #         the read class. Remark: one read class could have multiple match with the annotation. 
  
  rcInfo <- data.frame(ID = seq_along(rc),rc)
  tbl <- data.frame(findOverlaps(cutStartEndFromGrangesList(rowRanges(rcFile)[rc]), 
                                 cutStartEndFromGrangesList(annotations), 
                                 type = "equal")) %>% 
    left_join(rcInfo, by = c("queryHits" = "ID")) %>% 
    mutate(annotatedStarts = as.vector(unlist(min(start(annotations[subjectHits])))), 
           annotatedEnds = as.vector(unlist(max(end(annotations[subjectHits])))), 
           annotatedStrands = as.vector(unlist(runValue(strand(annotations[subjectHits])))), 
           annotatedTss = ifelse(annotatedStrands == "+", annotatedStarts, annotatedEnds))
  
  return(tbl)
}

getFeatureData <- function(metadata, rcFile, neighbourhoodUpstreamWS, neighbourhoodDownstreamWS, 
                           labelWS, upstreamWS, downstreamWS, upstreamCpGWS, downstreamCpGWS, 
                           upstreamkMerWS, downstreamkMerWS, kmer){
  # Input: metadata = the metadata information for each reads 
  #        rcFile = read class file.   
  #        neighbourhoodUpstreamWS = k, the number of peaks k bp upstream of the selected peaks.
  #        neighbourhoodDownstreamWS = k, the number of peaks k bp downstream of the selected peaks.
  #        labelWS = k, label = 1 if the distance of the selected peaks to the annotated tss is within k bp.
  #        upstreamWS = k, the number of sequence k bp upstream of the selected peak in the genome.
  #        downstreamWS = k, the number of sequence k bp downstream of the selected peak in the genome. 
  #        upstreamCpGWS = k <= upstreamWS, the number of sequence k bp upstream of the selected peak in the genome to be determined as CpG island region.
  #        downstreamCpGWS = k <= downstreamWS, the number of sequence k bp downstream of the selected peak in the genome to be determined as CpG island region.
  #        upstreamkMerWS = k <= upstreamWS, the number of sequence k bp upstream of the selected peak in the genome to be considered for computation of kmer distribution. 
  #        downstreamkMerWS = k <= downstreamWS, the number of sequence k bp downstream of the selected peak in the genome to be considered for computation of kmer distribution.
  #        kmer = a vector of integers, where each integer k indicates the kmer distribution to be computeed (eg. kmer = c(2,3) means 2-mer and 3-mer distribution)
  # Output: A feature matrix with labels for each peaks.
  
  # Read information
  ## Read class level 
  featureData.rc <- as.data.frame(metadata) %>%
    group_by(peaks, rc, readStrands) %>%
    summarize(peaks = peaks[1], rc = rc[1], GENEID = GENEID[1], chr = chr[1], readStrands = readStrands[1], 
              numReads.rc = n(), .groups = "drop") %>%
    ungroup() %>%
    mutate(ID = seq(n()), .before = "peaks") %>% 
    group_by(rc) %>%
    mutate(numReadsRelStrength.rc = numReads.rc / sum(numReads.rc), 
           neighbourhoodRatio.rc = map_dbl(peaks, ~(log(sum(numReads.rc[.x <= peaks & peaks <= .x + neighbourhoodDownstreamWS]) / sum(numReads.rc[.x - neighbourhoodUpstreamWS <= peaks & peaks <= .x])))), 
           neighbourhoodSum.rc = map_dbl(peaks, ~ sum(numReads.rc[.x - neighbourhoodUpstreamWS <= peaks & peaks <= .x + neighbourhoodDownstreamWS])), 
           neighbourhoodSumRelStrength.rc = neighbourhoodSum.rc / sum(neighbourhoodSum.rc)) %>% 
    ungroup() 
  
  ## Gene level 
  featureData.gene <- metadata %>%
    group_by(peaks, GENEID, readStrands) %>%
    summarize(peaks = peaks[1], GENEID = GENEID[1], readStrands = readStrands[1], 
              numReads.gene = n(), .groups = "drop") %>%
    ungroup() %>% 
    group_by(GENEID) %>% 
    mutate(numReadsRelStrength.gene = numReads.gene / sum(numReads.gene), 
           neighbourhoodRatio.gene = map_dbl(peaks, ~(log(sum(numReads.gene[.x <= peaks & peaks <= .x + neighbourhoodDownstreamWS]) / sum(numReads.gene[.x - neighbourhoodUpstreamWS <= peaks & peaks <= .x])))), 
           neighbourhoodSum.gene = map_dbl(peaks, ~ sum(numReads.gene[.x - neighbourhoodUpstreamWS <= peaks & peaks <= .x + neighbourhoodDownstreamWS])), 
           neighbourhoodSumRelStrength.gene = neighbourhoodSum.gene / sum(neighbourhoodSum.gene)) %>% 
    ungroup()
  
  ## Genome level 
  featureData.chr <- metadata %>%
    group_by(peaks, chr, readStrands) %>%
    summarize(peaks = peaks[1], chr = chr[1], readStrands = readStrands[1], 
              numReads.chr = n(),.groups = "drop") %>%
    ungroup() %>%
    group_by(chr) %>% 
    mutate(numReadsRelStrength.chr = numReads.chr / sum(numReads.chr)) %>% 
    ungroup()
  
  # Sequence information
  ## Get the upstream and downstream sequence information from genome file (fasta)
  featureData.sequence5Prime <- featureData.rc %>% filter(readStrands == TRUE)
  sequenceSet5Prime <- DNAStringSet(sapply(seq_along(featureData.sequence5Prime$peaks), ## **** this chunk of code takes a high computational time, hope can be improved ****
                                           function(x){return(as.character(subseq(fastaFile[grepl(paste0("GRCh38:", 
                                           featureData.sequence5Prime$chr[x], ":"), names(fastaFile))], 
                                           featureData.sequence5Prime$peaks[x] - upstreamWS, 
                                           featureData.sequence5Prime$peaks[x] + downstreamWS)))}))
  
  featureData.sequence3Prime <- featureData.rc %>% filter(readStrands == FALSE)
  sequenceSet3Prime <- DNAStringSet(sapply(seq_along(featureData.sequence3Prime$peaks), ## **** this chunk of code takes a high computational time, hope can be improved ****
                                           function(x){return(as.character(reverse(subseq(fastaFile[grepl(paste0("GRCh38:", 
                                           featureData.sequence3Prime$chr[x], ":"), names(fastaFile))], featureData.sequence3Prime$peaks[x] - upstreamWS, 
                                           featureData.sequence3Prime$peaks[x] + downstreamWS))))}))
  
  featureData.sequence <- featureData.rc 
  sequenceSet <- c(sequenceSet5Prime, sequenceSet3Prime)
  
  ## CpG island feature 
  ### upstream GC content feature
  
  upstreamStart <- upstreamWS - upstreamCpGWS + 1
  upstreamEnd <- upstreamWS
  
  downstreamStart <- upstreamWS + 2 
  downstreamEnd <- upstreamWS + downstreamCpGWS + 1
  
  newfeature <- as.vector(letterFrequency(subseq(sequenceSet, upstreamStart, upstreamEnd), "CG") / upstreamCpGWS) 
  featureData.sequence <- cbind(featureData.sequence, GCContent.up = newfeature)
  
  ### upstream observed to expected CpG ratio feature 
  observedCpG <- as.numeric(dinucleotideFrequency(subseq(sequenceSet, upstreamStart, upstreamEnd))[,"CG"])
  expectedCpG <- (1 + as.numeric(letterFrequency(subseq(sequenceSet, 
                 upstreamStart, upstreamEnd), "C")) * as.numeric(letterFrequency(subseq(sequenceSet, upstreamStart, upstreamEnd), "G")))/ upstreamCpGWS
  
  newfeature <- observedCpG / expectedCpG
  
  featureData.sequence <- cbind(featureData.sequence, OtECpGRatio.up = newfeature)
  
  ### downstream GC content feature 
  newfeature <- as.vector(letterFrequency(subseq(sequenceSet, downstreamStart, downstreamEnd), "CG") / downstreamCpGWS) 
  featureData.sequence <- cbind(featureData.sequence, GCContent.down = newfeature)
  
  ### observed to expected CpG ratio feature 
  observedCpG <- as.numeric(dinucleotideFrequency(subseq(sequenceSet, downstreamStart, downstreamEnd))[,"CG"])
  expectedCpG <- (1 + as.numeric(letterFrequency(subseq(sequenceSet, 
                 downstreamStart, downstreamEnd), "C")) * as.numeric(letterFrequency(subseq(sequenceSet, downstreamStart, downstreamEnd), "G")))/ downstreamCpGWS
  
  newfeature <- observedCpG / expectedCpG
  
  featureData.sequence <- cbind(featureData.sequence, OtECpGRatio.down = newfeature)
  
  ## Sequence Distribution Feature
  
  upstreamStart <- upstreamWS - upstreamkMerWS + 1
  upstreamEnd <- upstreamWS
  
  downstreamStart <- upstreamWS + 2 
  downstreamEnd <- upstreamWS + downstreamkMerWS + 1
  
  
  for (i in kmer){
      newfeature <- oligonucleotideFrequency(subseq(sequenceSet, upstreamStart, upstreamEnd), i)
      colnames(newfeature) <- paste0(colnames(newfeature), ".up")
      featureData.sequence <- cbind(featureData.sequence, newfeature)
      
      newfeature <- oligonucleotideFrequency(subseq(sequenceSet, downstreamStart, downstreamEnd), i)
      colnames(newfeature) <- paste0(colnames(newfeature), ".down")
      featureData.sequence <- cbind(featureData.sequence, newfeature)
  }

  featureData <- featureData.rc %>% 
    full_join(featureData.sequence, by = names(featureData.rc)) %>% 
    left_join(featureData.gene, by = c("peaks", "GENEID", "readStrands")) %>% 
    left_join(featureData.chr, by = c("peaks", "chr", "readStrands"))
  
  featureDataWithLabel <- featureData %>% 
    right_join(getAnnotatedStartSite(rcFile, featureData$rc), by = c("ID" = "queryHits", "rc")) %>%
    mutate(diff = abs(peaks - annotatedTss), .after = peaks) %>% 
    relocate(annotatedTss, .before = diff) %>% 
    group_by(ID) %>% 
    arrange(diff) %>% 
    filter(row_number() == 1) %>% 
    ungroup() %>% 
    mutate(label = ifelse(diff <= labelWS, 1, 0)) %>% 
    select(-c(ID, subjectHits, annotatedStarts, annotatedEnds, annotatedStrands)) 
  
  return(featureDataWithLabel)
}


# Main workflow (startendprediction)
# Train test split 
RCFile <- onlyAnnotatedRCWithMoreThanOneExon(RCFile)

set.seed(42)
trainingIndex <- sample.int(length(RCFile), size = 0.8 * length(RCFile))
testingIndex <- setdiff(seq(length(RCFile)), trainingIndex)

trainingRCFile <- RCFile[trainingIndex,]
testingRCFile <- RCFile[testingIndex]

# get the important metadata for training and testing using bam file and their corresponding rcFile 
trainingMetadata <- getMetadata(trainingRCFile)
testingMetadata <- getMetadata(testingRCFile)

# model training

## List of hyperparameters 

neighbourhoodUpstreamWS <- 15
neighbourhoodDownstreamWS <- 15
labelWS <- 10
upstreamWS <- 1500
downstreamWS <- 1500
upstreamCpGWS <- 500
downstreamCpGWS <- 500
upstreamkMerWS <- 500
downstreamkMerWS <- 500
kmer <- 2

trainingFeatureDataWithLabel <- getFeatureData(trainingMetadata, trainingRCFile, neighbourhoodUpstreamWS, neighbourhoodDownstreamWS, 
                                               labelWS, upstreamWS, downstreamWS, upstreamCpGWS, downstreamCpGWS, 
                                               upstreamkMerWS, downstreamkMerWS, kmer)

# trainingFeatureDataWithLabel <- trainingFeatureDataWithLabel %>% # Turn this off if read class labeling strat is used. This is labeling strat based on peaks
#     group_by(peaks) %>%
#     mutate(label = ifelse(any(label == 1), 1, 0)) %>%
#     ungroup()

validFeature <- names(which(apply(trainingFeatureDataWithLabel[,8:(length(trainingFeatureDataWithLabel) - 1)],2, sd) != 0)) # exclude features with no variance

scaledTrainingFeatureData <- scale(trainingFeatureDataWithLabel[,validFeature]) # normalise the data 
scaledTrainingFeatureDataWithLabel <- data.frame(trainingFeatureDataWithLabel[,1:7], scaledTrainingFeatureData, label = as.factor(trainingFeatureDataWithLabel$label))

set.seed(42)
model <- fitXGBoostModel(labels.train = as.numeric(as.character(scaledTrainingFeatureDataWithLabel$label)),  ## **** this chunk of code takes a high computational time, hope can be improved ****
                          data.train = as.matrix(scaledTrainingFeatureDataWithLabel[,-c(1:7,length(scaledTrainingFeatureDataWithLabel))]))

y.training <- predict(model, as.matrix(scaledTrainingFeatureDataWithLabel[,8:(length(scaledTrainingFeatureDataWithLabel) - 1)]))

training.final <- cbind(scaledTrainingFeatureDataWithLabel, y.training) %>% # final dataframe with labels and prediction score from the model for each read class 
    group_by(rc) %>% 
    arrange(desc(y.training)) %>% 
    filter(row_number() == 1) %>% # we choose the potential peaks (in the read class) that has the highest score as the peaks of the read class 
    ungroup() %>% 
    relocate(c(label,y.training), .before = peaks)

baselineTraining <- data.frame(rc = rownames(trainingRCFile), peaks.baseline = ifelse(rowData(trainingRCFile)$strand.rc == "+", # the peaks for each read class using Bambu's baseline
                                                                           as.vector(min(start(trainingRCFile))), 
                                                                           as.vector(max(end(trainingRCFile)))))

comparisonWithBaselineTraining <- training.final %>% # a dataframe to compare between baseline and model approach. 
    full_join(baselineTraining, by = "rc") %>% 
    mutate(diff.baseline = abs(annotatedTss - peaks.baseline)) %>% 
    select(rc, peaks, diff, peaks.baseline, diff.baseline, annotatedTss) %>% 
    filter(!is.na(peaks))

# model testing 
testingFeatureDataWithLabel <-  getFeatureData(testingMetadata, testingRCFile, neighbourhoodUpstreamWS, neighbourhoodDownstreamWS, 
                                               labelWS, upstreamWS, downstreamWS, upstreamCpGWS, downstreamCpGWS, 
                                               upstreamkMerWS, downstreamkMerWS, kmer)

# testingFeatureDataWithLabel <- testingFeatureDataWithLabel %>%  Turn this off if read class labeling strat is used. This is labeling strat based on peaks
#     group_by(peaks) %>% 
#     mutate(label = ifelse(any(label == 1), 1, 0)) %>% 
#     ungroup() 

scaledTestingFeatureData <- scale(testingFeatureDataWithLabel[,validFeature], # normalise the data 
                                  center = attr(scaledTrainingFeatureData, "scaled:center"), 
                                  scale = attr(scaledTrainingFeatureData, "scaled:scale"))
scaledTestingFeatureDataWithLabel <- data.frame(testingFeatureDataWithLabel[,1:7], 
                                                scaledTestingFeatureData, 
                                                label = as.factor(testingFeatureDataWithLabel$label))

y_pred <- predict(model, as.matrix(scaledTestingFeatureDataWithLabel[,8:(length(scaledTestingFeatureDataWithLabel) - 1)]))

final <- cbind(scaledTestingFeatureDataWithLabel, y_pred) %>% # final dataframe with labels and prediction score from the model for each read class 
    group_by(rc) %>% 
    arrange(desc(y_pred)) %>% 
    filter(row_number() == 1) %>% 
    ungroup() %>% 
    relocate(c(label,y_pred), .before = peaks)

baseline <- data.frame(rc = rownames(testingRCFile), peaks.baseline = ifelse(rowData(testingRCFile)$strand.rc == "+",  # the peaks for each read class using Bambu's baseline
                                                                           as.vector(min(start(testingRCFile))), 
                                                                           as.vector(max(end(testingRCFile)))))

comparisonWithBaseline <- final %>% # a dataframe to compare between baseline and model approach. 
    full_join(baseline, by = "rc") %>% 
    mutate(diff.baseline = abs(annotatedTss - peaks.baseline)) %>% 
    select(rc, peaks, diff, peaks.baseline, diff.baseline, annotatedTss, y_pred) %>% 
    filter(!is.na(peaks))




# Supplementary section 
## A function to get the top 3 local max peaks, it can be incorporated into the getFeatureData function 
## in the future if it shows similar or better performance to including all peaks. This resampling
## will definitely help in improving the computational efficiency. 

slidingBoxFreq <- function(windowSize, peakData, k){ 
  # Input: windowSize = the window size to consider the local maximum 
  #        peakData = a dataframe that contains the position of all peaks 
  #        k = the top number of (local max) peaks to be retained 
  # Output: a filtered peaks dataframe that include only top k local max peaks  
  
  peakData$slideBoxFrequency <- rep(0, nrow(peakData))
  for(i in seq(-windowSize, windowSize)){
    peakData <- peakData %>%
      mutate(interval = findInterval(peaks, seq(min(peaks) + i, max(peaks) + i, windowSize))) %>%
      group_by(rc, interval) %>% # read class level
      mutate(slideBoxFrequency = ifelse(numReads.rc == max(numReads.rc), slideBoxFrequency + 1, slideBoxFrequency)) %>% 
      ungroup() %>%
      group_by(rc) %>%
      filter(slideBoxFrequency == max(slideBoxFrequency)) %>%
      arrange(desc(numReads.rc)) %>%
      top_n(3, wt = numReads.rc) 
  }   
  
  return(peakData %>% select(-c("slideBoxFrequency", "interval")))
}

## Useful code to plot the ROC or PR curve 
readClassMinimumDiff <- testingMetadata %>% 
  mutate(ID = seq(n())) %>% 
  select(ID, peaks, rc) %>% 
  right_join(getAnnotatedStartSite(testingRCFile, testingMetadata$rc), by = c("ID" = "queryHits", "rc")) %>%
  mutate(diff = abs(peaks - annotatedTss), .after = peaks) %>% 
  relocate(annotatedTss, .before = diff) %>% 
  group_by(ID) %>% 
  arrange(diff) %>% 
  filter(row_number() == 1) %>% 
  ungroup() %>% 
  group_by(rc) %>% 
  summarize(min.diff = min(diff))

readClassMinimumDiffLessThanLabelWS <- (readClassMinimumDiff %>% filter(min.diff <= labelWS))$rc # read class with at least one read with distanceToAnnotation <= labelWS
totalReadClassLessThanLabelWS <- nrow(readClassMinimumDiffLessThanLabelWS)

final.readClassMinimumDiffLessThanLabelWS  <- cbind(scaledTestingFeatureDataWithLabel, y_pred = y_pred) %>% 
    filter(rc %in% readClassMinimumDiffLessThanLabelWS) %>% 
    relocate()

rocpr<- measureit(score = final.readClassMinimumDiffLessThanLabelWS$y_pred, 
                  class = final.readClassMinimumDiffLessThanLabelWS$label, 
                  measure = c("TPR", "FPR", "PREC", "REC"))

par(mfrow=c(1,2))
plot(rocpr$TPR~rocpr$FPR, col = 1, xlab = "1-Specificity (FPR)", ylab = "Sensitivity (TPR)", type = "l", lwd = 2, main = "ROC Curve")

plot(rocpr$PREC~rocpr$REC, col = 1, xlab = "Recall", ylab = "Precision", type = "l", lwd = 2, main = "PR Curve")
abline(h = table(final.readClassMinimumDiffLessThanLabelWS)[2], lty = 2, lwd = 2)

