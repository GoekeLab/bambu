context("Isoform reconstruction")


test_that("isore.constructReadClasses completes successfully", {
    readGrgList <- readRDS(system.file("extdata",
                                       "readGrgList_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",
                                       package = "bambu"
    ))
    annotations <- readRDS(system.file("extdata",
                                       "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds",
                                       package = "bambu"
    ))
    genomeSequence <- system.file("extdata",
                                  "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa",
                                  package = "bambu"
    )
    
    
    seReadClassUnstrandedExpected <- readRDS(system.file("extdata",
                                                         "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",
                                                         package = "bambu"
    ))
    seReadClassStrandedExpected <- readRDS(system.file("extdata",
                                                       "seReadClassStranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",
                                                       package = "bambu"
    ))
    
    
    seqlevelCheckReadsAnnotation(readGrgList, annotations)
    genomeSequence <- checkInputSequence(genomeSequence)
    #check seqlevels for consistency, drop ranges not present in genomeSequence
    refSeqLevels <- seqlevels(genomeSequence)
    if (!all(seqlevels(readGrgList) %in% refSeqLevels)) {
        message("not all chromosomes from reads present in reference genome 
            sequence, reads without reference chromosome sequence are dropped")
        refSeqLevels <- intersect(refSeqLevels, seqlevels(readGrgList))
        readGrgList <- keepSeqlevels(readGrgList, value =  refSeqLevels,
                                     pruning.mode = "coarse")
        # reassign Ids after seqlevels are dropped
        mcols(readGrgList)$id <- seq_along(readGrgList) 
    }
    if (!all(seqlevels(annotations) %in% refSeqLevels)) {
        message("not all chromosomes from annotations present in reference genome 
    sequence, annotations without reference chrosomomse sequence are dropped")
        annotations <- keepSeqlevels(annotations, value = refSeqLevels,
                                     pruning.mode = "coarse")
    }
    # create error and strand corrected junction tables
    unlisted_junctions <- unlistIntrons(readGrgList, use.ids = TRUE)
    uniqueJunctions <- isore.constructJunctionTables(unlisted_junctions, 
        annotations,genomeSequence, stranded = FALSE, verbose = FALSE)
    # create SE object with reconstructed readClasses
    se <- isore.constructReadClasses(readGrgList, unlisted_junctions, 
        uniqueJunctions, 
        runName =  "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000",
        annotations,  stranded = FALSE, verbose = FALSE)
    GenomeInfoDb::seqlevels(se) <- refSeqLevels
    seReadClassUnstranded <- scoreReadClasses(se,genomeSequence, annotations, 
                                              defaultModels = defaultModels, fit = TRUE,
                                              min.readCount = 2, verbose = FALSE)
    ## in case of testing on Mac
    names(seReadClassUnstranded@rowRanges@elementMetadata@listData$intronStarts) <-
        names(seReadClassUnstrandedExpected@rowRanges@elementMetadata@listData$intronStarts) <- NULL
    names(seReadClassUnstranded@rowRanges@elementMetadata@listData$intronEnds) <- 
        names(seReadClassUnstrandedExpected@rowRanges@elementMetadata@listData$intronEnds) <- NULL
    expect_equal(seReadClassUnstranded, seReadClassUnstrandedExpected)
    
    se <- isore.constructReadClasses(readGrgList = readGrgList,
                                     unlisted_junctions, 
                                     uniqueJunctions,
                                     runName = "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000_Stranded",
                                     annotations,
                                     stranded = TRUE,
                                     verbose = FALSE)
    GenomeInfoDb::seqlevels(se) <- refSeqLevels
    seReadClassStranded <- scoreReadClasses(se,genomeSequence, annotations, 
                                            defaultModels = defaultModels, fit = TRUE,
                                            min.readCount = 2, verbose = FALSE)
    names(seReadClassStranded@rowRanges@elementMetadata@listData$intronStarts) <- 
        names(seReadClassStrandedExpected@rowRanges@elementMetadata@listData$intronStarts) <- NULL
    names(seReadClassStranded@rowRanges@elementMetadata@listData$intronEnds) <- 
        names(seReadClassStrandedExpected@rowRanges@elementMetadata@listData$intronEnds) <- NULL
    expect_equal(seReadClassStranded, seReadClassStrandedExpected)
})

test_that("isore.combineTranscriptCandidates completes successfully", {
    seReadClass1 <- readRDS(system.file("extdata",
                                        "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",
                                        package = "bambu"
    ))
    seReadClass2 <- readRDS(system.file("extdata",
                                        "seReadClassStranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",
                                        package = "bambu"
    ))
    
    seIsoReRefExpected <- readRDS(system.file("extdata",
                                              "seIsoReRef_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",
                                              package = "bambu"
    ))
    seIsoReCombinedExpected <- readRDS(system.file("extdata",
                                                   "seIsoReCombined_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",
                                                   package = "bambu"
    ))
    rcFileList <- system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu")
    
    bpParameters <- setBiocParallelParameters(reads = NULL, readClass.file = rcFileList,
                                              ncore = 1, verbose = FALSE)
    seIsoReRef <- isore.combineTranscriptCandidates(readClassList = list(seReadClass1),
                                                    stranded = FALSE,min.readCount = 2,
                                                    min.txScore.multiExon = 0,
                                                    min.txScore.singleExon = 1,
                                                    min.readFractionByGene = 0.05,
                                                    verbose = FALSE, bpParameters = bpParameters)
    expect_equal(seIsoReRef, seIsoReRefExpected)
    
    rcFileList <- c(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"),
                    system.file("extdata", "seReadClassStranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
    bpParameters <- setBiocParallelParameters(reads = NULL, readClass.file = rcFileList,
                                              ncore = 1, verbose = FALSE)
    seIsoReCombined <- isore.combineTranscriptCandidates(readClassList = list(seReadClass1,seReadClass2),
                                                         stranded = FALSE,
                                                         min.readCount = 2,
                                                         min.txScore.multiExon = 0,
                                                         min.txScore.singleExon = 1,
                                                         min.readFractionByGene = 0.05,
                                                         verbose = FALSE, bpParameters = bpParameters)
    
    expect_equal(seIsoReCombined, seIsoReCombinedExpected)
    expect_named(seIsoReCombined,
                 c('intronStarts', 'intronEnds', 'chr', 'strand', 'maxTxScore', 
                   'NSampleReadCount', 'NSampleReadProp', 'NSampleTxScore', 'start', 'end', 'readCount', 'confidenceType') 
    )
})


test_that("isore.extendAnnotations completes successfully", {
    seIsoReCombined <- readRDS(system.file("extdata",
                                           "seIsoReCombined_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",
                                           package = "bambu"
    ))
    gr <- readRDS(system.file("extdata",
                              "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds",
                              package = "bambu"
    ))
    
    extendedAnnotationsExpected <- readRDS(system.file("extdata",
                                                       "extendedAnnotationGranges_txdbGrch38_91_chr9_1_1000000.rds",
                                                       package = "bambu"
    ))
    
    extendedAnnotations <- isore.extendAnnotations(
        combinedTranscripts=seIsoReCombined,
        annotationGrangesList=gr,
        remove.subsetTx = TRUE, min.sampleNumber = 1, NDR = 0.1, 
        min.exonDistance = 35, min.exonOverlap = 10,
        min.primarySecondaryDist = 5, min.primarySecondaryDistStartEnd = 5, 
        prefix='',
        verbose=FALSE
    )
    expect_equal(extendedAnnotations, extendedAnnotationsExpected)
})


test_that("isore.estimateDistanceToAnnotations completes successfully", {
    seReadClass1 <- readRDS(system.file("extdata",
                                        "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",
                                        package = "bambu"
    ))
    extendedAnnotations <- readRDS(system.file("extdata",
                                               "extendedAnnotationGranges_txdbGrch38_91_chr9_1_1000000.rds",
                                               package = "bambu"
    ))
    
    seWithDist <- isore.estimateDistanceToAnnotations(
        seReadClass = seReadClass1,
        annotationGrangesList = extendedAnnotations,
        min.exonDistance = 35
    )
    names(seWithDist@metadata$distTable$readCount) <- NULL
    names(seWithDistExpected@metadata$distTable$readCount) <- NULL
    expect_equal(seWithDist, seWithDistExpected)
})
