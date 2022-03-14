context("bambu-processReads_utilityConstructReadClasses")


test_that("reads are correctly assigned to read classes",{
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
    expect_equal(rowData(se)$readIds,
    rowData(seReadClassUnstrandedExpected)$readIds)
})

test_that("reads are correctly assigned to transcripts",{

    seExpected <- readRDS(system.file("extdata", 
                                    "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", 
                                    package = "bambu"))
    annotations <- readRDS(system.file("extdata",
                                    "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds",
                                    package = "bambu"
    ))
    seReadClassUnstrandedExpected <- readRDS(system.file("extdata",
                                                        "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",
                                                        package = "bambu"
    ))
    emParameters <- setEmParameters(emParameters = NULL)
    bpParameters <- setBiocParallelParameters("sample1", readClass.file = NULL, 1, FALSE)
    isoreParameters <- setIsoreParameters(isoreParameters = NULL)

    countsSe <- bambu.quantify(seReadClassUnstrandedExpected,
                            annotations = annotations, isoreParameters = isoreParameters,
                            emParameters = emParameters, ncore = 1, verbose = FALSE)

    metadata(seReadClassUnstrandedExpected)$distTable = metadata(countsSe)$distTable
    readModelMap = generateReadModelMap(seReadClassUnstrandedExpected, trackReads=FALSE)
    expect_equal(metadata(seExpected)$readModelMap,
                readModelMap) 
})
