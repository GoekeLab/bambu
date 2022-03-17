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

    test.bam <- system.file("extdata", 
        "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam",
        package = "bambu")
    fa.file <- system.file("extdata", 
        "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", 
        package = "bambu")


    txdb <- AnnotationDbi::loadDb(system.file("extdata", 
        "Homo_sapiens.GRCh38.91.annotations-txdb_chr9_1_1000000.sqlite", 
        package = "bambu"))

    seExpected <- readRDS(system.file("extdata", 
        "seOutput_trackReads_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", 
        package = "bambu"))

    # test case 1: bambu with single bam file, only using annotations (default option)
    set.seed(1234)
    se <- bambu(reads = test.bam, annotations = txdb, genome = fa.file,
        opt.em = list(degradationBias = FALSE), trackReads = TRUE)
    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(se, seExpected)
})

test_that("returnDistTable returns the correct distTable",{
    test.bam <- system.file("extdata", 
        "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam",
        package = "bambu")
    fa.file <- system.file("extdata", 
        "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", 
        package = "bambu")


    txdb <- AnnotationDbi::loadDb(system.file("extdata", 
        "Homo_sapiens.GRCh38.91.annotations-txdb_chr9_1_1000000.sqlite", 
        package = "bambu"))

    seExpected <- readRDS(system.file("extdata", 
        "seOutput_distTable_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", 
        package = "bambu"))

    # test case 1: bambu with single bam file, only using annotations (default option)
    set.seed(1234)
    se <- bambu(reads = test.bam, annotations = txdb, genome = fa.file,
        opt.em = list(degradationBias = FALSE), discovery = FALSE, returnDistTable = TRUE)
    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(se, seExpected)
})

test_that("low Memory mode does not change results",{
    test.bam <- system.file("extdata", 
        "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam",
        package = "bambu")
    fa.file <- system.file("extdata", 
        "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", 
        package = "bambu")


    txdb <- AnnotationDbi::loadDb(system.file("extdata", 
        "Homo_sapiens.GRCh38.91.annotations-txdb_chr9_1_1000000.sqlite", 
        package = "bambu"))

    seExpected <- readRDS(system.file("extdata", 
        "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", 
        package = "bambu"))

    # test case 1: bambu with single bam file, only using annotations (default option)
    set.seed(1234)
    se <- bambu(reads = test.bam, annotations = txdb, genome = fa.file,
        opt.em = list(degradationBias = FALSE), discovery = FALSE, lowMemory = TRUE)
    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(se, seExpected)
})

test_that("Running bambu without annotations works",{
    test.bam <- system.file("extdata", 
        "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam",
        package = "bambu")
    fa.file <- system.file("extdata", 
        "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", 
        package = "bambu")

    seExpected <- readRDS(system.file("extdata", 
        "seOutput_denovo_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", 
        package = "bambu"))

    # test case 1: bambu with single bam file, only using annotations (default option)
    set.seed(1234)
    qbambu = purrr::quietly(bambu)
    se <- qbambu(reads = test.bam, annotations = NULL, genome = fa.file,
            opt.em = list(degradationBias = FALSE), lowMemory = TRUE, NDR = 1)
    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(se$result, seExpected)
})