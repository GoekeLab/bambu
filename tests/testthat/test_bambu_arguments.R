context("bambu-processReads_utilityConstructReadClasses")


test_that("reads are correctly assigned to read classes",{
    test.bam <- system.file("extdata", "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam", package = "bambu")
    annotations <- readRDS(system.file("extdata","annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))
    genomeSequence <- system.file("extdata","Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", package = "bambu")
    
    seReadClassUnstrandedExpected <- readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
    
    se <- bambu(reads = test.bam, annotations = annotations, genome = genomeSequence, discovery = FALSE, quant = FALSE)[[1]]
    
    expect_equal(rowData(se)$readIds,
                 rowData(seReadClassUnstrandedExpected)$readIds)
})

test_that("reads are correctly assigned to transcripts",{
    
    test.bam <- system.file("extdata", "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam", package = "bambu")
    fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", package = "bambu")
    annotations <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))
    
    seExpected <- readRDS(system.file("extdata", 
                                      "seOutput_trackReads_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", 
                                      package = "bambu"))
    
    # test case 1: bambu with single bam file, only using annotations (default option)
    se <- bambu(reads = test.bam, annotations = annotations, genome = fa.file, trackReads = TRUE)
    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(se, seExpected)
})

test_that("returnDistTable returns the correct distTable",{
    test.bam <- system.file("extdata", "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam", package = "bambu")
    fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", package = "bambu")
    annotations <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))
    
    seExpected <- readRDS(system.file("extdata", "seOutput_distTable_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
    
    se <- bambu(reads = test.bam, annotations = annotations, genome = fa.file, discovery = FALSE, returnDistTable = TRUE)
    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(se, seExpected)
})

test_that("low Memory mode does not change results",{
    test.bam <- system.file("extdata", "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam", package = "bambu")
    fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", package = "bambu")
    annotations <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))
    
    seExpected <- readRDS(system.file("extdata", "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
    
    se <- bambu(reads = test.bam, annotations = annotations, genome = fa.file, lowMemory = TRUE)
    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(se, seExpected)
})

test_that("Running bambu without annotations works",{
    test.bam <- system.file("extdata", "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam", package = "bambu")
    fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", package = "bambu")
    
    seExpected <- readRDS(system.file("extdata", "seOutput_denovo_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
    
    qbambu = purrr::quietly(bambu)
    se <- qbambu(reads = test.bam, annotations = NULL, genome = fa.file, NDR = 1)
    expect_s4_class(se$result, "SummarizedExperiment")
    expect_equal(se$result, seExpected)
})