context("setNDR")

test_that("txRange generates a gene and transcript score",{
    readClasses <- readRDS(system.file("extdata","readClassesUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds",package = "bambu"))
    annotations <- readRDS(system.file("extdata","annotationGranges_txdbGrch38_91_chr9_1_1000000.rds",package = "bambu"))
    genomeSequence <- system.file("extdata","Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa",package = "bambu")
    genomeSequence <- checkInputSequence(genomeSequence)
    
    se = scoreReadClasses(readClasses, genomeSequence, annotations, defaultModels)
  
    seExpected = readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  
    expect_is(rowData(se)$txScore, class = 'numeric')
    expect_equal(rowData(se)$txScore, rowData(seExpected)$txScore)
})

test_that("Correctly moves transcripts with NDR above the threshold into lowCofidenceTranscripts",{
    extendedAnnotationsExpected <- readRDS(system.file("extdata", "extendedAnnotationGranges_txdbGrch38_91_chr9_1_1000000.rds",package = "bambu"))

    expect_equal(length(extendedAnnotationsExpected), 106)
    expect_equal(length(metadata(extendedAnnotationsExpected)$lowConfidenceTranscripts),1)

    extendedAnnotationsExpected_0.1 = setNDR(extendedAnnotationsExpected, 0.1)
    expect_equal(length(extendedAnnotationsExpected_0.1), 105)
    expect_equal(length(metadata(extendedAnnotationsExpected_0.1)$lowConfidenceTranscripts),2)
})

test_that("Correctly moves transcripts with NDR below the threshold into the extendedAnnotations",{
    extendedAnnotationsExpected <- readRDS(system.file("extdata", "extendedAnnotationGranges_txdbGrch38_91_chr9_1_1000000.rds",package = "bambu"))

    expect_equal(length(extendedAnnotationsExpected), 106)
    expect_equal(length(metadata(extendedAnnotationsExpected)$lowConfidenceTranscripts),1)

    extendedAnnotationsExpected_1 = setNDR(extendedAnnotationsExpected, 1)
    expect_equal(length(extendedAnnotationsExpected_1), 107)
    expect_equal(length(metadata(extendedAnnotationsExpected_1)$lowConfidenceTranscripts),0)
})

test_that("Reference annotations are not moved unless includeRef is used", {
    extendedAnnotationsExpected <- readRDS(system.file("extdata", "extendedAnnotationGranges_txdbGrch38_91_chr9_1_1000000.rds",package = "bambu"))

    expect_equal(length(extendedAnnotationsExpected), 106)
    expect_equal(length(metadata(extendedAnnotationsExpected)$lowConfidenceTranscripts),1)

    extendedAnnotationsExpected_0.5 = setNDR(extendedAnnotationsExpected, 0.5, includeRef = TRUE)
    expect_equal(length(extendedAnnotationsExpected_0.5), 104)
    expect_equal(length(metadata(extendedAnnotationsExpected_0.5)$lowConfidenceTranscripts),3)
})

test_that("setNDR only effects transcripts with the prefix", {
    extendedAnnotationsExpected <- readRDS(system.file("extdata", "extendedAnnotationGranges_txdbGrch38_91_chr9_1_1000000.rds",package = "bambu"))

    expect_equal(length(extendedAnnotationsExpected), 106)
    expect_equal(length(metadata(extendedAnnotationsExpected)$lowConfidenceTranscripts),1)

    extendedAnnotationsExpected_1 = setNDR(extendedAnnotationsExpected, 1, prefix = "Stringtie2")
    expect_equal(length(extendedAnnotationsExpected_1), 106)
    expect_equal(length(metadata(extendedAnnotationsExpected_1)$lowConfidenceTranscripts),1)
})

test_that("setNDR saves the used NDR threshold correctly", {
    extendedAnnotationsExpected <- readRDS(system.file("extdata", "extendedAnnotationGranges_txdbGrch38_91_chr9_1_1000000.rds",package = "bambu"))
    metadata(extendedAnnotationsExpected)$NDR = 0.7
    extendedAnnotationsExpected_1 = setNDR(extendedAnnotationsExpected, 1)
    metadata(extendedAnnotationsExpected_1)$NDR = 1
})

test_that("setNDR works when no lowConfidenceTranscripts are present", {
    seIsoReCombined <- readRDS(system.file("extdata", "seIsoReCombined_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
    gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))

    extendedAnnotations <- isore.extendAnnotations(combinedTranscripts=seIsoReCombined,
                                                annotationGrangesList=gr,
                                                remove.subsetTx = TRUE, min.sampleNumber = 1, NDR = 1, 
                                                min.exonDistance = 35, min.exonOverlap = 10,
                                                min.primarySecondaryDist = 5, min.primarySecondaryDistStartEnd = 5, 
                                                prefix='Bambu', verbose=FALSE, defaultModels = defaultModels)

    extendedAnnotationsExpected_0.7 = setNDR(extendedAnnotationsExpected, 0.7)
    expect_equal(length(extendedAnnotationsExpected), 106)
    expect_equal(length(metadata(extendedAnnotationsExpected)$lowConfidenceTranscripts),1)
})

test_that("setNDR recommends the correct NDR", {
    extendedAnnotationsExpected <- readRDS(system.file("extdata", "extendedAnnotationGranges_txdbGrch38_91_chr9_1_1000000.rds",package = "bambu"))
    
    extendedAnnotationsExpected_rec = setNDR(extendedAnnotationsExpected)
    expect_equal(unname(metadata(extendedAnnotationsExpected_rec)$NDR), -0.013)
})

test_that("setNDR handles annotations with no NDR", {
    annotations <- readRDS(system.file("extdata","annotationGranges_txdbGrch38_91_chr9_1_1000000.rds",package = "bambu"))
    expect_warning(setNDR(annotations, 0.5), "Annotations were not extended by Bambu")
})