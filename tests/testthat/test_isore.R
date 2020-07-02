context("Isoform reconstruction")


test_that("isore.constructReadClasses completes successfully", {

  readGrgList <- readRDS(system.file("extdata", "readGrgList_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))
  genomeSequence <- Rsamtools::FaFile(system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", package = "bambu"))


  seReadClassUnstrandedExpected <- readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  seReadClassStrandedExpected <- readRDS(system.file("extdata", "seReadClassStranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  seReadClassFromBsgenomeExpected <- readRDS(system.file("extdata", "seReadClassBsgenomeStranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))

  seReadClassUnstranded <- isore.constructReadClasses(readGrgList = readGrgList,
                                                      runName ='SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000_unstranded',
                                                      annotationGrangesList = gr,
                                                      genomeSequence = genomeSequence,
                                                      stranded = FALSE,
                                                      ncore = 1,
                                                      verbose = FALSE)
  expect_equal(seReadClassUnstranded, seReadClassUnstrandedExpected)


  seReadClassStranded <- isore.constructReadClasses(readGrgList = readGrgList,
                                                    runName = 'SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000_stranded',
                                                    annotationGrangesList = gr,
                                                    genomeSequence = genomeSequence,
                                                    stranded=TRUE,
                                                    ncore = 1,
                                                    verbose = FALSE)
  expect_equal(seReadClassStranded, seReadClassStrandedExpected)

  seReadClassFromBsgenome <- isore.constructReadClasses(readGrgList = readGrgList,
                                                        runName = 'SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000_stranded',
                                                        annotationGrangesList = gr,
                                                        genomeSequence = "BSgenome.Hsapiens.NCBI.GRCh38",
                                                        stranded = TRUE,
                                                        ncore = 1,
                                                        verbose = FALSE)

  expect_equal(seReadClassFromBsgenome,seReadClassFromBsgenomeExpected)

})

test_that("isore.combineTranscriptCandidates completes successfully", {

  seReadClass1 <- readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  seReadClass2 <- readRDS(system.file("extdata", "seReadClassStranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))

  seIsoReRefExpected <- readRDS(system.file("extdata", "seIsoReRef_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  seIsoReCombinedExpected <- readRDS(system.file("extdata", "seIsoReCombined_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))

  seIsoReRef <- isore.combineTranscriptCandidates(readClassSe = seReadClass1,
                                                  readClassSeRef = NULL,
                                                  stranded = FALSE,
                                                  verbose = FALSE)
  expect_equal(seIsoReRef, seIsoReRefExpected)

  seIsoReCombined <- isore.combineTranscriptCandidates(readClassSe = seReadClass2,
                                                       readClassSeRef = seIsoReRef,
                                                       stranded = FALSE,
                                                       verbose = FALSE)

  expect_equal(seIsoReCombined, seIsoReCombinedExpected)
  expect_named(assays(seIsoReCombined), c("counts", "start",  "end"))
  expect_named(rowData(seIsoReCombined), c("chr","start","end","strand","intronStarts", "intronEnds",  "confidenceType"))
})


test_that("isore.extendAnnotations completes successfully", {

  seIsoReCombined <- readRDS(system.file("extdata", "seIsoReCombined_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))

  extendedAnnotationsExpected <- readRDS(system.file("extdata", "extendedAnnotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))

  extendedAnnotations <- isore.extendAnnotations(se = seIsoReCombined,
                                                 annotationGrangesList = gr,
                                                 remove.subsetTx = TRUE,
                                                 min.readCount = 2,
                                                 min.readFractionByGene = 0.05,
                                                 min.sampleNumber = 1,
                                                 min.exonDistance = 35,
                                                 min.exonOverlap = 10,
                                                 prefix = '',
                                                 verbose=FALSE)
  expect_equal(extendedAnnotations, extendedAnnotationsExpected)
})


test_that("isore.estimateDistanceToAnnotations completes successfully", {

  seReadClass1 <- readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  extendedAnnotations <- readRDS(system.file("extdata", "extendedAnnotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))

  seWithDist <- isore.estimateDistanceToAnnotations(seReadClass = seReadClass1,
                                      annotationGrangesList = extendedAnnotations,
                                      min.exonDistance = 35)

  expect_equal(seWithDist, seWithDistExpected)
})
