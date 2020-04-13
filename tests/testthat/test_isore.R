context("Isoform reconstruction")
library(bamboo)

test_that("isore.constructReadClasses completes successfully", {

  readGrgList <- readRDS('/mnt/ont/github/testdata/testthat_bamboo/readGrgList_GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.rds')
  annotationGrangesList <- readRDS('/mnt/ont/github/testdata/testthat_bamboo/annotationGranges_txdbGrch38_91.rds')
  genomeSequence <- Rsamtools::FaFile(system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo"))

  seReadClassUnstrandedExpected <- readRDS(file='/mnt/ont/github/testdata/testthat_bamboo/seReadClassUnstranded.rds')
  seReadClassStrandedExpected <- readRDS(file='/mnt/ont/github/testdata/testthat_bamboo/seReadClassStranded.rds')

  seReadClassUnstranded <- isore.constructReadClasses(readGrgList=readGrgList,
                                                      runName='GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097_unstranded',
                                                      annotationGrangesList=annotationGrangesList,
                                                      genomeSequence=genomeSequence,
                                                      stranded=FALSE,
                                                      quickMode=FALSE,
                                                      verbose=FALSE)
  expect_equal(seReadClassUnstranded, seReadClassUnstrandedExpected)


  seReadClassStranded <- isore.constructReadClasses(readGrgList=readGrgList,
                                                    runName='GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097_stranded',
                                                    annotationGrangesList=annotationGrangesList,
                                                    genomeSequence=genomeSequence,
                                                    stranded=TRUE,
                                                    quickMode=FALSE,
                                                    verbose=FALSE)
  expect_equal(seReadClassStranded, seReadClassStrandedExpected)

  seReadClassFromBsgenome <- isore.constructReadClasses(readGrgList=readGrgList,
                                                        runName='GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097_stranded',
                                                        annotationGrangesList=annotationGrangesList,
                                                        genomeSequence="BSgenome.Hsapiens.NCBI.GRCh38",
                                                        stranded=TRUE,
                                                        quickMode=FALSE,
                                                        verbose=FALSE)

  expect_equal(seReadClassStranded, seReadClassFromBsgenome)
})

test_that("isore.combineTranscriptCandidates completes successfully", {

  seReadClass1 <- readRDS(file='/mnt/ont/github/testdata/testthat_bamboo/seReadClassUnstranded.rds')
  seReadClass2 <- readRDS(file='/mnt/ont/github/testdata/testthat_bamboo/seReadClassStranded.rds')

  seIsoReRefExpected <- readRDS(file='/mnt/ont/github/testdata/testthat_bamboo/seIsoReRef.rds')
  seIsoReCombinedExpected <- readRDS(file='/mnt/ont/github/testdata/testthat_bamboo/seIsoReCombined.rds')

  seIsoReRef <- isore.combineTranscriptCandidates(readClassSe=seReadClass1,
                                                  readClassSeRef = NULL,
                                                  stranded = FALSE,
                                                  verbose = FALSE)
  expect_equal(seIsoReRef, seIsoReRefExpected)

  seIsoReCombined <- isore.combineTranscriptCandidates(readClassSe=seReadClass2,
                                                       readClassSeRef = seIsoReRef,
                                                       stranded = FALSE,
                                                       verbose = FALSE)

  expect_equal(seIsoReCombined, seIsoReCombinedExpected)
  expect_named(assays(seIsoReCombined), c("counts", "start",  "end"))
  expect_named(rowData(seIsoReCombined), c("chr","start","end","strand","intronStarts", "intronEnds",  "confidenceType"))
})


test_that("isore.extendAnnotations completes successfully", {

  seIsoReCombined <- readRDS(file='/mnt/ont/github/testdata/testthat_bamboo/seIsoReCombined.rds')
  annotationGrangesList <- readRDS('/mnt/ont/github/testdata/testthat_bamboo/annotationGranges_txdbGrch38_91.rds')
  extendedAnnotationsExpected <- readRDS(file='/mnt/ont/github/testdata/testthat_bamboo/extendedAnnotations.rds')

  extendedAnnotations <- isore.extendAnnotations(se=seIsoReCombined,
                                                 annotationGrangesList=annotationGrangesList,
                                                 remove.subsetTx=TRUE,
                                                 min.readCount=2,
                                                 min.readFractionByGene=0.05,
                                                 min.sampleNumber=1,
                                                 min.exonDistance=35,
                                                 min.exonOverlap=10,
                                                 prefix='',
                                                 verbose=FALSE)
  expect_equal(extendedAnnotations, extendedAnnotationsExpected)
})


test_that("isore.estimateDistanceToAnnotations completes successfully", {

  seReadClass1 <- readRDS(file='/mnt/ont/github/testdata/testthat_bamboo/seReadClassUnstranded.rds')
  extendedAnnotations <- readRDS(file='/mnt/ont/github/testdata/testthat_bamboo/extendedAnnotations.rds')
  seWithDistExpected <- readRDS(file='/mnt/ont/github/testdata/testthat_bamboo/seWithDist.rds')
  seWithDist <- isore.estimateDistanceToAnnotations(seReadClass=seReadClass1,
                                      annotationGrangesList=extendedAnnotations,
                                      min.exonDistance = 35)

  expect_equal(seWithDist, seWithDistExpected)
})
