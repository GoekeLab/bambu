context("Prepare annotations")
library(bambu)
library(RSQLite)
library(testthat)

test_that("prepareAnnotations of txdb object is a GRangesList",{
  txdb <- loadDb(system.file("extdata", "Homo_sapiens.GRCh38.91.annotations-txdb_chr9_108865774_109014097.sqlite", package = "bambu"))

  expectedGR <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_108865774_109014097.rds", package = "bambu"))

  gr <- prepareAnnotations(txdb)

  expect_equal(gr, expectedGR)
  expect_s4_class(gr, class = 'CompressedGRangesList')
  expect_named(mcols(gr), c("TXNAME", "GENEID","eqClass"))
})


test_that("prepareAnnotations of genome library is a GRangesList",{
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  gr <- prepareAnnotations(TxDb.Hsapiens.UCSC.hg38.knownGene)

  expect_s4_class(gr, class = 'CompressedGRangesList')
  expect_named(mcols(gr), c("TXNAME", "GENEID","eqClass"))
})




test_that("prepareAnnotationsFromGTF is GRangesList",{
  gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_108865774_109014097.gtf", package = "bambu")

  expectedGR <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_108865774_109014097.rds", package = "bambu"))

  gr <- prepareAnnotationsFromGTF(gtf.file)

  expect_equal(gr, expectedGR)
  expect_s4_class(gr, class = 'CompressedGRangesList')
  expect_named(mcols(gr), c("TXNAME", "GENEID","eqClass"))
})


