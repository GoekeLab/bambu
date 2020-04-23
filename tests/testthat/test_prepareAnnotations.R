context("Prepare annotations")
library(bamboo)
library(RSQLite)

test_that("prepareAnnotations of txdb object is a GRangesList",{
  samplefile <- "http://s3.ap-southeast-1.amazonaws.com/ucsc-trackdata.store.genome.sg/chenying/bamboo_exampleDataset/Homo_sapiens.GRCh38.91.annotations-txdb.sqlite"
  # download remote sqlite file as a temporary object
  dbfile <- tempfile(fileext=".sqlite")
  download.file(samplefile, dbfile)
  txdb <- AnnotationDbi::loadDb(dbfile)


  expectedGR <- readRDS(system.file("extdata", "annotationGrangesList_txdbGrch38_91.rds", package = "bamboo"))


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
  gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9.gtf", package = "bamboo")

  expectedGR <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9.rds", package = "bamboo"))

  gr <- prepareAnnotationsFromGTF(gtf.file)

  expect_equal(gr, expectedGR)
  expect_s4_class(gr, class = 'CompressedGRangesList')
  expect_named(mcols(gr), c("TXNAME", "GENEID","eqClass"))
})
