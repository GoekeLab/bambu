context("Prepare annotations")
library(bamboo)
library(RSQLite)

test_that("Prepare annotations of txdb is List",{
  samplefile <- "http://s3.ap-southeast-1.amazonaws.com/ucsc-trackdata.store.genome.sg/chenying/bamboo_exampleDataset/Homo_sapiens.GRCh38.91.annotations-txdb.sqlite"
  # download remote sqlite file as a temporary object
  dbfile <- tempfile(fileext=".sqlite")
  download.file(samplefile, dbfile)
  txdb <- AnnotationDbi::loadDb(dbfile)

  expectedFile <- readRDS('/mnt/ont/github/testdata/testthat_bamboo/annotationGranges_txdbGrch38_91.rds')

  gr <- prepareAnnotations(txdb)

  expect_equal(gr, expectedFile)
  expect_s4_class(gr, class = 'CompressedGRangesList')
  expect_named(mcols(gr), c("TXNAME", "GENEID","eqClass"))
})
