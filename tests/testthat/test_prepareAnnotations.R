context("Prepare annotations")
library(bamboo)


test_that("Prepare annotations of txdb is List",{
  samplefile <- system.file("extdata","Homo_sapiens.GRCh38.91.annotations-txdb.sqlite", package = "bamboo")
  txdb <- AnnotationDbi::loadDb(samplefile)
  expect_output(prepareAnnotations(txdb), "test that rows are ordered by tx name")
  expect_type(prepareAnnotations(txdb),"list")
  expect_named(prepareAnnotations(txdb), c("intronsByTxEns" ,"exonsByTx","txIdToGeneIdTable","unlisted_introns" ))
  expect_equal(prepareAnnotations(txdb), example_txdbTablesList)
})
