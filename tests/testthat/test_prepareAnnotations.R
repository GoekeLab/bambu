context("Prepare annotations")
library(bamboo)
library(RSQLite)

test_that("Prepare annotations of txdb is List",{
  samplefile <- "http://s3.ap-southeast-1.amazonaws.com/ucsc-trackdata.store.genome.sg/chenying/bamboo_exampleDataset/Homo_sapiens.GRCh38.91.annotations-txdb.sqlite"
  # download remote sqlite file as a temporary object
  dbfile <- tempfile(fileext=".sqlite")
  download.file(samplefile, dbfile)
  txdb <- AnnotationDbi::loadDb(dbfile)
  expect_output(prepareAnnotations(txdb), "test that rows are ordered by tx name")
  expect_type(prepareAnnotations(txdb),"list")
  expect_named(prepareAnnotations(txdb), c("intronsByTxEns" ,"exonsByTx","txIdToGeneIdTable","unlisted_introns" ))
  example_txdbTablesList <- readRDS(url("http://s3.ap-southeast-1.amazonaws.com/ucsc-trackdata.store.genome.sg/chenying/bamboo_exampleDataset/HomoSapiens_ensembl91_Grch38_txdbTablesList.rds"))
  expect_equal(prepareAnnotations(txdb), example_txdbTablesList)
})
