context("Isoform quantification")
library(bamboo)

test_that("generic function of isoform quantification of data.table is list of 2",{
  expect_type(bamboo(example_data[[1]], algo.control=list(ncore = 1)), "list")
  expect_output(bamboo(example_data[[1]], algo.control=list(ncore = 1)),"Finished EM estimation in")
  expect_type(bamboo(example_data[[2]], algo.control=list(ncore = 1)), "list")
  expect_output(bamboo(example_data[[2]], algo.control=list(ncore = 1)),"Finished EM estimation in")
  expect_type(bamboo(example_data[[3]], algo.control=list(ncore = 1)), "list")
  expect_output(bamboo(example_data[[3]], algo.control=list(ncore = 1)),"Finished EM estimation in")
  expect_type(bamboo(example_data[[4]], algo.control=list(ncore = 1)), "list")
  expect_output(bamboo(example_data[[4]], algo.control=list(ncore = 1)),"Finished EM estimation in")
})


test_that("isoform quantification of summarizedExperiment is summrizedExperiment",{
   expect_s4_class(bamboo(object = sysdata_se,algo.control=list(ncore = 1)), "summarizedExperiment")
   expect_output(bamboo(sysdata_se,algo.control=list(ncore = 1)),"Finished EM estimation in")
})

test_that("isoform quantification of bam file is summrizedExperiment",{
  samplefile <- system.file("extdata","Homo_sapiens.GRCh38.91.annotations-txdb.sqlite", package = "bamboo")
  txdb <- AnnotationDbi::loadDb(samplefile)
  test.bam <- system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo")
  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo")
  expect_output(bamboo(test.bam, fa.file = fa.file, txdb = txdb), "summarizedExperiment")
  expect_s4_class(bamboo(test.bam, fa.file = fa.file, txdb = txdb), "summarizedExperiment")
  expect_output(bamboo(test.bam, fa.file = fa.file, txdb = txdb),"Finished EM estimation in")
})
