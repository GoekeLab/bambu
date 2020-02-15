context("Isoform quantification")
library(bamboo)

test_that("generic function of isoform quantification of data.table is list of 2",{
  expect_type(bamboo(dt = example_data[[1]],  algo.control=list(ncore = 1)), "list")
  expect_output(bamboo(dt = example_data[[1]], algo.control=list(ncore = 1)),"Finished EM estimation in")
  expect_type(bamboo(dt = example_data[[2]], algo.control=list(ncore = 1)), "list")
  expect_output(bamboo(dt = example_data[[2]],  algo.control=list(ncore = 1)),"Finished EM estimation in")
  expect_type(bamboo(dt = example_data[[3]], algo.control=list(ncore = 1)), "list")
  expect_output(bamboo(dt = example_data[[3]], algo.control=list(ncore = 1)),"Finished EM estimation in")
  expect_type(bamboo(dt = example_data[[4]], algo.control=list(ncore = 1)), "list")
  expect_output(bamboo(dt = example_data[[4]], algo.control=list(ncore = 1)),"Finished EM estimation in")
})



test_that("isoform quantification of summarizedExperiment is summrizedExperiment",{
   expect_s4_class(bamboo(se = example_se,
                          txdb = txdb,
                          #txdbTablesList = example_txdbTablesList,
                          algo.control=list(ncore = 1)),
                          "SummarizedExperiment")
   expect_output(bamboo(se = example_se,txdbTablesList = example_txdbTablesList,algo.control=list(ncore = 1)),"Finished EM estimation in")
   expect_equal(bamboo(se = example_se,txdbTablesList = example_txdbTablesList,algo.control=list(ncore = 1)),example_quantOutput)
})

test_that("isoform quantification of bam file is summrizedExperiment",{

  test.bam <- system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo")
  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo")
  expect_s4_class(bamboo(bam.file = test.bam, fa.file = fa.file,
                         txdb = txdb,
                         #txdbTablesList = example_txdbTablesList,
                         algo.control=list(ncore = 1)), "SummarizedExperiment")
  expect_output(bamboo(bam.file = test.bam, fa.file = fa.file, txdbTablesList = example_txdbTablesList,algo.control=list(ncore = 1)),"Finished EM estimation in")
  expect_equal(bamboo(bam.file = test.bam, fa.file = fa.file, txdbTablesList = example_txdbTablesList,algo.control=list(ncore = 1)),example_quantOutput)
})
