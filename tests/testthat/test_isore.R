context("Isoform reconstruction")
library(bamboo)


test_that("isoform reconstruction runs without error on example data", {
  bamFile <- BamFile(system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo"),yieldSize=1000000)
  genomeFA <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo")
  txdbTablesList <- readRDS(url("http://s3.ap-southeast-1.amazonaws.com/ucsc-trackdata.store.genome.sg/chenying/bamboo_exampleDataset/HomoSapiens_ensembl91_Grch38_txdbTablesList.rds"))
  stranded=TRUE
  protocol='cDNA'
  prefix=''
  minimumReadSupport=2
  minimumTxFraction = 0.02
  testData <- isore(bamFile = bamFile,  txdbTablesList = txdbTablesList, genomeFA = genomeFA, stranded=stranded, protocol=protocol, prefix=prefix, minimumReadSupport=minimumReadSupport, minimumTxFraction=minimumTxFraction)
  expect_s4_class(testData, "SummarizedExperiment")
})

test_that("isoform reconstruction of bamFile is summarizedExperiment", {
  test.bam <- system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo")
  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo")
  example_txdbTablesList <- readRDS(url("http://s3.ap-southeast-1.amazonaws.com/ucsc-trackdata.store.genome.sg/chenying/bamboo_exampleDataset/HomoSapiens_ensembl91_Grch38_txdbTablesList.rds"))
  expect_output(isore(bamFile = test.bam,  txdbTablesList = example_txdbTablesList,genomeFA = fa.file), "### prepare annotations ###")
  expect_s4_class(isore(bamFile = test.bam,  txdbTablesList = example_txdbTablesList,genomeFA = fa.file), "SummarizedExperiment")
  expect_equal(isore(bamFile = test.bam,  txdbTablesList = example_txdbTablesList,genomeFA = fa.file), example_se)
})


test_that("isoform reconstruction of missing bam file is error", {

  test.bam <- system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo")
  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo")
  example_txdbTablesList <- readRDS(url("http://s3.ap-southeast-1.amazonaws.com/ucsc-trackdata.store.genome.sg/chenying/bamboo_exampleDataset/HomoSapiens_ensembl91_Grch38_txdbTablesList.rds"))
  expect_error(isore(bamFile = '',  txdbTablesList = example_txdbTablesList,genomeFA = fa.file),"Bam file is missing from arguments.")
})

test_that("isoform reconstruction of missing both txdb and txdbTablesList object is error", {

  test.bam <- system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo")
  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo")
  expect_error(isore(bamFile = test.bam,  txdb = NULL, txdbTablesList = NULL, genomeFA = fa.file),"txdb object is missing.")
})

test_that("isoform reconstruction of missing fa file is error", {
  test.bam <- system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo")
  example_txdbTablesList <- readRDS(url("http://s3.ap-southeast-1.amazonaws.com/ucsc-trackdata.store.genome.sg/chenying/bamboo_exampleDataset/HomoSapiens_ensembl91_Grch38_txdbTablesList.rds"))
  expect_error(isore(bamFile = test.bam,  txdbTablesList = example_txdbTablesList,genomeFA = ""))
  expect_error(isore(bamFile = test.bam,  txdbTablesList = example_txdbTablesList,genomeFA = NULL),"GenomeFA file is missing.")
})
