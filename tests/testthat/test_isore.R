context("Isoform reconstruction")
library(bamboo)


test_that("isoform reconstruction runs without error on example data", {
  bamFile <- BamFile(system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo"),yieldSize=1000000)
  genomeFA <- FaFile(system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo"))
  txdbTablesList <- readRDS(url("http://s3.ap-southeast-1.amazonaws.com/ucsc-trackdata.store.genome.sg/chenying/bamboo_exampleDataset/HomoSapiens_ensembl91_Grch38_txdbTablesList.rds"))
  stranded=FALSE
  protocol='cDNA'
  prefix=''
  minimumReadSupport=2
  minimumTxFraction = 0.02
  yieldSize=1000000
  testData <- isore(bamFile = bamFile,  txdbTablesList = txdbTablesList, genomeFA = genomeFA, stranded=stranded, protocol=protocol, prefix=prefix, minimumReadSupport=minimumReadSupport, minimumTxFraction=minimumTxFraction, yieldSize = yieldSize)
  expect_s4_class(testData, "SummarizedExperiment")
})

test_that("combining isoform reconstruction runs without error on example data", {
  bamFile <- BamFile(system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo"),yieldSize=1000000)

  bfListSmall <- BamFileList(bamFile, bamFile)

  bamFile <- BamFile('/mnt/ont/s3.ontdata.store.genome.sg/Nanopore/03_Mapping/Grch38/minimap2-2.17-cDNA/GIS_Hct116_directcDNA_Rep2-Run2/GIS_Hct116_directcDNA_Rep2-Run2.bam',yieldSize=10000000)
  bamFile3 <- BamFile('/mnt/ont/s3.ontdata.store.genome.sg/Nanopore/03_Mapping/Grch38/minimap2-2.17-directRNA/GIS_Hct116_directRNA_3/GIS_Hct116_directRNA_3.bam',yieldSize=1000000)
  bamFile4 <- BamFile('/mnt/ont/s3.ontdata.store.genome.sg/Nanopore/03_Mapping/Grch38/minimap2-2.17-directRNA/GIS_Hct116_directRNA_4/GIS_Hct116_directRNA_4.bam',yieldSize=1000000)
  genomeFA = FaFile('/mnt/ont/s3.ontdata.store.genome.sg/annotations/Grch38/ensembl-91/Homo_sapiens.GRCh38.dna_sm.primary_assembly_wtChrIS.fa')

  txdbTablesList <- readRDS(url("http://s3.ap-southeast-1.amazonaws.com/ucsc-trackdata.store.genome.sg/chenying/bamboo_exampleDataset/HomoSapiens_ensembl91_Grch38_txdbTablesList.rds"))
  stranded=FALSE
  protocol='directRNA'
  prefix=''
  minimumReadSupport=2
  minimumTxFraction = 0.02
  yieldSize=1000000
  testData <- isore(bamFile = bamFile,  txdbTablesList = txdbTablesList, genomeFA = genomeFA, stranded=stranded, protocol=protocol, prefix=prefix, minimumReadSupport=minimumReadSupport, minimumTxFraction=minimumTxFraction, yieldSize = yieldSize)
  expect_s4_class(testData, "SummarizedExperiment")

})

test_that("constructSplicedReadClassTables runs without error", {

  uniqueJunctions <- readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directcDNA_Rep2-Run2_uniqueJunctions.rds')
  unlisted_junctions <- readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directcDNA_Rep2-Run2_unlisted_junctions.rds')
  readData <- readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directcDNA_Rep2-Run2_readData.rds')
  readClassListSpliced <- constructSplicedReadClassTables(uniqueJunctions, unlisted_junctions, readData$readGrglist, readData$readNames)
  readClassListSplicedOriginal <- readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directcDNA_Rep2-Run2_readClassListSpliced.rds')
  expect_equal(readClassListSpliced,readClassListSplicedOriginal)

  readClassTableTemp_2<-readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directcDNA_Rep2-Run2_readClassTableIntermediate.rds')



  uniqueJunctions <- readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directRNA_3_uniqueJunctions.rds')
  unlisted_junctions <- readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directRNA_3_unlisted_junctions.rds')
  readData <- readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directRNA_3_readData.rds')
  readClassListSpliced <- constructSplicedReadClassTables(uniqueJunctions, unlisted_junctions, readData$readGrglist, readData$readNames)
  readClassListSplicedOriginal <- readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directRNA_3_readClassListSpliced.rds')
  expect_equal(readClassListSpliced,readClassListSplicedOriginal)

  readClassTableTemp_3<- readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directRNA_3_readClassTableIntermediate.rds')


  uniqueJunctions <- readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directRNA_4_uniqueJunctions.rds')
  unlisted_junctions <- readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directRNA_4_unlisted_junctions.rds')
  readData <- readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directRNA_4_readData.rds')
  readClassListSpliced <- constructSplicedReadClassTables(uniqueJunctions, unlisted_junctions, readData$readGrglist, readData$readNames)
  readClassListSplicedOriginal <- readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directRNA_4_readClassListSpliced.rds')
  expect_equal(readClassListSpliced,readClassListSplicedOriginal)

   readClassTableTemp_4<-readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directRNA_4_readClassTableIntermediate.rds')


  ## with single exon read classes
  exonsByReadClass_3 <- readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directRNA_3_exonsByReadClass.rds')
  readClassTable_3 <- readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directRNA_3_readClassTable.rds')
  readTable_3 <- readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directRNA_3_readTable.rds')

  exonsByReadClass_2 <- readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directcDNA_Rep2-Run2_exonsByReadClass.rds')
  readClassTable_2 <- readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directcDNA_Rep2-Run2_readClassTable.rds')
  readTable_2 <- readRDS(file='/mnt/ont/github/testdata/GIS_Hct116_directcDNA_Rep2-Run2_readTable.rds')


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
