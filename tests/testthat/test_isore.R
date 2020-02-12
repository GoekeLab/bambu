context("Isoform reconstruction")
library(bamboo)

test_that("Prepare annotations of txdb is List",{
  samplefile <- system.file("extdata","Homo_sapiens.GRCh38.91.annotations-txdb.sqlite", package = "bamboo")
  txdb <- AnnotationDbi::loadDb(samplefile)
  expect_output(prepareAnnotations(txdb), "test that rows are ordered by tx name")
  expect_type(prepareAnnotations(txdb),"list")
  expect_named(prepareAnnotations(txdb), c("intronsByTxEns" ,"exonsByTx","txIdToGeneIdTable","unlisted_introns" ))

})


test_that("isoform reconstruction of bamFile is summarizedExperiment", {
  samplefile <- system.file("extdata","Homo_sapiens.GRCh38.91.annotations-txdb.sqlite", package = "bamboo")
  txdb <- AnnotationDbi::loadDb(samplefile)
  test.bam <- system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo")
  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo")
  expect_output(isore(bamFile = test.bam,  txdb = txdb,genomeFA = fa.file), "### prepare annotations ###")
  expect_s4_class(isore(bamFile = test.bam,  txdb = txdb,genomeFA = fa.file), "summarizedExperiment")
  expect_equal(isore(bamFile = test.bam,  txdb = txdb,genomeFA = fa.file), sysdata_se)
})


test_that("isoform reconstruction of missing bam file is error", {
  samplefile <- system.file("extdata","Homo_sapiens.GRCh38.91.annotations-txdb.sqlite", package = "bamboo")
  txdb <- AnnotationDbi::loadDb(samplefile)
  test.bam <- system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo")
  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo")
  expect_error(isore(bamFile = '',  txdb = txdb,genomeFA = fa.file),"Bam file is missing from arguments.")
})

test_that("isoform reconstruction of missing both txdb and txdbTablesList object is error", {

  test.bam <- system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo")
  fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bamboo")
  expect_error(isore(bamFile = test.bam,  txdb = NULL, txdbTablesList = NULL, genomeFA = fa.file),"txdb object is missing.")
})

test_that("isoform reconstruction of missing fa file is error", {
  samplefile <- system.file("extdata","Homo_sapiens.GRCh38.91.annotations-txdb.sqlite", package = "bamboo")
  txdb <- AnnotationDbi::loadDb(samplefile)
  test.bam <- system.file("extdata", "GIS_HepG2_cDNAStranded_Rep5-Run4_chr9_108865774_109014097.bam", package = "bamboo")
  expect_error(isore(bamFile = test.bam,  txdbTablesList = txdb,genomeFA = ""))
  expect_error(isore(bamFile = test.bam,  txdbTablesList = txdb,genomeFA = NULL),"GenomeFA file is missing.")
})
