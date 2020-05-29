context("Generate GTF file from summarizedExperiment object")

test_that("readGTF can generate a GRangesList from a GTF file",{
  se <- readRDS(system.file("extdata", "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
  path <- tempdir()
  outputGtfFile <- tempfile()
  expect_null(writeBambuOutput(se,path))
  gr <- readFromGTF(gtf.file)
  expect_null(writeToGTF(gr, outputGtfFile))
  
  expect_s4_class(gr, class = 'CompressedGRangesList')
  expect_named(mcols(gr), c("TXNAME", "GENEID"))
})
