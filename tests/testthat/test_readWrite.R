context("Generate GTF file from summarizedExperiment object")

test_that("read.gtf can generate a GRangesList from a GTF file",{
  se <- readRDS(system.file("extdata", "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  path <- "test_writebambu_gtf"
  gtf <- write.bambu(se,path)
  gr <- read.gtf(gtf)
  
  expect_s4_class(gr, class = 'CompressedGRangesList')
  expect_named(mcols(gr), c("TXNAME", "GENEID"))
})
