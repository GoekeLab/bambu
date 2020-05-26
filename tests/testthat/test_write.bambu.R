context("Generate GTF file from summarizedExperiment object")

test_that("write.bambu can generate a GTF, transcript count file, and a gene count file by running write.gtf",{
  se <- readRDS(system.file("extdata", "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  path <- "test_writebambu_gtf"
  gtf <- write.bambu(se,path)
  expect_is(gtf, "data.frame")
  expect_equal(ncol(gtf), 9)
  })
