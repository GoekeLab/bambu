context("Generate GTF from summarizedExperiment object")

test_that("SEtoGTF can generate a GTF file",{
  se <- readRDS(system.file("extdata", "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  gtf <- SEtoGTF(se)
  expect_is(gtf, "data.frame")
  expect_named(ncol(gtf), 9)
  })
