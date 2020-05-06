context("Gene quantification")

test_that("transcriptToGeneExpression (reducing transcript to gene expression) runs successfully",{
  se <- readRDS(system.file("extdata", "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  seExtended <- readRDS(system.file("extdata", "seOutputExtended_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  seCombined <- readRDS(system.file("extdata", "seOutputCombined_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  seCombinedExtended <- readRDS(system.file("extdata", "seOutputCombinedExtended_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))

  # 1. test single sample case, non extended
  seGene <- transcriptToGeneExpression(se)
  expect_s4_class(seGene, "SummarizedExperiment")
  expect_equal(seGene,seGeneExpected)

  # 2. test single sample case, extended
  seExtendedGene <- transcriptToGeneExpression(seExtended)
  expect_s4_class(seExtendedGene, "SummarizedExperiment")
  expect_equal(seExtendedGene,seExtendedGeneExpected)

  # 3. test mutiple sample case, non extended
  seCombinedGene <- transcriptToGeneExpression(seCombined)
  expect_s4_class(seCombinedGene, "SummarizedExperiment")
  expect_equal(seCombinedGene,seCombinedGeneExpected)

  # 4. test multiple sample case, extended
  seCombinedExtendedGene <- transcriptToGeneExpression(seCombinedExtended)
  expect_s4_class(seCombinedExtendedGene, "SummarizedExperiment")
  expect_equal(seCombinedExtendedGene,seCombinedExtendedGeneExpected)

})
