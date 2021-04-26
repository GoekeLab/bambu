context("bambu-processReads_utilityConstructReadClasses")

test_that("assignGeneId returns correct gene ids",{
  se <- readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
  gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))
  expected = read.table(system.file("extdata", "rowData_seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.csv", package = "bambu"), header = TRUE)
  expect_equal(assignGeneIds(rowRanges(se), gr), expected$GENEID)
})
