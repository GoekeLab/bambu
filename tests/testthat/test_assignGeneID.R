context("bambu-processReads_utilityConstructReadClasses")

test_that("assignGeneId returns correct gene ids",{
  se <- readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/test_se.rds")
  annotations <- readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/test_annotations.rds")
  seExpected = readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000_readClassSe.rds")
  expect_equal(assignGeneIds(rowRanges(se), annotations), rowData(seExpected)$GENEID)
})