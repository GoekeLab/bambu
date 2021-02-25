context("bambu-processReads_utilityConstructReadClasses")

test_that("assignGeneId returns correct gene ids",{
  se <- readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/test_se.rds")
  annotations <- readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/test_annotations.rds")
  seExpected = readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000_readClassSe.rds")
  expect_equal(assignGeneIds(rowRanges(se), annotations), rowData(seExpected)$GENEID)
})

#for test purpose
#grl <- IRangesList(IRanges(c(1,90),c(5,95)), IRanges(c(10,20),c(15,25)), IRanges(c(20,30,40), c(25,35,45)),IRanges(c(40,50),c(45,55)), IRanges(c(50,60),c(55,65)), IRanges(c(70,80),c(75,85)),IRanges(c(80,90),c(85,95)), IRanges(c(200,210),c(205,215)), IRanges(c(300),c(305)))
  