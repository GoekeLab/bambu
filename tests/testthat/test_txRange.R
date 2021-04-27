context("txRange, full length transcript prediction")

test_that("txRange generates a gene and transcript score",{
  se <- readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/test_se.rds")
  genomeSequence <- readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/test_genomeSequence.rds")
  annotations <- readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/test_annotations.rds")
  se = txrange.scoreReadClasses(se, genomeSequence, annotations)
  expect_is(rowData(se)$geneScore, class = 'numeric')
  expect_is(rowData(se)$txScore, class = 'numeric')
})

test_that("addRowData adds all the correct rowData features",{
  se <- readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/test_se.rds")
  genomeSequence <- readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/test_genomeSequence.rds")
  annotations <- readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/test_annotations.rds")
  se = addRowData(se, genomeSequence, annotations)
  expect_is(rowData(se)$numExons, class = 'numeric')
  expect_is(rowData(se)$equal, class = 'logical')
  expect_is(rowData(se)$GENEID, class = 'character')
  expect_is(rowData(se)$novel, class = 'logical')
  expect_is(rowData(se)$totalGeneReadProp, class = 'numeric')
  expect_is(rowData(se)$numAstart, class = 'integer')
  expect_is(rowData(se)$numAend, class = 'integer')
  expect_is(rowData(se)$numTstart, class = 'integer')
  expect_is(rowData(se)$numTend, class = 'integer')
})

test_that("calculateGeneProportion",{
  se <- readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/test_se.rds")
  annotations <- readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/test_annotations.rds")
  seExpected = readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000_readClassSe.rds")
  rowData(se)$GENEID = assignGeneIds(rowRanges(se), annotations)
  se = calculateGeneProportion(se)
  expect_equal(rowData(se)$totalGeneReadProp, rowData(seExpected)$totalGeneReadProp)
  
})

test_that("isReadClassEqual",{
  se <- readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/test_se.rds")
  annotations <- readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/test_annotations.rds")
  seExpected = readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000_readClassSe.rds")
  expect_equal(isReadClassEqual(rowRanges(se), annotations),rowData(seExpected)$equal)
})

test_that("countPolyATerminals",{
  se <- readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/test_se.rds")
  genomeSequence <- readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/test_genomeSequence.rds")
  seExpected = readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000_readClassSe.rds")
  se = countPolyATerminals(se, genomeSequence)
  expect_equal(rowData(se)$numAstart,rowData(seExpected)$numAstart)
  expect_equal(rowData(se)$numAend,rowData(seExpected)$numAend)
  expect_equal(rowData(se)$numTstart,rowData(seExpected)$numTstart)
  expect_equal(rowData(se)$numTend,rowData(seExpected)$numTend)
})

test_that("getAgnosticFeatures",{
  #TODO make a expected features dataset to compare to
  se <- readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/test_se.rds")
  features = getAgnosticFeatures(se)
  features = as.data.frame(features)
  expect_is(features$numReads, class = 'numeric')
  expect_is(features$SD, class = 'numeric')
  expect_is(features$SDend, class = 'numeric')
  expect_is(features$geneReadProp, class = 'numeric')
  expect_is(features$tx_strand_bias, class = 'numeric')
  expect_is(features$numAstart, class = 'numeric')
  expect_is(features$numAend, class = 'numeric')
  expect_is(features$numTstart, class = 'numeric')
  expect_is(features$numTend, class = 'numeric')
  expect_is(features$transcriptProp, class = 'numeric')
  expect_is(features$transcriptPropMin, class = 'numeric')
})

test_that("assignGeneIds",{
  se <- readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/test_se.rds")
  annotations <- readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/test_annotations.rds")
  seExpected = readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000_readClassSe.rds")
  expect_equal(assignGeneIds(rowRanges(se), annotations),rowData(seExpected)$GENEID)
})

test_that("prepareGeneModelFeatures",{
  se <- readRDS("C:/Users/simandred/Documents/GitHub/bambu/tests/testData/test_se.rds")
  geneFeatures = prepareGeneModelFeatures(se)
  geneFeatures = as.data.frame(geneFeatures)
  expect_is(features$numReadsLog, class = 'numeric')
  expect_is(features$strand_bias, class = 'numeric')
  expect_is(features$numRCs, class = 'numeric')
  expect_is(features$numNonSubsetRCs, class = 'numeric')
  expect_is(features$numExons, class = 'numeric')
  expect_is(features$isSpliced, class = 'numeric')
  expect_is(features$highConfidence, class = 'numeric')
})
