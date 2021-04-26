context('Preparing model features')

test_that("prepareGeneModelFeature returns expected tibble of features",{
  se_rowData = read.csv(system.file('extdata', 'rowData_seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.csv', package='bambu'), sep=' ')
  # Filter to keep only RCs with >= 2 reads
  se_rowData = se_rowData[se_rowData$readCount>=2,]
  
  geneFeatures = prepareGeneModelFeatures(se_rowData)
  
  # Load reference tibble
  ref_features = readRDS(system.file('extdata', 'prepareGeneModelFeatures.rds', package='bambu'))
 
  expect_equal(geneFeatures, ref_features)
})

test_that('prepareTranscriptModelFeature returns expected tibble of features', {
  se_rowData = read.csv(system.file('extdata', 'rowData_seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.csv', package='bambu'), sep=' ')
  # Filter to keep only RCs with >= 2 reads
  se_rowData = se_rowData[se_rowData$readCount>=2,]
  
  transcriptFeatures = prepareTranscriptModelFeatures(se_rowData)
  
  # Load reference tibble
  ref_features = readRDS(system.file('extdata', 'prepareTranscriptModelFeatures.rds', package='bambu'))
  
  expect_equal(transcriptFeatures, ref_features)
})