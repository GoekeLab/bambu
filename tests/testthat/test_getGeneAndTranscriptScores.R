context('Get gene/transcript score and feature preparation')

test_that('getTranscriptScore returns expected dataframe of txScore and txFDR', {
  se_rowData = read.csv(system.file('extdata', 'rowData_seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.csv', package='bambu'), sep=' ')
  # Filter to keep only RCs with >= 2 reads
  se_rowData = se_rowData[se_rowData$readCount>=2,]
  
  txScore_df = getTranscriptScore(se_rowData)
  
  # Load in reference dataframe
  ref_txScore_df = read.csv(system.file('extdata', 'txScore_df.csv', package='bambu'))
  
  expect_equal(txScore_df, ref_txScore_df)
})

test_that('getGeneScore returns expected dataframe of geneScore and geneFDR',{
  se_rowData = read.csv(system.file('extdata', 'rowData_seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.csv', package='bambu'), sep=' ')
  # Filter to keep only RCs with >= 2 reads
  se_rowData = se_rowData[se_rowData$readCount>=2,]
  
  geneScore_df = getGeneScore(se_rowData)
  
  # Load in reference dataframe
  ref_geneScore_df = read.csv(system.file('extdata', 'geneScore_df.csv', package='bambu'))
  
  expect_equal(geneScore_df, ref_geneScore_df)
})

test_that('prepareTranscriptModelFeature returns expected tibble of features', {
  se_rowData = read.csv(system.file('extdata', 'rowData_seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.csv', package='bambu'), sep=' ')
  # Filter to keep only RCs with >= 2 reads
  se_rowData = se_rowData[se_rowData$readCount>=2,]
  
  features = prepareTranscriptModelFeatures(se_rowData)
  
  # Load reference tibble
  ref_features = readRDS(system.file('extdata', 'prepareTranscriptModelFeatures.rds', package='bambu'))

  expect_equal(features, ref_features)
})