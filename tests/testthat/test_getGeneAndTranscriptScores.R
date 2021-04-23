context('Get gene/transcript score and feature preparation')

test_that('getTranscriptScore returns expected dataframe of txScore and txFDR', {
  se = readRDS(system.file('extdata', 'se_forCalculatingGeneAndTranscriptScore.rds', package='bambu'))

  txScore_df = getTranscriptScore(rowData(se))
  
  # Load in reference dataframe
  ref_txScore_df = read.csv(system.file('extdata', 'txScore_df.csv', package='bambu'))
  
  expect_equal(txScore_df, ref_txScore_df)
})

test_that('getGeneScore returns expected dataframe of geneScore and geneFDR',{
  se = readRDS(system.file('extdata', 'se_forCalculatingGeneAndTranscriptScore.rds', package='bambu'))
  
  geneScore_df = getGeneScore(rowData(se))
  
  # Load in reference dataframe
  ref_geneScore_df = read.csv(system.file('extdata', 'geneScore_df.csv', package='bambu'))
  
  expect_equal(geneScore_df, ref_geneScore_df)
})

test_that('prepareTranscriptModelFeature returns expected tibble of features', {
  se = readRDS(system.file('extdata', 'se_forCalculatingGeneAndTranscriptScore.rds', package='bambu'))

  features = prepareTranscriptModelFeatures(rowData(se))
  
  # Load reference tibble
  ref_features = readRDS(system.file('extdata', 'prepareTranscriptModelFeatures.rds', package='bambu'))

  expect_equal(features, ref_features)
})