context('Calculation of False Discovery Rate (FDR)')

test_that('Calculated FDR scores are as expected', {
  scores = as.numeric(readLines(system.file('extdata','calculateFDR_arg_scores.txt', package='bambu')))
  labels = as.numeric(readLines(system.file('extdata','calculateFDR_arg_labels.txt', package='bambu')))
  FDR = calculateFDR(scores, labels)
  
  ref_FDR = unname(unlist(read.csv(system.file('extdata','FDR.csv', package='bambu'))))
  
  expect_equal(FDR, ref_FDR)
})