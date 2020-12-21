context("annotate splice overlap by distance")

test_that("annotateSpliceOverlapByDist can generate a dataframe with annotation",{
  query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
  subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
  refTab <- readRDS(system.file("extdata", "annotateSpliceOverlapsByDist_refoutput.rds", package = "bambu"))
  tab <- annotateSpliceOverlapsByDist(query, subject)
  expect_is(tab, class = 'data.frame')
  expect_equal(tab, refTab[, colnames(tab)])
})
