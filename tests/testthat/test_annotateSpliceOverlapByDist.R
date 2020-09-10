context("annotate splice overlap by distance")

test_that("annotateSpliceOverlapByDist can generate a dataframe with annotation",{
  query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
  subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
  tab <- annotateSpliceOverlapByDist(query, subject)
  expect_is(tab, class = 'data.frame')
  expect_equal(names(tab), c("queryId","subjectId","strand","intronRetention.subject","intronRetention.query",
                             "exonSkipping.query","exonSkipping.subject","alternativeFirstExon","alternativeLastExon",     
                             "alternativeTSS","alternativeTES","internalFirstExon.query","internalLastExon.query",
                             "internalFirstExon.subject","internalLastExon.subject","exon5Prime","exon3Prime"))
})
