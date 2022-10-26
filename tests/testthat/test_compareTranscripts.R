context("annotate splice overlap by distance")

test_that("compareTranscripts between the query transcripts and subject transcripts match the expectations",{
    
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    refTab <- readRDS(system.file("extdata", "annotateSpliceOverlapsByDist_refoutput.rds", package = "bambu"))
  
    expect_is(tab, class = 'data.frame')
    expect_equal(refTab, tab)
})

# For the remaining of the test, we select a few test case to validate the output from the 
# compareTranscripts function. They are selected in a way to capture as many scenarios as possible.
# These test cases are shown below. 

### examples for test purposes
## Expected annotations of transcripts used in test query
# 'ENST00000344579', # exon skipping, alternative TSS (-35), +, ENSG00000158109
# 'ENST00000270792', # intron retention subject 1(last exon),alt.TSS,alt.TES, +,
# 'ENST00000410032', # alternative first exon, exon skipping query: 2, 
#                    # exon skipping subject: 0, alternative TES (2bp only), 
#                    # internalFirstExon.subject +
# 'ENST00000468178', # alternative last exon +
# 'ENST00000485956', # alternative first exon, alternative last exon,
# #exon skipping subject = 1, internal first exon query, +
# 'ENST00000530807', # exon skipping query 1, alternative TSS (-17),  -
# 'ENST00000409894', # alternative 3' exon splice site, exon skipping query 2,
# #alternative TSS, alterantive TES, +, ENSG00000125630
# 'ENST00000524447',  # alternative TSS, alternative last exon (internal), 
# #alternative exon 3' end,-, ENSG00000165916
# 'ENST00000591696' # alternative TSS, alternative 3' exon (2), 
# #alternative 5' exon (1) alternative TES, ,+,ENSG00000141349

test_that("the strand column of compareTranscripts matches the expectations",{
    
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    expect_equal(tab$strand, c("+", "+", "+", "+", "+", "-", "+", "-", "+"))
    
})


test_that("the alternativeFirstExon column of compareTranscripts matches the expectations",{
  
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    expect_equal(tab$alternativeFirstExon, c(FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE))
    
})


test_that("the alternativeTSS column of compareTranscripts matches the expectations",{
    
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    expect_equal(tab$alternativeTSS, c(-35, 36, 0, 0, 0, -17, -311, -41, -141))
    
})


test_that("the internalFirstExon.query column of compareTranscripts matches the expectations",{
    
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    expect_equal(as.logical(tab$internalFirstExon.query), c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE))
    
})


test_that("the internalFirstExon.subject column of compareTranscripts matches the expectations",{
    
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    expect_equal(as.logical(tab$internalFirstExon.subject), c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
    
})


test_that("the alternativeLastExon column of compareTranscripts matches the expectations",{
    
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    expect_equal(tab$alternativeLastExon, c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE))
    
})


test_that("the alternativeTES column of compareTranscripts matches the expectations",{
    
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    expect_equal(tab$alternativeTES, c(0, 72, 2, 0, 0, 0, 759, 0, -459))
    
})


test_that("the internalLastExon.query column of compareTranscripts matches the expectations",{
    
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    expect_equal(as.logical(tab$internalLastExon.query), c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE))
    
})


test_that("the internalLastExon.subject column of compareTranscripts matches the expectations",{
    
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    expect_equal(as.logical(tab$internalLastExon.subject), c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
    
})


test_that("the intronRetention.subject column of compareTranscripts matches the expectations",{
    
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    expect_equal(as.numeric(tab$intronRetention.subject), c(0, 1, 0, 0, 0, 0, 0, 0, 0))
    
})


test_that("the intronRetention.query column of compareTranscripts matches the expectations",{
    
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    expect_equal(as.numeric(tab$intronRetention.query), c(0, 0, 0, 0, 0, 0, 0, 0, 0))
    
})


test_that("the exonSkipping.query column of compareTranscripts matches the expectations",{
    
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    expect_equal(tab$exonSkipping.query, c(1, 0, 2, 0, 0, 1, 2, 0, 0))
    
})


test_that("the exonSkipping.subject column of compareTranscripts matches the expectations",{
    
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    expect_equal(tab$exonSkipping.subject, c(0, 0, 0, 0, 1, 0, 0, 0, 0))
    
})


test_that("the exon5Prime column of compareTranscripts matches the expectations",{
    
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    expect_equal(as.numeric(tab$exon5Prime), c(0, 0, 0, 0, 0, 0, 0, 0, 1))
    
})


test_that("the exon3Prime column of compareTranscripts matches the expectations",{
    
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    expect_equal(as.numeric(tab$exon3Prime), c(0, 0, 0, 0, 0, 0, 1, 1, 2))
    
})
