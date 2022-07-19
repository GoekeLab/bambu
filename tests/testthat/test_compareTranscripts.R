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
# compareTranscripts function. These test cases are selected in a way to cover as
# many scenarios as possible. 

test_that("compareTranscripts gives correct output about exon skipping (query) 
          and alternative TSS for a transcript along positive strand", {
    
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    # transcripts to inspect & compare
    t1 <- query[1]
    t2 <- subject[1]
    
    # check alternative TSS
    expect_equal(tab[1,]$alternativeTSS, -(3625050 - 3625015)) 
    
    # check exon skipping (query)
    expect_equal(tab[1,]$exonSkipping.query, 1) 
    
})


test_that("compareTranscripts gives correct output about intron retention (subject), 
          alternative TSS and alternative TES for a transcript along positive strand", {
    
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    # transcripts to inspect & compare
    t1 <- query[2]
    t2 <- subject[2]
    
    # check alternative TSS
    expect_equal(tab[2,]$alternativeTSS, -(26280086 - 26280122))
    
    # check alternative TES 
    expect_equal(tab[2,]$alternativeTES, 26281522 - 26281450)
    
    # check intron retention (subject)
    expect_equal(tab[2,]$intronRetention.subject, 1)
    
})


test_that("compareTranscripts gives correct output about alternative first exon, exon skipping (query), 
          exon skipping (subject), alternative TSS and internalFirstExon (subject) for a transcript along positive strand", {
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    # transcripts to inspect & compare
    t1 <- query[3]
    t2 <- subject[3]

    # check alternative first exon 
    expect_equal(tab[3,]$alternativeFirstExon, TRUE)
    
    # check alternative TES 
    expect_equal(tab[3,]$alternativeTES, 237553996 - 237553994)
    
    # check internal first exon (subject)
    expect_equal(tab[3,]$internalFirstExon.subject, TRUE)

    # check exon skipping (query)
    expect_equal(tab[3,]$exonSkipping.query, 2)
    
    # check exon skipping (subject)
    expect_equal(tab[3,]$exonSkipping.subject, 0)
    
})


test_that("compareTranscripts gives correct output about alternative last exons and intron retention (query) for a transcript along positive strand", {
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    # transcripts to inspect & compare
    t1 <- query[4]
    t2 <- subject[4]
    
    # check alternative last exons 
    expect_equal(tab[4,]$alternativeLastExon, TRUE)
    
    # check intron retention (query)
    expect_equal(tab[4, ]$intronRetention.query, 0)
    
})


test_that("compareTranscripts gives correct output about alternative first exon, alternative last exon, 
          exon skipping (subject) and internal first exon (query) for a transcript along positive strand", {
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    # transcripts to inspect & compare
    t1 <- query[5]
    t2 <- subject[5]
    
    # check alternative first exon 
    expect_equal(tab[5,]$alternativeFirstExon, TRUE)
    
    # check alternative last exon 
    expect_equal(tab[5,]$alternativeLastExon, TRUE)
    
    # check exon skipping (subject)
    expect_equal(tab[5,]$exonSkipping.subject, 1)
    
})


test_that("compareTranscripts gives correct output about exon skipping (query) and 
          alternative TSS for a transcript along negative strand", {
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    # transcripts to inspect & compare
    t1 <- query[6]
    t2 <- subject[6]

    # check exon skipping (query)
    expect_equal(tab[6,]$exonSkipping.query, 1)
    
    # check alternative TSS
    expect_equal(tab[6,]$alternativeTSS, 144413569 - 144413586)
    
})


test_that("compareTranscripts gives correct output about alternative 3' exon splice site, exon skipping (query), 
          alternative TSS and alternative TES for a transcript along positive strand", {
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    # transcripts to inspect & compare
    t1 <- query[7]
    t2 <- subject[7]

    # check exon skipping (query)
    expect_equal(tab[7,]$exonSkipping.query, 2)
    
    # check alternative TSS 
    expect_equal(tab[7,]$alternativeTSS, -(112542226 - 112541915))
    
    # check alternative TES
    expect_equal(tab[7,]$alternativeTES, 112577058 - 112576299)
    
    # check alternative 3' exon splice site
    expect_equal(tab[7,]$exon3Prime, 1)
    
})


test_that("compareTranscripts gives correct output about alternative TSS, alternative last exon, internal last exon (query), 
          alternative exon 3' end and alternative exon 5' end for a transcript along negative strand", {
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    # transcripts to inspect & compare
    t1 <- query[8]
    t2 <- subject[8]
    
    # check alternative TSS
    expect_equal(tab[8,]$alternativeTSS, 47426225 - 47426266)
    
    # check alternative last exon
    expect_equal(tab[8,]$alternativeLastExon, TRUE)
    
    # check internal last exon
    expect_equal(tab[8,]$internalLastExon.query, TRUE)
    
    # check alternative 3' exon splice site
    expect_equal(tab[8,]$exon3Prime, 1)
    
    # check alternative 5' exon splice site
    expect_equal(tab[8,]$exon5Prime, 0)
    
})


test_that("compareTranscripts gives correct output about alternative TSS, alternative exon 3' end, 
          alternative exon 5' end and alternative TES for a transcript along positive strand", {
    query <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testQuery.rds", package = "bambu"))
    subject <- readRDS(system.file("extdata", "annotateSpliceOverlapByDist_testSubject.rds", package = "bambu"))
    tab <- compareTranscripts(query, subject)
    
    # transcripts to inspect & compare
    t1 <- query[9]
    t2 <- subject[9]
    
    # check alternative TSS 
    expect_equal(tab[9,]$alternativeTSS, -(44070892 - 44070751))
    
    # check alternative TES 
    expect_equal(tab[9,]$alternativeTES, 44075841 - 44076300)
    
    # check alternative 3' exon splice site 
    expect_equal(tab[9,]$exon3Prime, 2)
    
    # check alternative 5' exon splice site
    expect_equal(tab[9,]$exon5Prime, 1)
    
})
    