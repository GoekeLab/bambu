context("Prepare annotations")

# This test additionally saves the output of prepareAnnotations so other unit test can reuse it.
test_that("prepareAnnotations works properly on a txdb object", {
    txdb <- AnnotationDbi::loadDb(system.file("extdata", "Homo_sapiens.GRCh38.91.annotations-txdb_chr9_1_1000000.sqlite", package = "bambu"))
    grTXDB <- prepareAnnotations(txdb) 
  
    expect_s4_class(grTXDB, class = "CompressedGRangesList")

    saveRDS(grTXDB, test_path("fixtures", "grTXDB.rds"))
})


# This test additionally saves the output of prepareAnnotations so other unit test can reuse it.
test_that("prepareAnnotations works properly on a path to gtf file", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
    grGTF <- prepareAnnotations(gtf.file)
  
    expect_s4_class(grGTF, class = "CompressedGRangesList")
  
    saveRDS(grGTF, test_path("fixtures", "grGTF.rds"))
})


test_that("prepareAnnotations of txdb object (seqnames, ranges, strand) matches the expectation", {
    gr <- readRDS(test_path("fixtures", "grTXDB.rds"))
    
    expectedGR <- readRDS(system.file("extdata", "annotationGranges_txdb2Grch38_91_chr9_1_1000000.rds", package = "bambu"))
    
    expect_equal(lapply(gr, granges), lapply(expectedGR, granges))
})


test_that("prepareAnnotations of a path to gtf file (seqnames, ranges, strand) matches the expectation", {
    gr <- readRDS(test_path("fixtures", "grGTF.rds"))
    
    expectedGR <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))
    
    expect_equal(lapply(gr, granges), lapply(expectedGR, granges))
})


test_that("prepareAnnotations of txdb object (metadata) matches the expectation", {
    gr <- readRDS(test_path("fixtures", "grTXDB.rds"))
    
    expectedGR <- readRDS(system.file("extdata", "annotationGranges_txdb2Grch38_91_chr9_1_1000000.rds", package = "bambu"))
    
    expect_named(mcols(gr), c("TXNAME", "GENEID", "txid", "eqClass", "eqClassById"))
    expect_equal(mcols(gr), mcols(expectedGR))
})


test_that("prepareAnnotations of a path to gtf file (metadata) matches the expectation", {
    gr <- readRDS(test_path("fixtures", "grGTF.rds"))
    
    expectedGR <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))
  
    expect_named(mcols(gr), c("TXNAME", "GENEID", "txid", "eqClass", "eqClassById"))
    expect_equal(mcols(gr), mcols(expectedGR))
})


test_that("positive strand gives ascending exon_rank", {
    gr <- readRDS(test_path("fixtures", "grGTF.rds"))
    
    unlisted_gr <- unlist(gr)
    
    positive_check <- data.frame(txname = names(unlisted_gr), unlisted_gr) %>% 
        filter(strand == "+") %>% 
        group_by(txname) %>% 
        summarise(validate = all(exon_rank == seq(n())))
    
    expect_true(all(positive_check$validate))
})


test_that("positive strand gives descending exon_endRank", {
    gr <- readRDS(test_path("fixtures", "grGTF.rds"))
    
    unlisted_gr <- unlist(gr)
    
    positive_check <- data.frame(txname = names(unlisted_gr), unlisted_gr) %>% 
        filter(strand == "+") %>% 
        group_by(txname) %>% 
        summarise(validate = all(exon_endRank == seq(n(),1)))
  
    expect_true(all(positive_check$validate))
})


test_that("negative strand gives descending exon_rank", {
    gr <- readRDS(test_path("fixtures", "grGTF.rds"))
    
    unlisted_gr <- unlist(gr)
    
    negative_check <- data.frame(txname = names(unlisted_gr), unlisted_gr) %>% 
        filter(strand == "-") %>% 
        group_by(txname) %>% 
        summarise(validate = all(exon_rank == seq(n(),1)))
  
    expect_true(all(negative_check$validate))
})


test_that("negative strand gives ascending exon_endRank", {
    gr <- readRDS(test_path("fixtures", "grGTF.rds"))
    
    unlisted_gr <- unlist(gr)
    
    negative_check <- data.frame(txname = names(unlisted_gr), unlisted_gr) %>% 
        filter(strand == "-") %>% 
        group_by(txname) %>% 
        summarise(validate = all(exon_endRank == seq(n())))
  
    expect_true(all(negative_check$validate))
})


test_that("txid must be in EqClassById", {
    gr <- readRDS(test_path("fixtures", "grGTF.rds"))
    
    check <- data.frame(mcols(gr)) %>% 
        dplyr::select(txid, eqClassById) %>% 
        rowwise() %>% 
        mutate(validate = txid %in% eqClassById)
    
    expect_true(all(check$validate))
})


# This function will get the transcripts, and verify the 
# eqClassById (transcript equivalence class) for each selected transcript.
test_that("eqClassById is correct", {
  gr <- readRDS(test_path("fixtures", "grGTF.rds"))  
  splitEqClassById <- data.frame(mcols(gr)) %>% 
      dplyr::select(txid, eqClassById) %>% 
      tidyr::unnest_longer(eqClassById)
  
  check <- compareTranscripts(gr[splitEqClassById$txid,], 
                              gr[splitEqClassById$eqClassById,]) %>% 
      filter(queryId != subjectId) %>% # no need to compare identical transcript
      mutate(validate = (alternativeFirstExon == FALSE | internalFirstExon.query == TRUE) 
           & (alternativeLastExon == FALSE | internalLastExon.query == TRUE)
           & (exonSkipping.query == 0)
           & (exonSkipping.subject == 0)
           & (exon5Prime == 0)
           & (exon3Prime == 0)
           & (intronRetention.subject == 0))
  #note: intronRetention.query == 0 is allowed as this might be due to longer first and last exons, which are not considered for the equivalent class calculation.
  expect_true(all(check$validate))
})


test_that("eqClass and eqClassById matches", {
    gr <- readRDS(test_path("fixtures", "grGTF.rds"))    

    # convert the eqClassById to eqClass
    convert <- data.frame(mcols(gr)) %>% 
        dplyr::select(TXNAME, eqClass, eqClassById) %>% 
        mutate(validate=sapply(eqClassById, function(idlist){paste(sort(TXNAME[idlist]), collapse=".")}))

    expect_equal(convert$validate,convert$eqClass)

})

# delete the test file
unlink(test_path("fixtures", "grTXDB.rds"))
unlink(test_path("fixtures", "grGTF.rds"))
