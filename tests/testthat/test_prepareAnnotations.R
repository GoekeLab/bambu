context("Prepare annotations")

test_that("prepareAnnotations of txdb object (seqnames, ranges, strand) matches the expectation", {
    txdb <- AnnotationDbi::loadDb(system.file("extdata", "Homo_sapiens.GRCh38.91.annotations-txdb_chr9_1_1000000.sqlite", package = "bambu"))
    
    expectedGR <- readRDS(system.file("extdata", "annotationGranges_txdb2Grch38_91_chr9_1_1000000.rds", package = "bambu"))
    
    gr <- prepareAnnotations(x = txdb)
    
    expect_equal(lapply(gr, granges), lapply(expectedGR, granges), 
                 info = "Test failed: prepareAnnotations of txdb object (seqnames, ranges, strand) does not match the expectation")
})


test_that("prepareAnnotations of a path to gtf file (seqnames, ranges, strand) matches the expectation", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
    
    expectedGR <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))
    
    gr <- prepareAnnotations(x = gtf.file)
    
    expect_equal(lapply(gr, granges), lapply(expectedGR, granges), 
                 info="Test failed: prepareAnnotations of a path to gtf file (seqnames, ranges, strand) does not match the expectation")
})


test_that("prepareAnnotations of txdb object (metadata) matches the expectation", {
    txdb <- AnnotationDbi::loadDb(system.file("extdata", "Homo_sapiens.GRCh38.91.annotations-txdb_chr9_1_1000000.sqlite", package = "bambu"))
  
    expectedGR <- readRDS(system.file("extdata", "annotationGranges_txdb2Grch38_91_chr9_1_1000000.rds", package = "bambu"))
  
    gr <- prepareAnnotations(x = txdb)
    
    expect_named(mcols(gr), c("TXNAME", "GENEID", "txid", "eqClass", "eqClassById"), 
                 info="Test failed: prepareAnnotations of txdb object (metadata column name) does not match the expectation")
    expect_equal(mcols(gr), mcols(expectedGR), 
                 info="Test failed: prepareAnnotations of txdb object (metadata) does not match the expectation")
})


test_that("prepareAnnotations of a path to gtf file (metadata) matches the expectation", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
  
    expectedGR <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))
  
    gr <- prepareAnnotations(x = gtf.file)
  
    expect_named(mcols(gr), c("TXNAME", "GENEID", "txid", "eqClass", "eqClassById"), 
                 info="Test failed: prepareAnnotations of a path to gtf file (metadata column name) does not match the expectation")
    expect_equal(mcols(gr), mcols(expectedGR), 
                 info="Test failed: prepareAnnotations of a path to gtf file (metadata) does not match the expectation")
})


test_that("positive strand gives ascending exon_rank", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
  
    gr <- prepareAnnotations(x = gtf.file)
    
    unlisted_gr <- unlist(gr)
    
    positive_check <- data.frame(txname = names(unlisted_gr), unlisted_gr) %>% 
        filter(strand == "+") %>% 
        group_by(txname) %>% 
        summarise(validate = all(exon_rank == seq(n())))
    
    expect_true(all(positive_check$validate), 
                info="Test failed: positive strand does not give ascending exon_rank")
})


test_that("positive strand gives descending exon_endRank", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
  
    gr <- prepareAnnotations(x = gtf.file)
    
    unlisted_gr <- unlist(gr)
    
    positive_check <- data.frame(txname = names(unlisted_gr), unlisted_gr) %>% 
        filter(strand == "+") %>% 
        group_by(txname) %>% 
        summarise(validate = all(exon_endRank == seq(n(),1)))
  
    expect_true(all(positive_check$validate), 
                info = "Test failed: positive strand does not give descending exon_endRank" )
})


test_that("negative strand gives descending exon_rank", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
  
    gr <- prepareAnnotations(x = gtf.file)
    unlisted_gr <- unlist(gr)
    
    negative_check <- data.frame(txname = names(unlisted_gr), unlisted_gr) %>% 
        filter(strand == "-") %>% 
        group_by(txname) %>% 
        summarise(validate = all(exon_rank == seq(n(),1)))
  
    expect_true(all(negative_check$validate), 
                info="Test failed: negative strand does not give descending exon_rank")
})


test_that("negative strand gives ascending exon_endRank", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")

    gr <- prepareAnnotations(x = gtf.file)
    
    unlisted_gr <- unlist(gr)
    
    negative_check <- data.frame(txname = names(unlisted_gr), unlisted_gr) %>% 
        filter(strand == "-") %>% 
        group_by(txname) %>% 
        summarise(validate = all(exon_endRank == seq(n())))
  
    expect_true(all(negative_check$validate), 
                info="Test failed: negative strand does not give ascending exon_endRank")
})


test_that("txid must be in EqClassById", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
  
    gr <- prepareAnnotations(x = gtf.file)
    
    check <- data.frame(mcols(gr)) %>% 
        dplyr::select(txid, eqClassById) %>% 
        rowwise() %>% 
        mutate(validate = txid %in% eqClassById)
    
    expect_true(all(check$validate), 
                info = "Test failed: txid for some transcripts are not in EqClassById")
})


# This function will pick a few genes, get their transcripts, and verify the 
# eqClassById (transcript equivalence class) for each selected transcript.
test_that("eqClassById is correct (tested for a few genes)", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
    
    gr <- prepareAnnotations(x = gtf.file)  
    
    set.seed(42) # Ensure the test is consistent
    gene <- sample(na.omit(mcols(gr)$GENEID), 10) # Pick a few genes 
    startendremoved <- cutStartEndFromGrangesList(gr)
    transcript_to_test <- data.frame(mcols(gr)) %>% 
        filter(GENEID %in% gene) %>% # filter according to selected gene
        dplyr::select(TXNAME, txid, eqClassById) %>% 
        rowwise() %>% 
        mutate(validate = list(all(sapply(eqClassById, # Check whether eqClassById is correct
                            function(id){length(findOverlaps(startendremoved[[txid]], 
                            startendremoved[[id]])) == length(startendremoved[[txid]])}))))
    expect_true(all(unlist(transcript_to_test$validate)), 
                info="Test failed: eqClassById for some transcripts are incorrect")
  
})


test_that("eqClass and eqClassById matches", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
    
    gr <- prepareAnnotations(x = gtf.file)    

    # convert the eqClassById to eqClass
    convert <- data.frame(mcols(gr)) %>% 
        dplyr::select(TXNAME, eqClass, eqClassById) %>% 
        mutate(validate=sapply(eqClassById, function(idlist){paste(sort(TXNAME[idlist]), collapse=".")}))

    expect_equal(convert$validate,convert$eqClass, 
                 info="Test failed: eqClass and eqClassById do not match")
})
