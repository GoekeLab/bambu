context("Prepare annotations")

test_that("prepareAnnotations of txdb object (seqnames, ranges, strand) matches the expectation", {
    txdb <- AnnotationDbi::loadDb(system.file("extdata", "Homo_sapiens.GRCh38.91.annotations-txdb_chr9_1_1000000.sqlite", package = "bambu"))
    
    expectedGR <- readRDS(system.file("extdata", "annotationGranges_txdb2Grch38_91_chr9_1_1000000.rds", package = "bambu"))
    
    gr <- prepareAnnotations(x = txdb)
    
    expect_equal(lapply(gr, granges), lapply(expectedGR, granges))
})


test_that("prepareAnnotations of a path to gtf file (seqnames, ranges, strand) matches the expectation", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
    
    expectedGR <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))
    
    gr <- prepareAnnotations(x = gtf.file)
    
    expect_equal(lapply(gr, granges), lapply(expectedGR, granges))
})


test_that("prepareAnnotations of txdb object (metadata) matches the expectation", {
    txdb <- AnnotationDbi::loadDb(system.file("extdata", "Homo_sapiens.GRCh38.91.annotations-txdb_chr9_1_1000000.sqlite", package = "bambu"))
  
    expectedGR <- readRDS(system.file("extdata", "annotationGranges_txdb2Grch38_91_chr9_1_1000000.rds", package = "bambu"))
  
    gr <- prepareAnnotations(x = txdb)
    
    expect_named(mcols(gr), c("TXNAME", "GENEID", "txid", "eqClass", "eqClassById"))
    expect_equal(mcols(gr), mcols(expectedGR))
})


test_that("prepareAnnotations of a path to gtf file (metadata) matches the expectation", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
  
    expectedGR <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))
  
    gr <- prepareAnnotations(x = gtf.file)
  
    expect_named(mcols(gr), c("TXNAME", "GENEID", "txid", "eqClass", "eqClassById"))
    expect_equal(mcols(gr), mcols(expectedGR))
})


test_that("positive strand gives ascending exon_rank", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
  
    gr <- prepareAnnotations(x = gtf.file)
  
    positive_check <- data.frame(txname = names(unlist(gr)), unlist(gr)) %>% 
        filter(strand == "+") %>% 
        group_by(txname) %>% 
        summarise(validate = all(exon_rank == seq(n())))
    
    expect_true(all(positive_check$validate))
})


test_that("positive strand gives descending exon_endRank", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
  
    gr <- prepareAnnotations(x = gtf.file)
  
    positive_check <- data.frame(txname = names(unlist(gr)), unlist(gr)) %>% 
        filter(strand == "+") %>% 
        group_by(txname) %>% 
        summarise(validate = all(exon_endRank == seq(n(),1)))
  
    expect_true(all(positive_check$validate))
})


test_that("negative strand gives descending exon_rank", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
  
    gr <- prepareAnnotations(x = gtf.file)
  
    negative_check <- data.frame(txname = names(unlist(gr)), unlist(gr)) %>% 
        filter(strand == "-") %>% 
        group_by(txname) %>% 
        summarise(validate = all(exon_rank == seq(n(),1)))
  
    expect_true(all(negative_check$validate))
})


test_that("negative strand gives ascending exon_endRank", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")

    gr <- prepareAnnotations(x = gtf.file)
  
    negative_check <- data.frame(txname = names(unlist(gr)), unlist(gr)) %>% 
        filter(strand == "-") %>% 
        group_by(txname) %>% 
        summarise(validate = all(exon_endRank == seq(n())))
  
    expect_true(all(negative_check$validate))
})


test_that("txid must be in EqClassById", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
  
    gr <- prepareAnnotations(x = gtf.file)
    
    check <- data.frame(mcols(gr)) %>% 
        dplyr::select(txid, eqClassById) %>% 
        rowwise() %>% 
        mutate(validate = txid %in% eqClassById)
    
    expect_true(all(check$validate))
})


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
    expect_true(all(unlist(transcript_to_test$validate)))
  
})


test_that("eqClass and eqClassById matches", {
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
    
    gr <- prepareAnnotations(x = gtf.file)    

    # convert the eqClassById to eqClass
    convert <- data.frame(mcols(gr)) %>% 
        dplyr::select(TXNAME, eqClass, eqClassById) %>% 
        mutate(validate=sapply(eqClassById, function(idlist){paste(sort(TXNAME[idlist]), collapse=".")}))

    expect_equal(convert$validate,convert$eqClass)
})
