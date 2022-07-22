context("Generate GTF file from summarizedExperiment object")

# test for the function: writeBambuOutput
test_that("the output files have correct prefix, file format and all of them are placed correctly in the given output path", {
    se <- readRDS(system.file("extdata", "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
    path <- test_path("fixtures")
    prefix <- "replicate5_run1_"
    
    writeBambuOutput(se, path, prefix)
    
    outputFileName <- c("replicate5_run1_extended_annotations.gtf", "replicate5_run1_counts_transcript.txt",
                        "replicate5_run1_counts_transcript.txt", "replicate5_run1_CPM_transcript.txt", 
                        "replicate5_run1_fullLengthCounts_transcript.txt", "replicate5_run1_partialLengthCounts_transcript.txt", 
                        "replicate5_run1_uniqueCounts_transcript.txt", "replicate5_run1_counts_gene.txt")
    
    checkOutput <- sapply(outputFileName, function(name){return(file.exists(test_path("fixtures", name)))}) 
    
    expect_true(all(checkOutput)) # check prefix and file format 
    
    expect_equal(length(list.files(test_path("fixtures"))), 8) # check number of output files at desired output path.
    
    unlink(test_path("fixtures", "*"))
    
})


# test for the function: readFromGTF
test_that("readGTF can generate a GRangesList from a GTF file", {
    se <- readRDS(system.file("extdata", "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
    gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
    path <- tempdir()
    outputGtfFile <- tempfile()
    expect_null(writeBambuOutput(se, path))
    gr <- readFromGTF(gtf.file)
    expect_null(writeToGTF(gr, outputGtfFile))

    expect_s4_class(gr, class = "CompressedGRangesList")
    expect_named(mcols(gr), c("TXNAME", "GENEID"))
})
