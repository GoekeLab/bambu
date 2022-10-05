context("Generate GTF file from summarizedExperiment object")

# test for the function: writeBambuOutput
test_that("the output files of writeBambuOutput have correct prefix, file format and all of them are placed correctly in the given output path", {
    se <- readRDS(system.file("extdata", "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
    path <- test_path("fixtures")
    prefix <- "replicate5_run1_"
    
    writeBambuOutput(se, path, prefix)
    
    outputFileName <- c("replicate5_run1_extended_annotations.gtf", "replicate5_run1_counts_transcript.txt",
                        "replicate5_run1_counts_transcript.txt", "replicate5_run1_CPM_transcript.txt", 
                        "replicate5_run1_fullLengthCounts_transcript.txt", "replicate5_run1_partialLengthCounts_transcript.txt", 
                        "replicate5_run1_uniqueCounts_transcript.txt", "replicate5_run1_counts_gene.txt")
    
    checkOutput <- outputFileName %in% list.files(path)
    
    expect_true(all(checkOutput))
    
    unlink(test_path("fixtures", "*"))
    
})


test_that("the gene count output file from writeBambuOutput has GENEID column", {
    se <- readRDS(system.file("extdata", "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
    path <- test_path("fixtures")
    
    writeBambuOutput(se, path)
    
    geneCountOutputFilesPath <- test_path("fixtures", "counts_gene.txt")
    
    expect_true("GENEID" %in% colnames(read.table(geneCountOutputFilesPath, header = TRUE)))
    
    unlink(test_path("fixtures", "*"))
    
})


test_that("the transcript count output files from writeBambuOutput have TXNAME and GENEID columns", {
    se <- readRDS(system.file("extdata", "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
    path <- test_path("fixtures")
    
    writeBambuOutput(se, path)
    
    outputFilesName <- list.files(path)
    transcriptCountOutputFilesPath <- test_path("fixtures", outputFilesName[!grepl("counts_gene.txt|extended_annotations.gtf", outputFilesName)])
    
    check <- sapply(seq_along(transcriptCountOutputFilesPath), function(i){
      df <- read.table(transcriptCountOutputFilesPath[i], header = TRUE)
      all(c("TXNAME", "GENEID") %in% colnames(df))
    })
    
    expect_true(all(check))
    
    unlink(test_path("fixtures", "*"))
    
})


test_that("the sample columns for count output files from writeBambuOutput are numeric", {
    se <- readRDS(system.file("extdata", "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
    path <- test_path("fixtures")
    
    writeBambuOutput(se, path)
    
    outputFilesName <- list.files(path)
    countOutputFilesPath <- test_path("fixtures", outputFilesName[!grepl("extended_annotations.gtf", outputFilesName)])
    
    check <- sapply(seq_along(countOutputFilesPath), function(i){
        df <- read.table(countOutputFilesPath[i], header = TRUE)
        sampleColumn <- setdiff(colnames(df), c("GENEID", "TXNAME"))
        numericColumn <- colnames(select_if(df, is.numeric))
        all(sampleColumn %in% numericColumn)
    })
    
    expect_true(all(check))
    
    unlink(test_path("fixtures", "*"))
    
})


# test for the function: readGTF
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
