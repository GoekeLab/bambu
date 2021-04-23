context('Get gene/transcript score and feature preparation')

library(SummarizedExperiment)
library(devtools)
library('BiocParallel')
library('Rsamtools')
library('GenomicAlignments')
library('dplyr')
library('ROCR')
library('stringr')
library('GenomicRanges')

load_all('~/Downloads/FYP/bambu_new/bambu')

test_that('getTranscriptScore returns expected dataframe of txScore and txFDR', {
  
  test.bam <- system.file("extdata",
                          "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam",
                          package = "bambu")
  fa.file <- system.file("extdata",
                         "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa",
                         package = "bambu")
  gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
  
  bambuAnnotations <- prepareAnnotations(gtf.file)
  
  bpParameters <- setBiocParallelParameters(test.bam, readClass.file = NULL,
                                            1, F)
  se_list <- bambu.processReads(test.bam, bambuAnnotations,
                           genomeSequence = fa.file,
                           yieldSize = NULL,
                           bpParameters = bpParameters, stranded = F, verbose =F)
  se <- se_list[[1]]

  # Filter to keep only read classes with > 1 read
  se = se[assays(se)$count>1,]
  
  txScore_df = getTranscriptScore(rowData(se))
  
  # Load in reference dataframe
  ref_txScore_df = read.csv('~/Downloads/FYP/bambu_new/bambu/inst/extdata/txScore_df.csv')

  #write.csv(txScore_df, '~/Downloads/FYP/bambu_new/bambu/inst/extdata/txScore_df.csv', row.names = F)
  
  expect_equal(txScore_df, ref_txScore_df)
})

test_that('getGeneScore returns expected dataframe of geneScore and geneFDR',{
  test.bam <- system.file("extdata",
                          "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam",
                          package = "bambu")
  fa.file <- system.file("extdata",
                         "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa",
                         package = "bambu")
  gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
  
  bambuAnnotations <- prepareAnnotations(gtf.file)
  
  bpParameters <- setBiocParallelParameters(test.bam, readClass.file = NULL,
                                            1, F)
  se_list <- bambu.processReads(test.bam, bambuAnnotations,
                                genomeSequence = fa.file,
                                yieldSize = NULL,
                                bpParameters = bpParameters, stranded = F, verbose =F)
  se <- se_list[[1]]
  
  # Filter to keep only read classes with > 1 read
  se = se[assays(se)$count>1,]
  
  geneScore_df = getGeneScore(rowData(se))
  
  # Load in reference dataframe
  ref_geneScore_df = read.csv('~/Downloads/FYP/bambu_new/bambu/inst/extdata/geneScore_df.csv')
  
  #write.csv(geneScore_df, '~/Downloads/FYP/bambu_new/bambu/inst/extdata/geneScore_df.csv', row.names = F)
  
  expect_equal(geneScore_df, ref_geneScore_df)
})