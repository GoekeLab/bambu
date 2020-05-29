## code to prepare `sysdata.rda` dataset goes here

a <- 500

data1 <-  data.frame(tx_id = c(1,2,3,1,2,1,2,3),
                                 read_class_id = c(1,2,3,rep(4,2),rep(5,3)),
                                 nobs = c(0,0,0,rep(10,2),rep(a,3)),
                                 gene_id = 1)

data2 <- data.frame(tx_id = c(1,2,1,2,3),
                                read_class_id = c(rep(4,2),rep(5,3)),
                                nobs = c(rep(10,2),rep(a,3)),
                                gene_id = 2)

data3 <- data.frame(tx_id = c(1,2),
                                read_class_id = c(1,2),
                                nobs = c(5,10),
                                gene_id = 3)

data4 <- data.frame(tx_id = c(1,2,3,1,2,1,2,3),
                                read_class_id = c(1,2,3,rep(4,2),rep(5,3)),
                                nobs = c(100,5,500,rep(20,2),rep(5,3)),
                                gene_id = 4)

data5 <- data.frame(tx_id = c(1:6,2,3,1,2,5,1,2,3,5,1:5),
                                read_class_id = c(1:6,rep(7,2),rep(8,3),rep(9,4),rep(10,5)),
                                nobs = c(2,0,1,0,0,1,rep(2,2),rep(282,3),rep(64,4),rep(5,5)),
                                gene_id = 5)



estOutput_woBC <- lapply(1:5, function(s){
  est <- bambu.quantDT(readClassDt = get(paste0("data",s)),  emParameters=list(bias = FALSE))
})

estOutput_wBC <- lapply(1:5, function(s){
  est <- bambu.quantDT(readClassDt = get(paste0("data",s)))
})



## expected output for test isore

seReadClass1 <- readRDS(system.file("extdata", "seReadClassUnstranded_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
extendedAnnotations <- readRDS(system.file("extdata", "extendedAnnotationGranges_txdbGrch38_91_chr9_1_1000000.rds", package = "bambu"))


seWithDistExpected <- isore.estimateDistanceToAnnotations(seReadClass=seReadClass1,
                                                  annotationGrangesList=extendedAnnotations,
                                                  min.exonDistance = 35)


## expected seGeneOutput

se <- readRDS(system.file("extdata", "seOutput_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
seExtended <- readRDS(system.file("extdata", "seOutputExtended_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
seCombined <- readRDS(system.file("extdata", "seOutputCombined_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))
seCombinedExtended <- readRDS(system.file("extdata", "seOutputCombinedExtended_SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.rds", package = "bambu"))


seGeneExpected <- transcriptToGeneExpression(se)
seExtendedGeneExpected <- transcriptToGeneExpression(seExtended)
seCombinedGeneExpected <- transcriptToGeneExpression(seCombined)
seCombinedExtendedGeneExpected <- transcriptToGeneExpression(seCombinedExtended)




usethis::use_data(data1,data2,data3,data4,data5,
                  estOutput_woBC,
                  estOutput_wBC,
                  standardJunctionModels_temp,
                  seWithDistExpected,
                  seGeneExpected,seExtendedGeneExpected,
                  seCombinedGeneExpected,seCombinedExtendedGeneExpected,
                  internal = TRUE,overwrite = TRUE)




