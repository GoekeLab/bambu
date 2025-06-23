#download the bam file and make sure you set your own paths for the reference annotations
#change to bambu directory so that outputs are correct

#aws s3 cp --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/bam/genome/SGNex_HepG2_directRNA_replicate5_run1/SGNex_HepG2_directRNA_replicate5_run1.bam


devtools::load_all("bambu")
annotations = readRDS("Homo_sapiens.GRCh38.91.sorted.rds")
fa.file = "hg38_sequins_SIRV_ERCCs_longSIRVs.fa"
sample = "SGNex_HepG2_directRNA_replicate5_run1.bam"
  
#Get transcript discovery model
rcf <- bambu(reads = sample, annotations = annotations, genome = fa.file, discovery = FALSE, assignDist = FALSE, quant = FALSE, verbose = TRUE)

defaultModels = trainBambu(rcf[[1]])
xgb.save(defaultModels$transcriptModelME, "./inst/extdata/read_class_ME.model")
xgb.save(defaultModels$transcriptModelSE, "./inst/extdata/read_class_SE.model")
defaultModels$transcriptModelME = NULL
defaultModels$transcriptModelSE = NULL
saveRDS(defaultModels, "./inst/extdata/defaultModels.rds")


#Get junction model
readGrgList = prepareDataFromBam(sample, verbose = TRUE)
genomeSequence <- checkInputSequence(fa.file)
mcols(readGrgList)$id <- seq_along(readGrgList) 
unlisted_junctions <- unlistIntrons(readGrgList, use.ids = TRUE)
uniqueJunctions <- isore.constructJunctionTables(unlisted_junctions, 
                                                 annotations,genomeSequence, 
                                                 stranded = FALSE, verbose = TRUE,
                                                 returnModel = TRUE)
junctionModel = metadata(uniqueJunctions)$junctionModel
xgb.save(junctionModel$spliceSitePredictionStart.start, "./inst/extdata/spliceSitePredictionStart.start.model")
xgb.save(junctionModel$spliceSitePredictionStart.end, "./inst/extdata/spliceSitePredictionStart.end.model")
xgb.save(junctionModel$spliceSitePredictionEnd.start, "./inst/extdata/spliceSitePredictionEnd.start.model")
xgb.save(junctionModel$spliceSitePredictionEnd.end, "./inst/extdata/spliceSitePredictionEnd.end.model")