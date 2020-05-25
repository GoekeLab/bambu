<a href="https://raw.githubusercontent.com/GoekeLab/bambu/master/figures/transparent-bambu.png?token=AGA7DTCQ2VT5ILG3R6ORKUK6WP424"><img src="https://raw.githubusercontent.com/GoekeLab/bambu/master/figures/transparent-bambu.png?token=AGA7DTCQ2VT5ILG3R6ORKUK6WP424" title="Bambu" alt="Bambu"></a>

# bambu: reference-guided transcript discovery and quantification for long read RNA-Seq data


***bambu*** is a R package for multi-sample transcript discovery and quantification using long read RNA-Seq data. You can use ***bambu*** after read alignment to obtain expression estimates for known and novel transcripts and genes. The output from ***bambu*** can directly be used for visualisation and downstream analysis such as differential gene expression or transcript usage.



[![bambu](https://img.shields.io/badge/bambu-v0.9.0-brightgreen)](https://github.com/GoekeLab/bambu) [![Maintained?](https://img.shields.io/badge/Maintained%3F-Yes-blue)](https://gemnasium.com/badges/badgerbadgerbadger)  [![Install](https://img.shields.io/badge/Install-Github-brightgreen)](https://github.com/badges/badgerbadgerbadger/issues) 
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

---


## Installation

You can install ***bambu*** from github:

```rscript
install.packages("devtools")
devtools::install_github("GoekeLab/bambu")
```
---

## General Usage 

The default mode to run ***bambu*** is using a set of aligned reads (bam files), reference genome annotations (gtf file, TxDb object, or bambuAnnotation object), and reference genome sequence (fasta file or BSgenome). ***bambu*** will return a summarizedExperiment object with the genomic coordinates for annotated and new transcripts and genes and with the transcript expression estimates: 
 ```rscript
 
bambu(reads, annotations, genomeSequence,...)

library(bambu)

test.bam <- system.file("extdata", "SGNex_HepG2_cDNAStranded_replicate5_run4_chr9_108865774_109014097.bam", package = "bambu")
  

se <- bambu(reads = test.bam, annotations = "TxDb.Hsapiens.UCSC.hg38.knownGene", genomeSequence = "BSgenome.Hsapiens.NCBI.GRCh38")
       
```


We highly recommend to use the same annotations that were used for genome alignment. If you have a gtf file and fasta file you can run ***bambu*** with the following options:

```rscript
test.bam <- system.file("extdata", "SGNex_HepG2_cDNAStranded_replicate5_run4_chr9_108865774_109014097.bam", package = "bambu")
  
fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9.fa.gz", package = "bambu")

gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_108865774_109014097.gtf", package = "bambu")

bambuAnnotations <- prepareAnnotationsFromGtf(gtf.file)

se <- bambu(reads = test.bam, annotations=bambuAnnotations, genomeSequence ="BSgenome.Hsapiens.NCBI.GRCh38")

```
---


## Other Use Cases
**Use precalculated annotation objects**

If you plan to run ***bambu*** more frequently, you can store a bambuAnnotations object.

```rscript
library(bambu)
test.bam <- system.file("extdata", "SGNex_HepG2_cDNAStranded_replicate5_run4_chr9_108865774_109014097.bam", package = "bambu")

gr <- readRDS(system.file("extdata", "annotationGranges_txdbGrch38_91_chr9_108865774_109014097.rds", package = "bambu"))

# or hg38 <- readRDS(url(‘...’))

se <- bambu(reads = test.bam, annotations = gr, genomeSequence = fa.file)

```

The bambuAnnotation object can be calculated from a gtf file:
```rscript
annotations <- prepareAnnotationFromGTF(gtf.file)
```

From TxDb object
```rscript
annotations <- prepareAnnotations(txdb)
```

You can also download annotations for recent genome releases here:
-  hg38:  ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens          
-  mm10:  ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus     

**Quantification of annotated transcripts and genes only (no transcript/gene discovery)**

```rscript
bambu(reads = test.bam, annotations = txdb, genomeSequence = fa.file, extendAnnotations = FALSE)
```

**Large sample number/ limited memory**     
For larger sample numbers we recommend to write the processed data to a file:
```rscript
bambu(reads, outputReadClassFolder, genomeSequence, annotations)
```



---

## Advanced Options

**More strigent filtering thresholds imposed on potential novel transcripts**    
- For example   
> Keep novel transcripts with min 5 read count in at least 1 sample:  
```rscript
bambu(reads, genomeSequence, annotations, ir.control = list(min.readCount = 5))
```

> Keep novel transcripts with min 5 samples having at least 2 counts:

```rscript
bambu(reads, genomeSequence, annotations, ir.control = list(min.sampleNumber = 5))
```

> Filter out transcripts with relative abundance within gene lower than 10%: 
```rscript
bambu(reads, genomeSequence, annotations, ir.control = list(min.readFractionByGene = 0.1))
```

**Quantification without bias correction**     
> The default estimation automatically does bias correction for expression estimates. However, you can choose to perform the quantification without bias correction.    
```rscript
bambu(reads, genomeSequence, annotations, algo.control(bias_correction = FALSE))
```

**Parallel computation**      
> ***bambu*** also allows parallel computation for EM.    
```rscript
bambu(reads, genomeSequence, annotations, ncore = 8)
```

See [manual]() for details to customize other conditions.

---

## Complementary functions

**Transcript expression to gene expression**

```rscript
transcriptToGeneExpression(se)
```

**Visualization**
> You can visualize the novel genes/transcripts using ***plot.bambu*** function 
```rscript
plot.bambu(se, type = "annotation", gene_id)

plot.bambu(se, type = "annotation", transcript_id)
```

> ***plot.bambu*** can also be used to visualize the clustering of input samples on gene/transcript expressions
```rscript
plot.bambu(se, type = "heatmap") # heatmap 

plot.bambu(se, type = "pca") # PCA visualization
```

> ***plot.bambu*** can also be used to visualize the clustering of input samples on gene/transcript expressions with grouping variable
```rscript
plot.bambu(se, type = "heatmap", group.var) # heatmap 

plot.bambu(se, type = "pca", group.var) # PCA visualization
```


---


## Contributors

This package is developed and maintained by [Jonathan Goeke](https://github.com/jonathangoeke) and [Ying Chen](https://github.com/cying111) at Genome Institute of Singapore. If you want to contribute, please leave an issue. Thank you.

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

- **[GPL v3](https://www.gnu.org/licenses/gpl-3.0)**
