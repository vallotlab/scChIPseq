# Analysis of scChIPseq datasets

Scripts to analyse high-throughput single-cell ChIP-seq experiments.

## How to run the analysis 

### Download datasets from GEO

Download the dataset of interest from GEO (Grosselin et al.) https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117309. To run the script up to the differential analysis step, you need the count matrices. The peak calling and gene set enrichment parts require BAM files. 

### Run script

First download the repository in the location of your choice, either with `git clone https://github.com/vallotlab/scChIPseq.git scChIPseq` or by clicking on 'Clone or Download' -> 'Download ZIP' and unzip.

Make sure to make the main script **R_scChIP_seq_analysis.R** executable :

```
chmod 755 R_scChIP_seq_analysis.R
```

Run the script using Rscript with the following commands :

```
Rscript R_scChIP_seq_analysis.R  <source_file_directory> \
        <name> \
        <annot = mm10 | hg38> \
        -1 <count_matrix_1.txt> \
        -2 <count_matrix_2.txt> \
        -b1 <file_1.bam> \
        -b2 <file_2.bam> \
        -n <nclust> \
        -p <percent> \
        -e <exclude.bed> \
        -h
```

The arguments are described below : 

* Mandatory arguments (in the given order):


```
source_file_directory   - path to script location to set working directory
name                  - Name of analysis
annot    - annotation to use ('mm10' or 'hg38')
-1 count_matrix_1.txt   - full path to first count matrix file (.tsv/.txt)
```

* Optional arguments: 

```
-[int] count_matrix_[int].txt   - full path to [int]th count matrix file (.tsv/.txt)
-b[int] file_[int].bam          - full path to [int]th bam file for peak calling (.bam)
-n nclust        - number of cluster to choose (optional)
-p percent [default = 1]         - percent (base 100) of cells to correlate with in correlation clustering and filtering step (optional) 
-e exclude.bed    -bed files containing regions to exclude (e.g. high CNV regions)
--help              - print this text
```
## Requirements
```
  #Bioinfo
  library(scater)
  library(scran)
  library(IRanges)
  library(GenomicRanges)
  library(ConsensusClusterPlus)
  library(Rtsne)
  
  #Data mining & utils
  library(tibble)
  library(dplyr)
  library(stringr)
  library(irlba)
  library(reshape2)
  library(DT)
  library(tidyr)
  library(splitstackshape)
  library(rlist)
  library(envDocument)
  library(rstudioapi)
  library(dplyr)
  
  #geco
  library(geco.utils)
  library(geco.visu)
  library(geco.unsupervised)
  library(geco.supervised)
  
  #Graphics
  library(RColorBrewer)
  library(colorRamps)
  library(colourpicker)
  library(kableExtra)
  library(knitr)
  library(viridis)
  library(ggplot2)
  library(gplots)
  library(png)
  library(grid)
  library(gridExtra)

  #Modules and functions
  source("Modules/geco.annotToCol2.R")
  source("Modules/geco.wilcox.R")
```



Note that the geco packages are custom packages that are embedded in the application (under packages/).
If the **installation script** didn't work, you can try to install them manually :
```
install.packages("packages/geco.utils.tar.gz",repos = NULL,type = "source")
install.packages("packages/geco.visu.tar.gz",repos = NULL,type = "source")
install.packages("packages/geco.unsupervised.tar.gz",repos = NULL,type = "source")
install.packages("packages/geco.supervised.tar.gz",repos = NULL,type = "source")
```


## Output

In the repo, the script should have created a directory **datasets** in which a new directory is created for each run with a different input name. Inside that directory are created a directory for each part of the analysis, containing RData and figures.
  
## Other

The config file **annotation/MSIGdb_classes** contains the MSIG predefined classes (one per line) used in the gene set enrichment step. You can modify this file to add or remove MSIG classes in your analysis. Check the MSIG db website :http://software.broadinstitute.org/gsea/msigdb .

The bash script **run_paper.sh** contains the command lines used to produce analysis and most of the figures present in the paper. To run the analysis for the 4 datasets in the paper first download all the matrices and bam files in the repo root. Then run: 

```
cd <scChiPseq_SOURCE_DIRECTORY>
chmod 755 run_paper.sh
./run_paper.sh
```

# Authors
Please do not hesitate to post an issue or contact the authors :

Celine Vallot : celine.vallot@curie.fr

Pacome Prompsy : pacome.prompsy@curie.fr


# Session Info
```
R version 3.5.2 (2018-12-20)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.6 LTS

Matrix products: default
BLAS: /usr/lib/libblas/libblas.so.3.6.0
LAPACK: /usr/lib/lapack/liblapack.so.3.6.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=fr_FR.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=fr_FR.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] grid      stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] gridExtra_2.3                     png_0.1-7                         gplots_3.0.1                     
 [4] viridis_0.5.1                     viridisLite_0.3.0                 knitr_1.22                       
 [7] kableExtra_1.1.0                  colourpicker_1.0                  colorRamps_2.3                   
[10] RColorBrewer_1.1-2                geco.supervised_1.0.0             geco.unsupervised_1.0.0          
[13] geco.visu_1.0.0                   geco.utils_1.0.0                  rstudioapi_0.10                  
[16] envDocument_2.4.0                 rlist_0.4.6.1                     splitstackshape_1.4.6            
[19] tidyr_0.8.3                       DT_0.5                            reshape2_1.4.3                   
[22] irlba_2.3.3                       Matrix_1.2-15                     stringr_1.4.0                    
[25] dplyr_0.8.0.1                     tibble_2.1.1                      edgeR_3.24.3                     
[28] limma_3.38.3                      Rtsne_0.15                        ConsensusClusterPlus_1.46.0      
[31] scran_1.10.2                      scater_1.10.1                     ggplot2_3.1.0                    
[34] SingleCellExperiment_1.4.1        SummarizedExperiment_1.12.0       DelayedArray_0.8.0               
[37] BiocParallel_1.16.6               matrixStats_0.54.0                Biobase_2.42.0                   
[40] seqCBS_1.2.1                      clue_0.3-57                       DNAcopy_1.56.0                   
[43] SCOPE_0.0.1                       CODEX2_1.3.0                      BSgenome.Hsapiens.UCSC.hg19_1.4.0
[46] Rsamtools_1.34.1                  BSgenome.Hsapiens.UCSC.hg38_1.4.1 BSgenome_1.50.0                  
[49] rtracklayer_1.42.2                Biostrings_2.50.2                 XVector_0.22.0                   
[52] GenomicRanges_1.34.0              GenomeInfoDb_1.18.2               IRanges_2.16.0                   
[55] S4Vectors_0.20.1                  BiocGenerics_0.28.0              

loaded via a namespace (and not attached):
 [1] ggbeeswarm_0.6.0         colorspace_1.4-1         dynamicTreeCut_1.63-1    BiocNeighbors_1.0.0      manipulate_1.0.1        
 [6] mvtnorm_1.0-10           xml2_1.2.0               cluster_2.0.7-1          shiny_1.2.0              HDF5Array_1.10.1        
[11] httr_1.4.0               readr_1.3.1              BiocManager_1.30.4       compiler_3.5.2           assertthat_0.2.1        
[16] lazyeval_0.2.2           later_0.8.0              htmltools_0.3.6          tools_3.5.2              igraph_1.2.4            
[21] gtable_0.3.0             glue_1.3.1               GenomeInfoDbData_1.2.0   Rcpp_1.0.1               gdata_2.18.0            
[26] DelayedMatrixStats_1.4.0 xfun_0.6                 rvest_0.3.2              mime_0.6                 miniUI_0.1.1.1          
[31] gtools_3.8.1             statmod_1.4.30           XML_3.98-1.19            zlibbioc_1.28.0          MASS_7.3-51.1           
[36] scales_1.0.0             hms_0.4.2                promises_1.0.1           rhdf5_2.26.2             expm_0.999-4            
[41] yaml_2.2.0               stringi_1.4.3            caTools_1.17.1.2         boot_1.3-20              rlang_0.3.3             
[46] pkgconfig_2.0.2          bitops_1.0-6             evaluate_0.13            lattice_0.20-38          purrr_0.3.2             
[51] Rhdf5lib_1.4.3           GenomicAlignments_1.18.1 htmlwidgets_1.3          tidyselect_0.2.5         plyr_1.8.4              
[56] magrittr_1.5             R6_2.4.0                 DescTools_0.99.28        pillar_1.3.1             foreign_0.8-71          
[61] withr_2.1.2              RCurl_1.95-4.12          crayon_1.3.4             KernSmooth_2.23-15       rmarkdown_1.12          
[66] locfit_1.5-9.1           data.table_1.12.0        webshot_0.5.1            digest_0.6.18            xtable_1.8-3            
[71] httpuv_1.5.0             munsell_0.5.0            beeswarm_0.2.3           vipor_0.4.5
```
