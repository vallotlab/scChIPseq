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
name                  - path to script location to set working directory
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
        
## Output

In the repo, the script should have created a directory **datasets** in which a new directory is created for each run with a different input name. Inside that directory are created a directory for each part of the analysis, containing RData and figures.
  
## Other

The config file **annotation/MSIGdb_classes** contains the MSIG predefined classes (one per line) used in the gene set enrichment step. You can modify this file to add or remove MSIG classes in your analysis. Check the MSIG db website :http://software.broadinstitute.org/gsea/msigdb .

The bash script **run.sh** contains the command lines used to produce analysis and most of the figures present in the paper. To run the analysis for the 4 datasets in the paper first download all the matrices and bam files in the repo root. Then run: 

```
cd <scChiPseq_SOURCE_DIRECTORY>
chmod 755 run.sh
./run.sh
```

# Authors
Please do not hesitate to post an issue or contact the authors :

Celine Vallot : celine.vallot@curie.fr

Pacome Prompsy : pacome.prompsy@curie.fr
