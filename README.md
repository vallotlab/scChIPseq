# Analysis of scChIPseq datasets

The scripts analyses and produce figures of high-throughput single-cell immuno-precitpitation followed by sequencing experiments.

## Running analysis 

### Download datasets from GEO

In a first time, download the dataset of interest from GEO (Grosselin et al.) https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117309. To run the script until the differential analysis step, you need count datamatrices from one cell line and 2 conditions (e.g. HBCx-95 human and HBCx-95 CapaR human). The peak calling and gene set enrichment parts require BAM files which are not accessible on GEO but can be privately provided by the authors (patient privacy protection). 

### Run script

First download the repository in the location of your choice, either with `git clone https://github.com/vallotlab/scChIPseq.git scChIPseq` or by clicking on 'Clone or Download' -> 'Download ZIP' and unzip.

Make sure to make the main script 'R_scChIP_seq_analysis' executable :

```
chmod 755 R_scChIP_seq_analysis.R
```

Run the script using Rscript with the following commands :

```
Rscript R_scChIP_seq_analysis.R  <source_file_directory> \
        <name> \
        <annot = mm10 | hg38> \
        <count_matrix_1.txt> \
        <count_matrix_2.txt> \
        -b1 <file_1.bam> \
        -b2 <file_2.bam> \
        -n <nclust> \
        -p <percent> \
        -e <exclude.bed> \
        -h ]
```

The arguments are described below : 

* Mandatory arguments (in the given order):


```
source_file_directory   - path to script location to set working directory
name                  - path to script location to set working directory
annot    - annotation to use ('mm10' or 'hg38')
count_matrix_1.txt   - full path to first count matrix file (.tsv/.txt)
count_matrix_2.txt   - full path to second count matrix file (.tsv/.txt)
```

* Optional arguments: 

```
-b1 file_1.bam          - full path to first bam file for peak calling (.bam)
-b2 file_2.bam          - full path to second bam file for peak calling (.bam)
-n nclust        - number of cluster to choose (optional)
-p percent [default = 1]         - percent (base 100) of cells to correlate with in correlation clustering and filtering step (optional) 
-e exclude.bed    -bed files containing regions to exclude (e.g. high CNV regions)
--help              - print this text
```
        
## Output

In the repo, the script should have created a directory 'datasets' in which a new directory is created for each run with a different input name. Inside that directory are created a directory for each part of the analysis, containing RData and figures.
  
## Other

The differential analysis and gene set enrichment analysis are quite sensitive and therefore changing even slightly the initial parameters might have a high impact on theses steps. Also, the correlation filtering step is based on a threshold calculated on a randomly sampled matrix, and is source of a tiny amount of randomness which might impact downstream steps. 

The bash script 'run.sh' contains the command lines used to produce analysis and most of the figures present in the paper. To run the analysis for the 4 datasets in the paper first download all the matrices and bam files in the repo root. Then run: 

```
cd <scChiPseq_SOURCE_DIRECTORY>
chmod 755 run.sh
./run.sh
```

# Authors
Please do not hesitate to post an issue or contact the authors :

Celine Vallot : celine.vallot@curie.fr

Pacome Prompsy : pacome.prompsy@curie.fr

Pia Kirchmeier

