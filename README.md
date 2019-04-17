# Analysis of scChIPseq datasets

The scripts reproduce analysis and figures of single-cell immuno-precitpitation followed by sequencing experiments done in the paper : Grosselin et al., 2019. 


## Running analysis 

### Download datasets from GEO

In a first time, download the dataset of interest from GEO (Grosselin et al.) https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117309. To run the script until the differential analysis step, you need count datamatrices from one cell line and 2 conditions (e.g. HBCx-95 human and HBCx-95 CapaR human). The peak calling and gene set enrichment parts require BAM files which are not accessible on GEO but can be privately provided by the authors (patient privacy protection). 

### Run script

First download the repository in the location of your choice, either with ```git clone https://github.com/vallotlab/scChIPseq.git scChIPseq ``` or by clicking on 'Clone or Download' -> 'Download ZIP' and unzip.

Make sure to make the main script 'R_scChIP_seq_analysis' executable :
```chmod 755 R_scChIP_seq_analysis.R```

Run the script using Rscript with the following commands 
```Rscript R_scChIP_seq_analysis.R <source_file_directory> <name> <annot = mm10 | hg38> <count_matrix_1.txt> <count_matrix_2.txt> <file_1.bam> <file_2.bam>```
  
Example with GEO sample HBCx-95 and HBCx-95-CapaR mouse: 
```Rscript R_scChIP_seq_analysis.R '~/scChIPseq' 'HBCx_95_human_order_1' 'hg38' HBCx_95_CapaR_original_hg38.txt HBCx_95_original_hg38.txt HBCx_95_CapaR_flagged_rmDup.bam HBCx_95_flagged_rmDup.bam -n 2 -p 2 -e annotation/hg38/exclude_regions_hg38.bed
```

## Output
In the repo, the script should have created a directory 'datasets' in which a new directory is created for each run with a different input name. Inside that directory are created a directory for each part of the analysis, containing RData and figures 
  

# Authors
Please do not hesitate to post an issue or contact the authors :
Celine Vallot : celine.vallot@curie.fr

Pacome Prompsy : pacome.prompsy@curie.fr

Pia Kirchmeier

