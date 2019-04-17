# Analysis of scChIP-seq datasets

From count matrices (50kb/5kbp regions x cells) to consensus clustering, differential analysis and gene set enrichment. 

Use the script : 
In your prefered command line run : 

```chmod 755 R_scChIP_seq_analysis.R```

```Rscript R_scChIP_seq_analysis.R <source_file_directory> <name> <annot = mm10 | hg38> <count_matrix_1.txt> <count_matrix_2.txt> <file_1.bam> <file_2.bam>```
  
Datasets : load the datasets of your choice from Grosselin et al. - downloadable from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117309

e.g. 
```Rscript R_scChIP_seq_analysis.R . 'HBCx_95_mouse' 'mm10' HBCx_95_original_mm10.txt HBCx_95_CapaR_original_mm10.txt HBCx_95_flagged.bam HBCx_95_CapaR_flagged_rmDup.bam 3```
  

# Authors
Celine Vallot : celine.vallot@curie.fr

Pacome Prompsy : pacome.prompsy@curie.fr

Pia Kirchmeier

