# Analysis of scChIPseq datasets

From count matrices (50kb/5kbp regions x cells) to consensus clustering, differential analysis and gene set enrichment. 

Use the script : 
In your prefered command line run : 

```chmod 755 R_scChIP_seq_analysis.R```

```Rscript R_scChIP_seq_analysis.R <source_file_directory> <name> <annot = mm10 | hg38> <count_matrix_1.txt> <count_matrix_2.txt> <file_1.bam> <file_2.bam>```
  
Datasets : load the datasets of your choice from Grosselin et al. - downloadable from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117309
 
Example with GEO sample HBCx-95 and HBCx-95-CapaR mouse: 
```Rscript R_scChIP_seq_analysis.R '~/scChIPseq' 'HBCx_95_human_order_1' 'hg38' HBCx_95_CapaR_original_hg38.txt HBCx_95_original_hg38.txt HBCx_95_CapaR_flagged_rmDup.bam HBCx_95_flagged_rmDup.bam -n 2 -p 2 -e annotation/hg38/exclude_regions_hg38.bed
```
  

# Authors
Celine Vallot : celine.vallot@curie.fr

Pacome Prompsy : pacome.prompsy@curie.fr

Pia Kirchmeier

