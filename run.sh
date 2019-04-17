#!/bin/bash

#Human 1
Rscript R_scChIP_seq_analysis.R . 'HBCx_95_human_order_1' 'hg38' test_set/HBCx_95_hg38/HBCx_95_CapaR_original_hg38.txt test_set/HBCx_95_hg38/HBCx_95_original_hg38.txt test_set/HBCx_95_hg38/HBCx_95_CapaR_flagged_rmDup.bam test_set/HBCx_95_hg38/HBCx_95_flagged_rmDup.bam -n 2 -p 2 -e annotation/hg38/exclude_regions_hg38.bed

#Human 2
Rscript R_scChIP_seq_analysis.R . 'HBCx_95_human_order_2' 'hg38' test_set/HBCx_95_hg38/HBCx_95_original_hg38.txt test_set/HBCx_95_hg38/HBCx_95_CapaR_original_hg38.txt  test_set/HBCx_95_hg38/HBCx_95_flagged_rmDup.bam test_set/HBCx_95_hg38/HBCx_95_CapaR_flagged_rmDup.bam -n 2 -p 2 -e annotation/hg38/exclude_regions_hg38.bed


#Mouse 1
Rscript R_scChIP_seq_analysis.R . 'HBCx_95_mouse_order_1' 'mm10' test_set/HBCx_95_mm10/HBCx_95_CapaR_original_mm10.txt test_set/HBCx_95_mm10/HBCx_95_original_mm10.txt  test_set/HBCx_95_mm10/HBCx_95_CapaR_flagged_rmDup.bam test_set/HBCx_95_mm10/HBCx_95_flagged_rmDup.bam -n 3 -p 1

#Mouse 2
Rscript R_scChIP_seq_analysis.R . 'HBCx_95_mouse_order_2' 'mm10' test_set/HBCx_95_mm10/HBCx_95_original_mm10.txt test_set/HBCx_95_mm10/HBCx_95_CapaR_original_mm10.txt test_set/HBCx_95_mm10/HBCx_95_flagged_rmDup.bam test_set/HBCx_95_mm10/HBCx_95_CapaR_flagged_rmDup.bam -n 3 -p 1
