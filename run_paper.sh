#!/bin/bash

#Move to source file directory root
#cd <SOURCE_FILE_DIRECTORY_ROOT>

#Mouse HBCx-95
echo -e 'hallmark\nc5_GO' > annotation/MSIGdb_classes
Rscript R_scChIP_seq_analysis.R . 'HBCx_95_mouse_paper' 'mm10' -1 MatCov_Sample1_merged123_mm10_rmdup_1000reads_mm10_50kb_chunks.txt -2 MatCov_Sample2_merged123_mm10_rmdup_1000reads_mm10_50kb_chunks.txt -b1 HBCx_95_mm10.bam -b2 HBCx_95_CapaR_mm10.bam -n 3 -p 1

#Human HBCx-22
echo -e 'hallmark\nc2_curated' > annotation/MSIGdb_classes
Rscript R_scChIP_seq_analysis.R . 'HBCx_22_human_paper' 'hg38' -1 MatCov_Sample4_merged_hg38_mapped_rmdup_1000reads_hg38_50kb_chunks.txt -2 MatCov_Sample5_merged_hg38_mapped_rmdup_1000reads_hg38_50kb_chunks.txt -b1 HBCx_22_hg38.bam -b2 HBCx_22_TamR_hg38.bam -n 2 -p 1

#Human HBCx-95
echo -e 'hallmark\nc2_curated' > annotation/MSIGdb_classes
Rscript R_scChIP_seq_analysis.R . 'HBCx_95_human_paper' 'hg38' -1 MatCov_Sample1_merged123_hg38_rmdup_1000reads_hg38_50kb_chunks.txt -2 MatCov_Sample2_merged123_hg38_rmdup_1000reads_hg38_50kb_chunks.txt -b1 HBCx_95_hg38.bam -b2 HBCx_95_CapaR_hg38.bam -n 2 -p 2 -e annotation/hg38/exclude_regions_hg38.bed


