# Single-Cell Chromatin Immuno-precipitation followed by sequencing
# Testing multiple minimum read coverage per cell and plotting pca 
# Authors:
# Celine Vallot : celine.vallot@curie.fr
# Pacome Prompsy : pacome.prompsy@curie.fr
# Pia Kirchmeier

#Read in user input commands 
args <- commandArgs(TRUE)
input = list()
## Default setting when no arguments passed
print(length(args))
if(length(args) < 7) {
  args <- c("--help")
} else{
  input$source_file_directory = as.character(args[1])
  input$name = as.character(args[2])
  input$count_matrix_1 = as.character(args[3])
  input$count_matrix_2 = as.character(args[4])
  input$from = as.integer(args[5])
  input$to = as.integer(args[6])
  input$by = as.integer(args[7])
  
  if( is.null(input$from) |  is.null(input$to)  | is.null(input$by)) {
    print("ERROR : Enter correct arguments for 'from', 'to' and 'by': e.g. : 1000 2000 200")
    args <- c("--help")
  } else if( (input$from >  input$to) | (input$by >  input$to) | ( input$to - input$from  < input$by) | ( (input$to - input$from) %% input$by != 0)) {
    print("ERROR : Enter correct arguments for 'from', 'to' and 'by': e.g. : 1000 2000 200 ")
    args <- c("--help")
  }
  if(!file.exists(input$count_matrix_1) | !file.exists(input$count_matrix_2) ){
    print("ERROR : Count Matrix file doesn't exists")
    args <- c("--help")
  }
}
if("--help" %in% args | "-h" %in% args) {
  cat("
      Single-cell ChIP seq analysis :
      
      Arguments:
      source_file_directory   - path to script location to set working directory
      name   - path to script location to set working directory
      count_matrix_1.txt   - full path to first count matrix file (.tsv/.txt)
      count_matrix_2.txt   - full path to second count matrix file (.tsv/.txt)
      from   - number of reads to start with for the minimum coverage per cell range
      to   - number of reads to end with for the minimum coverage per cell range
      by   - increment for the minimum coverage per cell range

      Example:
      Rscript R_scChIP_multiple_pca.R '~/scChIPseq/' 'HBCx_95'  HBCx_95_CapaR_original_mm10.txt  HBCx_95_original_mm10.txt 1000 2000 200 \n\n
      ")
  
  q(save="no")
}
cat("source_file_directory = ", input$source_file_directory, "\n") 
cat("name = ",input$name, "\n") 
cat("annotation_id = ",input$annotation_id, "\n")
cat("count_matrix_1 = ",input$count_matrix_1, "\n")
cat("count_matrix_2 = ",input$count_matrix_2, "\n" )
cat("Testing filtering for minimum coverages per cell = ", seq(from = input$from, to = input$to, by = input$by)  ,"\n")

####################################################################################
#If running from R : 
# input = list()
# input$source_file_directory = "/home/pprompsy/Documents/GitLab/single_cell_ChIPseq_github/"
# input$name = "HBCx_95_hg38_2"
# input$annotation_id = "hg38"
# input$count_matrix_1 = "test_set/HBCx_95_hg38/HBCx_95_CapaR_original_hg38.txt"
# input$count_matrix_2 = "test_set/HBCx_95_hg38/HBCx_95_original_hg38.txt"
# input$from =  1000
# input$to =  2000
# input$by = 200

####################################################################################
#Set working directory to source file location (need change :)
setwd(input$source_file_directory)

#User input parameters :

input$datafile_matrix = list(datapath=c(input$count_matrix_1,input$count_matrix_2),
                             name = c( sub(".txt|.tsv","",basename(input$count_matrix_1)),sub(".txt|.tsv","",basename(input$count_matrix_2))))

inputBams <- list(input$bam1,input$bam2)
print(input$datafile_matrix )
print(inputBams)
clust <- list(cc.col=NULL, consclust.mat=NULL, hc=NULL, tsne_corr=NULL, annot_sel2=NULL,
              clust_pdf=NULL, available_k=10, chi=NULL)
cf=list()

#Parameters

#QC & Filtering
input$min_coverage_cell = 1600
input$min_cells_window = 1
input$quant_removal = 95
#Correlation Clustering Filtering
input$corr_thresh=98
if(is.null(input$percent_corr )) input$percent_corr = 1

#Diff Analysis 
input$cdiff.th = 1
input$qval.th = 0.01

#0. Importing packages ####
#Bioinfo
library(scater)
library(scran)
library(IRanges)
library(GenomicRanges)
library(ConsensusClusterPlus)
library(Rtsne)
library(edgeR)

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


#Custom functions
gg_fill_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

choose_perplexity <- function(dataset){
  perplexity=30
  if(nrow(dataset) <= 200 ){
    perplexity = 20
  }
  if(nrow(dataset) <= 250 ){
    perplexity = 25
  }
  if (nrow(dataset) <= 150 ){
    perplexity= 15
  }
  if (nrow(dataset) <= 100 ){
    perplexity= 10
  } 
  if (nrow(dataset) <= 50 ) {
    perplexity= 5
  }
  perplexity
}

#Inititalization step
init <- list(data_folder=NULL, datamatrix=data.frame(), annot_raw=data.frame(), available_raw_datasets = NULL,
             available_reduced_datasets=NULL, available_filtered_datasets=NULL)

annotCol <- c("sample_id","total_counts")

init$data_folder = "."

print("Creating folders")
#Creating Dir
dir.create(file.path(init$data_folder,"datasets"), showWarnings = FALSE)
dir.create(file.path(init$data_folder,"datasets", input$name), showWarnings = FALSE)
dir.create(file.path(init$data_folder,"datasets", input$name, "reduced_data"), showWarnings = FALSE)
dir.create(file.path(init$data_folder,"datasets", input$name, "cor_filtered_data"), showWarnings = FALSE)
dir.create(file.path(init$data_folder,"datasets", input$name, "consclust"), showWarnings = FALSE)
dir.create(file.path(init$data_folder,"datasets", input$name, "supervised"), showWarnings = FALSE)

datamatrix=NULL
annot_raw = NULL

print("Loading count matrices")
for(i in 1:length(input$datafile_matrix$datapath)){
  datamatrix_single <- read.table(input$datafile_matrix$datapath[i], header=TRUE, stringsAsFactors=FALSE)
  
  datamatrix_single <- datamatrix_single[!duplicated(rownames(datamatrix_single)),] #put IN for new format
  rownames(datamatrix_single) <- gsub(":", "_", rownames(datamatrix_single))
  rownames(datamatrix_single) <- gsub("-", "_", rownames(datamatrix_single))
  total_cell <-length(datamatrix_single[1,])
  sample_name <- gsub('.{4}$', '', input$datafile_matrix$name[i])
  annot_single <- data.frame(barcode=colnames(datamatrix_single), cell_id=paste0(sample_name, "_c", 1:total_cell), sample_id=rep(sample_name, total_cell), batch_id=i)
  colnames(datamatrix_single) <- annot_single$cell_id
  if(is.null(datamatrix)){ datamatrix <- datamatrix_single
  }else{
    common_regions <- intersect(rownames(datamatrix), rownames(datamatrix_single))
    datamatrix <- cbind(datamatrix[common_regions,], datamatrix_single[common_regions,])
  }
  if(is.null(annot_raw)){ annot_raw <- annot_single} else{ annot_raw <- rbind(annot_raw, annot_single)}
}

save(datamatrix, annot_raw, file=file.path(init$data_folder, "datasets", input$name, "scChIP_raw.RData"))
datamatrix <- datamatrix
annot_raw <-annot_raw

#Filtering and QC
print("Filtering and QC")
exclude_regions = NULL
###############################################################
# 1. Data loading
###############################################################
assays = list(counts = as.matrix(datamatrix)) #, colData =  annot_raw

umi <- SingleCellExperiment(assays = assays, colData = annot_raw)
#umi <- umi[rowSums(counts(umi)>0)>0, ] # remove windows that do not have any read in any cells
umi <- scater::calculateQCMetrics(umi)
thresh <- quantile(colSums(counts(umi)), probs=seq(0, 1, 0.01))

###############################################################
# 2. Filtering & Window selection
###############################################################
#Cell selection based on number of total counts, between min_cov_cell and upper 5%
annot <- colData(umi)




sel1000 =  (colSums(counts(umi))>1000 & colSums(counts(umi))< thresh[input$quant_removal+1])
SelMatCov1000 <- counts(umi)[,sel1000]
bina_counts <- SelMatCov1000
bina_counts[bina_counts<2] <-0
bina_counts[bina_counts>1] <-1
fixedWin <- names(which((rowSums(bina_counts) > ( (input$min_cells_window/100.0)*(dim(bina_counts)[2])) ))) # window selection


Nmin <- c(1000,1200,1400,1600,1800,2000) # Minimum number of total_counts

for(NbCounts in Nmin){
  print("Running filtering and pca for n min =")
  print(NbCounts)

  annot <- colData(umi)

  #Cell selection based on number of total counts, between NbCounts and 10000 (to remove outliers)
  sel <- ( colSums(counts(umi))> NbCounts & colSums(counts(umi)) < thresh[input$quant_removal+1] )

  SelMatCov <- counts(umi)[,sel]

  annot <- annot[sel,]
  SelMatCov <- SelMatCov[fixedWin,]

  print("Filtered cells :")
  print(dim(SelMatCov))

  #Normalization, NormMatCov
  mat <-  mean(colSums(SelMatCov))*t(t(SelMatCov)/colSums(SelMatCov))
  colnames(mat) <- annot$cell_id

  #Result file
  resSUBdir <- file.path(init$data_folder, "datasets", input$name,"reduced_data",paste0("mouse_PCA_based_fixedWin_NbMinCounts",NbCounts)) ; if(!file.exists(resSUBdir)){dir.create(resSUBdir)}
  anocol <- geco.annotToCol2(annotS=annot[,annotCol],annotT=annot,plotLegend=T,plotLegendFile=file.path(resSUBdir,"Annotation_legends.pdf"), categCol=NULL)

  mat <- mat-apply(mat,1,mean)

  ################################################################################
  ## PCA
  ################################################################################
  print("Running pca ...")
  pca <- stats::prcomp(t(mat),center=F,scale.=F)
  print("PCA done !")


  ## Automatic 2D PCA plots all versus all (up to the PC number specified)
  significant.PCs <- 2
  combinations <- as.data.frame(combn(significant.PCs,2,simplify=TRUE))
  extendXlimLeft <- 5 # Extend left Xlim by x percent of min
  extendXlimRight <- 15 # Extend right Xlim by x percent of max
  extendYlimBottom <- 5 # Extend bottom Ylim by x percent of min
  extendYlimTop <- 5 # Extend top Ylim by x percent of max
  TextSize <- 0.4
  pcaText <- FALSE
  annotText <- "sample_id"
  hcText <- "sample_id"

    for (i in 1:ncol(combinations)) {
    pdf(file.path(resSUBdir,sprintf("PCA%sVs%s_plot_2D.pdf", combinations[1,i], combinations[2,i] )), height=5, width=5)
    for(j in 1:ncol(anocol))	{
      xlimits <- c((min(pca$x[,combinations[1,i]])-abs(min(pca$x[,combinations[1,i]])/100*extendXlimLeft)), (max(pca$x[,combinations[1,i]])+abs(max(pca$x[,combinations[1,i]])/100*extendXlimRight)))
      ylimits <- c((min(pca$x[,combinations[2,i]])-abs(min(pca$x[,combinations[2,i]])/100*extendYlimBottom)), (max(pca$x[,combinations[2,i]])+abs(max(pca$x[,combinations[2,i]])/100*extendYlimTop)))
      plot(pca$x[,combinations[1,i]],pca$x[,combinations[2,i]],col=alpha(anocol[,j],0.6),xlab=colnames(pca$x)[combinations[1,i]],ylab=colnames(pca$x)[combinations[2,i]],cex=0.8,lwd=1.5, main=colnames(anocol)[j], pch=19,
           xlim=xlimits,
           ylim=ylimits
      ) # changed type of point and xlim
      if(pcaText)	text(pca$x[,combinations[1,i]],pca$x[,combinations[2,i]],labels=annot[,annotText],pos=4,cex=TextSize)
    }
    dev.off()

  }
}