# Single-Cell Chromatin Immuno-precipitation followed by sequencing
# analysis : from count matrices (50kb regions x cells ) to consensus correlation clustering,
# differential analysis, gene set enrichment 
# Dataset : load the datasets of your choice - downloadable from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117309
# e.g. : HBCx_95 and HBCx_95_CapaR 
# Small differences compared to the paper can be due to the randomness in correlation filtering step
# Also, the class names (C1,C2,...) are assigned in random order which might not be the same from paper
# Authors:
# Celine Vallot : celine.vallot@curie.fr
# Pacome Prompsy : pacome.prompsy@curie.fr
# Pia Kirchmeier


##############################################################################################################################
### 0. Initialization
##############################################################################################################################

print("Initializing pipeline...")

## 0.1 Read in user input commands ## 

  args <- commandArgs(TRUE)
  input = list()
  print(args)  
  if(length(args) < 5) {
    args <- c("--help")
  } else{
    input$source_file_directory = as.character(args[1])
    input$name = as.character(args[2])
    input$annotation_id = as.character(args[3])

    for(i in 1:10){
      if( paste0('-',i) %in% args) {
        eval(parse(text = paste0('input$count_matrix_',i,' = as.character(args[which(args == paste0("-",',i,'))+1])' )))
        if(!file.exists(eval(parse(text = paste0('input$count_matrix_',i))))){
          print("ERROR :  Count Matrix file not found ")
          q(save="no")
        }
      }
      if( paste0('-b',i) %in% args) {
        eval(parse(text = paste0('input$bam',i,' = as.character(args[which(args == paste0("-b",',i,'))+1])' )))
        if(!file.exists(eval(parse(text = paste0('input$bam',i))))){
          print("ERROR :  Bam file not found ")
          q(save="no")
        }
      }
    }
    
    if( '-n' %in% args) {
      input$nclust = as.integer(args[which(args == '-n')+1])
    }
    if( '-p' %in% args) {
      input$percent_corr = as.integer(args[which(args == '-p')+1])
    }
    if ('-e' %in% args) {
      input$exclude = as.character(args[which(args == '-e')+1])
    }
    if(input$annotation_id != "mm10" && input$annotation_id != "hg38"){
      print("ERROR :  Annotation id in wrong format, please input 'mm10' or 'hg38' ")
      args <- c("--help")
    }
    
  }

  #If running from R, change and uncomment below  
  # input = list()
  # input$source_file_directory = "/home/pprompsy/Documents/GitLab/scChIPseq/"
  # input$name = "HBCx_95_human_paper"
  # input$annotation_id = "hg38"
  # input$count_matrix_1 = "test_set/HBCx_95_hg38_paper/MatCov_Sample1_merged123_hg38_rmdup_1000reads_hg38_50kb_chunks.txt"
  # input$count_matrix_2 = "test_set/HBCx_95_hg38_paper/MatCov_Sample2_merged123_hg38_rmdup_1000reads_hg38_50kb_chunks.txt"
  # input$bam1 =  "test_set/HBCx_95_hg38_paper/HBCx_95_flagged_rmDup.bam"
  # input$bam2 =  "test_set/HBCx_95_hg38_paper/HBCx_95_flagged_rmDup.bam"
  # input$nclust = 2
  # input$percent_corr = 2
  # input$exclude = "annotation/hg38/exclude_regions_hg38.bed"
  
   if("--help" %in% args | "-h" %in% args) {
    cat("
        Single-cell ChIP seq analysis :
   
        >Mandatory arguments (in the given order):

        source_file_directory   - path to script location to set working directory
        name                  - path to script location to set working directory
        annot= <mm10|hg38>   - annotation to use (Mouse or Human)
        -1 count_matrix_1.txt   - full path to first count matrix file (.tsv/.txt)

        >Optional arguments: 

        -[int] count_matrix_[int].txt   - full path to [int]th count matrix file (.tsv/.txt)
        -b[int] file_[int].bam          - full path to [int]th bam file for peak calling (.bam)
        -n nclust         - number of cluster to choose (optional)
        -p percent [default = 1]         - percent (base 100) of cells to correlate with in correlation clustering and filtering step (optional) 
        -e exclude.bed    -bed files containing regions to exclude (e.g. high CNV regions)
        --help              - print this text
   
        Example:
        Rscript R_scChIP_seq_analysis.R . 'HBCx_95' 'hg38'  -1 HBCx_95_hg38.txt -2 HBCx_95_CapaR_hg38.txt -b1 HBCx_95_flagged.bam -b2 HBCx_95_CapaR_hg38.bam -n 3 -p 1 -e regions.bed  \n\n
        ")
    
    q(save="no")
  }

## 0.2 Initializing fixed parameters ## 

  #Set working directory to source file location (need change :)
  setwd(input$source_file_directory)
  
  #User input parameters :
  
  input$datafile_matrix = list(datapath=c(input$count_matrix_1,input$count_matrix_2),
                               name = c( basename(input$count_matrix_1),basename(input$count_matrix_2)))
  
  inputBams <- list(input$bam1,input$bam2)
  print(input$datafile_matrix )
  print(inputBams)
  clust <- list(cc.col=NULL, consclust.mat=NULL, hc=NULL, tsne_corr=NULL, annot_sel2=NULL,
                                        clust_pdf=NULL, available_k=10, chi=NULL)
  cf=list()
  
  #QC & Filtering
  input$min_coverage_cell = 1600
  input$min_cells_window = 1
  input$quant_removal = 95

  
  #Correlation Clustering Filtering
  input$corr_thresh=99
  if(is.null(input$percent_corr )) input$percent_corr = 1
  
  #Diff Analysis 
  input$cdiff.th = 1
  input$qval.th = 0.01

  
  init <- list(data_folder=NULL, datamatrix=data.frame(), annot_raw=data.frame(), available_raw_datasets = NULL,
               available_reduced_datasets=NULL, available_filtered_datasets=NULL)
  
  annotCol <- c("sample_id","total_counts")
  
  init$data_folder = "."
  
  
## 0.3 Importing packages ##

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


## 0.4 Custom functions ##

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

## 0.5 Loading all previous data for trouble shooting ##  
  
  #Uncomment specific data loading for trouble shooting :
  
  # load( file=file.path(init$data_folder, "datasets", input$name, "scChIP_raw.RData"))
  # load(file.path(init$data_folder, "datasets", input$name, "reduced_data", paste(input$name, input$min_coverage_cell,input$min_cells_window, input$quant_removal, "uncorrected", sep="_", "normMat.RData")))
  # load(file=file.path(init$data_folder, "datasets", input$name, "reduced_data", paste0(paste(input$name, input$min_coverage_cell, input$min_cells_window, input$quant_removal, "uncorrected", sep="_"), "_annotFeat.RData")))
  # load(file.path(init$data_folder, "datasets", input$name, "reduced_data", paste0(paste(input$name, input$min_coverage_cell, input$min_cells_window, input$quant_removal, "uncorrected", sep="_"), ".RData")))
  # load(file=file.path(init$data_folder, "datasets", input$name, "cor_filtered_data", paste0(input$name, "_", input$corr_thresh, "_", input$percent_corr, ".RData")))
  # load(file=file.path(init$data_folder, "datasets", input$name, "consclust", paste0(input$name, ".RData")))
  # load(file=file.path(init$data_folder, "datasets", input$name, "consclust", paste0(input$name,'_affectation_k',input$nclust, ".RData")))
  # load(file=file.path(init$data_folder, "datasets", input$name, "supervised", paste0(input$name, "_", input$nclust, "_", input$qval.th, "_", input$cdiff.th, "_", input$de_type, ".RData")))
  # #For after diff analysis
  # diff = list()
  # my.res_save -> diff$my.res
  # summary_save -> diff$summary
  # groups_save -> diff$groups
  # refs_save -> diff$refs


## 0.6 Data folder initialization ##

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

    #If matrix in new format, comment two lines below :
    rownames(datamatrix_single) = as.character(datamatrix_single[,1])
    datamatrix_single = datamatrix_single [,-1]

    # & uncomment those
    # rownames(datamatrix_single) <- gsub("-", "_", rownames(datamatrix_single))
    # rownames(datamatrix_single) <- gsub(":", "_", rownames(datamatrix_single))


    total_cell <- length(datamatrix_single[1,])
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

  #Removing weird chromosomes
  splitID <- sapply(rownames(datamatrix), function(x) strsplit(as.character(x), split="_"))
  normalChr <- which(sapply(splitID, length) <= 3) # weird chromosomes contain underscores in the name
  datamatrix <- datamatrix[normalChr,]
  
  #Remove chrM from mat if it is inside
  if(length(grep("chrM",rownames(datamatrix)))>0)  datamatrix <- datamatrix[-grep("chrM",rownames(datamatrix)),]

  save(datamatrix, annot_raw, file=file.path(init$data_folder, "datasets", input$name, "scChIP_raw.RData"))


##############################################################################################################################
### 1. Filtering and QC
##############################################################################################################################

print("Running filtering and QC...")

## 1.1 Data loading ##

  assays = list(counts = as.matrix(datamatrix)) #, colData =  annot_raw

  umi <- SingleCellExperiment(assays = assays, colData = annot_raw)
  umi <- umi[rowSums(counts(umi)>0)>0, ] # remove windows that do not have any read in any cells
  umi <- scater::calculateQCMetrics(umi)
  thresh <- quantile(colSums(counts(umi)), probs=seq(0, 1, 0.01))

## 1.2 Filtering & Window selection  ##


  #Cell selection based on number of total counts, between min_cov_cell and upper 5%
  sel1000 =  (colSums(counts(umi))>1000 & colSums(counts(umi))< thresh[input$quant_removal+1])

  sel <- ( colSums(counts(umi))>input$min_coverage_cell & colSums(counts(umi)) < thresh[input$quant_removal+1] )

  annot <- colData(umi)


  SelMatCov1000 <- counts(umi)[,sel1000]
  bina_counts <- SelMatCov1000
  bina_counts[bina_counts<2] <-0
  bina_counts[bina_counts>1] <-1
  fixedWin <- names(which((rowSums(bina_counts) > ( (input$min_cells_window/100.0)*(dim(bina_counts)[2])) ))) # window selection
  length(fixedWin)

  SelMatCov <- counts(umi)[,sel]
  SelMatCov <- SelMatCov[fixedWin,]

  annot <- colData(umi)
  annot <- as.data.frame(annot[sel,])
  annot = cbind(annot,annot_raw[which(annot_raw$cell_id %in% rownames(annot)),])

  #Removing user specified regions
  if(!is.null(input$exclude)){
    exclude_regions <- setNames(read.table(input$exclude, header=FALSE, stringsAsFactors=FALSE), c("chr", "start", "stop"))
    dim(exclude_regions)
    regions <- data.frame(loc=rownames(SelMatCov))
    regions <- separate(regions, loc, into=c("chr", "start", "stop"), sep="_", convert=TRUE)
    reg_gr <- makeGRangesFromDataFrame(regions, ignore.strand=TRUE, seqnames.field=c("chr"), start.field=c("start"), end.field="stop")
    excl_gr <- makeGRangesFromDataFrame(exclude_regions, ignore.strand=TRUE, seqnames.field=c("chr"), start.field=c("start"), end.field="stop")
    ovrlps <- as.data.frame(findOverlaps(reg_gr, excl_gr))[, 1]
    SelMatCov <- SelMatCov[-unique(ovrlps), ]
  }
  

  mat <- NULL
  mat <- mean(colSums(SelMatCov))*t(t(SelMatCov)/colSums(SelMatCov))

  norm_mat = mat

  save(norm_mat, file=file.path(init$data_folder, "datasets", input$name, "reduced_data", paste(input$name, input$min_coverage_cell,input$min_cells_window, input$quant_removal, "uncorrected", sep="_", "normMat.RData"))) # used for supervised analysis

  mat <- mat-apply(mat, 1, mean)

## 1.3 Feature annotation ##

  print(paste0("Making annotation"))

  ID <- rownames(norm_mat)
  chr <- unlist(lapply(strsplit(ID, split= "_"), FUN=function(x) unlist(x)[1]))
  start <- unlist(lapply(strsplit(ID, split= "_"), FUN=function(x) unlist(x)[2]))
  end <- unlist(lapply(strsplit(ID, split= "_"), FUN=function(x) unlist(x)[3]))
  feature <- data.frame(ID=rownames(norm_mat), chr=chr, start=start, end=end)
  write.table(feature[, 2:4], file=("feature.bed"), sep="\t", row.names=F, col.names=F, quote=F)
  system("bedtools sort -i feature.bed > featuresort.bed")
  system(paste0("bedtools closest -a featuresort.bed -b ", file.path("~/Documents/GitLab/single_cell_ChIPseq/Shiny_scChIP/annotation", input$annotation_id, "Gencode_TSS_pc_lincRNA_antisense.bed"), " -wa -wb -d> out.bed"))
  annotFeat <- read.table("out.bed", header=F)
  unlink(c("feature.bed", "featuresort.bed", "out.bed"))
  annotFeat <- annotFeat[, c(1:3, 7:8)]
  colnames(annotFeat) <- c("chr", "start", "end", "Gene", "distance")
  annotFeat$ID <- paste(annotFeat$chr, annotFeat$start, annotFeat$end, sep="_")
  annotFeat <- annotFeat %>% group_by(ID) %>% summarise_all(funs(paste(unique(.), collapse = ', '))) %>% as.data.frame
  save(annotFeat, file=file.path(init$data_folder, "datasets", input$name, "reduced_data", paste0(paste(input$name, input$min_coverage_cell, input$min_cells_window, input$quant_removal, "uncorrected", sep="_"), "_annotFeat.RData"))) # used for supervised analysis

  tmp_meta <- data.frame(Sample=rownames(annot), sample_id=annot$sample_id) # modify if coloring should be possible for other columns
  anocol <- geco.annotToCol2(annotS=annot[, annotCol], annotT=annot, plotLegend=T, plotLegendFile=file.path(init$data_folder, "datasets", input$name, "Annotation_legends.pdf"), categCol=NULL)
  annotColors <- data.frame(sample_id=as.data.frame(anocol)$sample_id) %>% setNames(str_c(names(.), "_Color_Orig")) %>% rownames_to_column("Sample") %>% left_join(tmp_meta,. , by="Sample")

## 1.5 PCA ##

  print("Running pca ...")
  pca <- stats::prcomp(t(mat),center=F,scale.=F)
  pca = pca$x[,1:50]


  extendXlimLeft <- 5 # Extend left Xlim by x percent of min
  extendXlimRight <- 15 # Extend right Xlim by x percent of max
  extendYlimBottom <- 5 # Extend bottom Ylim by x percent of min
  extendYlimTop <- 5 # Extend top Ylim by x percent of max
  TextSize <- 0.4
  pcaText <- FALSE

  pdf(file.path(init$data_folder, "datasets", input$name, "reduced_data", paste0(paste(input$name, input$min_coverage_cell, input$min_cells_window, input$quant_removal, "PCA", sep="_"), ".pdf")), height=5, width=5)

  #PCA colored by sample_id
  plot(pca[,1],pca[,2],col=alpha(anocol[,'sample_id'],0.6),xlab='PC1',ylab='PC2',cex=0.8,lwd=1.5, main=paste0('PCA colored by ',colnames(anocol)[1]), pch=19,
   )

  #PCA colored by counts
  plot(pca[,1],pca[,2],col=alpha(anocol[,'total_counts'],0.6),xlab='PC1',ylab='PC2',cex=0.8,lwd=1.5, main=paste0('PCA colored by ',colnames(anocol)[2]), pch=19,
  )
  dev.off()

  print("PCA done !")

## 1.6 t-SNE ##

  print("Running t-sne ...")

  #Reduce the perplexity if the number of samples is too low to avoid perplexity error
  tsne <- Rtsne(pca[,1:50], dims=2, pca=FALSE, theta=0.0, perplexity=choose_perplexity(pca), verbose=FALSE, max_iter=1000)

  pdf(file.path(init$data_folder, "datasets", input$name, "reduced_data", paste0(paste(input$name, input$min_coverage_cell, input$min_cells_window, input$quant_removal, "Tsne", sep="_"), ".pdf")), height=5, width=5)

  #T-sne colored by sample_id
  p <- ggplot(as.data.frame(tsne$Y), aes(x=V1, y=V2)) + geom_point(alpha=0.6, aes(color=annot[,'sample_id'])) +
    labs(color='sample_id', x="t-SNE 1", y="t-SNE 2") +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          panel.background=element_blank(), axis.line=element_line(colour="black"),
          panel.border=element_rect(colour="black", fill=NA)) + ggtitle("T-SNE colored by sample id")
  p <- p + scale_color_manual(values = levels(as.factor(unique(anocol[,'sample_id'])))) + theme(legend.position = "none")
  p

  #T-sne colored by counts
  p <- ggplot(as.data.frame(tsne$Y), aes(x=V1, y=V2)) + geom_point(alpha=0.6, aes(color=annot[,'total_counts'])) +
    labs(color='sample_id', x="t-SNE 1", y="t-SNE 2") +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          panel.background=element_blank(), axis.line=element_line(colour="black"),
          panel.border=element_rect(colour="black", fill=NA)) +  ggtitle("T-SNE colored by total count")

  p <- p + scale_color_gradientn(colours = matlab.like(100))
  p

  dev.off()

  print("T-sne done !")

## 1.7 Save data ##

  save(pca,mat, annot, tsne, file=file.path(init$data_folder, "datasets", input$name, "reduced_data", paste0(paste(input$name, input$min_coverage_cell, input$min_cells_window, input$quant_removal, "uncorrected", sep="_"), ".RData")))

##############################################################################################################################
### 2. Correlation clustering & Filtering
##############################################################################################################################

print("Running correlation clustering & filtering... ")

## 2.1 Determining correlation threshold ##

  mati <- as.matrix(t(pca[,1:50]))
  hc_cor <- hclust(as.dist(1 - cor(mati)), method="ward.D")
  mat.so.cor <- mati[,hc_cor$order]


  correlation_values <- list(limitC=vector(length=100))
  corChIP <-  cor(mati)

  for(i in 1:500){
      random_mat <-  matrix(sample(mati), nrow=dim(mati)[1])
      thresh2 <- quantile(cor(random_mat), probs=seq(0,1,0.01))
      limitC <-  thresh2[input$corr_thresh+1]
      correlation_values$limitC[i] = limitC
  }

  correlation_values$limitC_mean = mean(correlation_values$limitC,na.rm = T)

## 2.2 Filtering cells based on correlation threshold ##


  cf$sel2 <- (apply(corChIP, 1, function(x) length(which(x>correlation_values$limitC_mean))) > (input$percent_corr*0.01)*dim(corChIP)[1])
  cf$mati2 <- mati[, cf$sel2]
  print(dim(cf$mati2))



  tmp_meta <- data.frame(Sample=rownames(annot), sample_id=annot$sample_id) # modify if coloring should be possible for other columns
  anocol <- geco.annotToCol2(annotS=annot[, annotCol], annotT=annot, plotLegend=T, plotLegendFile=file.path(init$data_folder, "datasets", input$name, "Annotation_legends.pdf"), categCol=NULL)
  annotColors <- data.frame(sample_id=as.data.frame(anocol)$sample_id) %>% setNames(str_c(names(.), "_Color_Orig")) %>% rownames_to_column("Sample") %>% left_join(tmp_meta,. , by="Sample")

  cf$hc_cor2 <- hclust(as.dist(1 - cor(cf$mati2)), method="ward.D")
  cf$mat.so.cor2 <- cf$mati2[, cf$hc_cor2$order]
  cf$annot_sel <- annot[cf$sel2,]


  cf$anocol_sel <- geco.annotToCol2(annotS=cf$annot_sel[, annotCol], annotT=cf$annot_sel, plotLegend=T, plotLegendFile=file.path(init$data_folder, "datasets", input$name,"Annotation_legends_reclustering.pdf"), categCol=NULL)
  anocol_sel <- cf$anocol_sel

  mati2 <- cf$mati2
  print(dim(mati2))

  annot_sel <- cf$annot_sel
  mat.so.cor2 <- cf$mat.so.cor2
  hc_cor2 <- cf$hc_cor2
  save(mati2,mati,correlation_values,corChIP, annot_sel, mat.so.cor2, hc_cor2, file=file.path(init$data_folder, "datasets", input$name, "cor_filtered_data", paste0(input$name, "_", input$corr_thresh, "_", input$percent_corr, ".RData")))

## 2.3 Plotting correlation distribution ##

  pdf(file.path(init$data_folder, "datasets", input$name, "cor_filtered_data","CorrDistrib.pdf"),height=4.5,width=4.5)
  hist(corChIP,prob=TRUE,col=alpha("red",0.8),breaks=50,ylim=c(0,4),main="Distribution of cell to cell correlation scores",xlab="Pearson Corr Scores")
  lines(density(corChIP),col="red",lwd=2)
  lines(density(cor(random_mat)),col=1,lwd=2)
  abline(v=correlation_values$limitC,lwd=2)
  dev.off()

## 2.4 Plotting correlation heatmap before filtering ##

  hmColors <- colorRampPalette(c("royalblue","white","indianred1"))(256)
  corColors <- colorRampPalette(c("royalblue","white","indianred1"))(256)
  png(file.path(init$data_folder, "datasets", input$name, "cor_filtered_data","Clustering_correlation_matrix_PCA_wo_lowcorr.png"), height=1000,width=1000,res=150)
  geco.hclustAnnotHeatmapPlot(x=cor(mat.so.cor),
                              hc=hc_cor,
                              hmColors=corColors,
                              anocol= anocol[hc_cor$order,],
                              xpos=c(0.15,0.9,0.164,0.885),
                              ypos=c(0.1,0.5,0.5,0.6,0.62,0.95),
                              dendro.cex=0.05,
                              xlab.cex=0.5,
                              hmRowNames=FALSE,
                              hmRowNames.cex=0.01
  )
  dev.off()

## 2.4 Plotting correlation heatmap after filtering ##

  png(file.path(init$data_folder, "datasets", input$name, "cor_filtered_data","Clustering_correlation_matrix_PCA_wo_highcorr.png"), height=1000,width=1000,res=150)

  geco.hclustAnnotHeatmapPlot(x=cor(cf$mat.so.cor2),
                                   hc=cf$hc_cor2,
                                   hmColors=corColors,
                                   anocol=cf$anocol_sel[cf$hc_cor2$order,],
                                   xpos=c(0.15, 0.9, 0.164, 0.885),
                                   ypos=c(0.1, 0.5, 0.5, 0.6, 0.62, 0.95),
                                   dendro.cex=0.04,
                                   xlab.cex=0.8,
                                   hmRowNames=FALSE
    )
  dev.off()

print("Correlation filtering done...")

##############################################################################################################################
### 3. Consensus clustering
##############################################################################################################################

print("Running consensus clustering ...")

## 3.1 Running consensus clustering ##

  if(is.null(input$nclust)) input$maxK = 10 else input$maxK = input$nclust+1

  consclust <- ConsensusClusterPlus(mati2, maxK=input$maxK, reps=1000, pItem=0.8, pFeature=1,
                                    title="Consensus_clustering_dir", clusterAlg="hc", distance="pearson",
                                    innerLinkage="ward.D", finalLinkage="ward.D", seed=3.14, verbose = F,plot="pdf")

  icl <- calcICL(consclust, plot="png", title="Consensus_clustering_dir")

## 3.2 Saving plot and data ##

  file.copy(from = "Consensus_clustering_dir", to = file.path(init$data_folder, "datasets", input$name,"consclust"), recursive=TRUE)
  unlink("Consensus_clustering_dir", recursive=TRUE)

  save(consclust, icl, file=file.path(init$data_folder, "datasets", input$name, "consclust", paste0(input$name, ".RData")))


## 3.3 Read in number of cluster to use for differential analysis ##

  if (is.null(input$nclust)){

  input$nclust = -1
  while(is.na(input$nclust) | as.integer(input$nclust)<2 | as.integer(input$nclust) > 10){
    cat (paste0("How much clusters do you want to set to perform differential analysis ? \n(Choose based on the plots created in 'datasets/",input$name,"/consclust/' )\n"))
    input$nclust = as.integer(readline(prompt = " k = "))
    if (is.na(input$nclust) | input$nclust<2 | input$nclust > 10) {
      cat("Please choose a number of clusters between 2 and 10 \n")
    }
  }
  }

## 3.4 Affecting cells to each class based on consensus clustering ##

  anocol_sel <- geco.annotToCol2(annotS=annot_sel[, annotCol], annotT=annot_sel, plotLegend=T, plotLegendFile=file.path(init$data_folder, "datasets", input$name,"Annotation_legends_reclustering.pdf"), categCol=NULL)

  clustCol <- "ChromatinGroup"
  so <- colnames(mat.so.cor2)
  conscol <- c("#1F78B4", "#33A02C","#A6CEE3", "#B2DF8A", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#bd18ea", "#2ef4ca", "#f4cced", "#f4cc03", "#05188a", "#e5a25a", "#06f106", "#85848f", "#000000", "#076f25", "#93cd7f", "#4d0776", "#ffffff")
  cc <- consclust[[as.integer(input$nclust)]]$consensusClass
  #cc = cc[sel_itemCons]
  cc. <- lapply(unique(cc), function(z) names(which(cc==z)))
  mat.cc <- geco.groupMat(mati2, margin=1, groups=cc., method="mean")
  hcc <- hclust(distPearson(t(mat.cc)), method="ward.D")
  clust$annot_sel2 <- annot_sel
  clust$annot_sel2[, clustCol] <- paste("C", match(cc, hcc$order), sep="")
  clust$cc.col <- cbind(ChromatinGroup=conscol[match(clust$annot_sel2[so, clustCol], paste("C", 1:as.integer(input$nclust), sep=""))], anocol_sel[so,])

  affectation <- clust$annot_sel2[,c("barcode", "cell_id", "ChromatinGroup", "sample_id")]
  save(affectation, file=file.path(init$data_folder, "datasets", input$name, "consclust", paste0(input$name, '_affectation_k', input$nclust, ".RData")))


print("Consensus clustering done...")


##############################################################################################################################
### 4. Differential Analysis
##############################################################################################################################

print("Running differential analysis...")

## 4.1 Calculating wilcoxon rank parameter for each cluster ##

  diff <- list(my.res=NULL, summary=NULL, groups=NULL, refs=NULL)
  Counts <- norm_mat[, rownames(affectation)]
  feature <- data.frame(ID=annotFeat$ID, chr=annotFeat$chr, start=annotFeat$start, end=annotFeat$end)
  # compare each cluster to all the rest
  mygps <- lapply(1:input$nclust, function(i){ affectation[which(affectation$ChromatinGroup==paste0("C", i)), "cell_id"]})
  names(mygps) <- paste0('C', 1:input$nclust)
  groups <- names(mygps)
  myrefs <- lapply(1:input$nclust, function(i){ affectation[which(affectation$ChromatinGroup!=paste0("C", i)), "cell_id"]})
  names(myrefs) <- paste0('notC', 1:input$nclust)
  refs <- names(myrefs)
  diff$my.res <- geco.CompareWilcox(dataMat=Counts, annot=affectation, ref=myrefs, groups=mygps, featureTab=feature)

## 4.2 Retrieving significatively enriched and depleted regions ##

  diff$summary <- matrix(nrow=3, ncol=length(groups), dimnames=list(c("differential", "over", "under"), groups))
  for(k in 1:length(groups)){
      gpsamp <- groups[k]
      diff$summary["differential", gpsamp] <- sum(diff$my.res[, paste("qval", gpsamp, sep=".")] <= input$qval.th & abs(diff$my.res[, paste("cdiff", gpsamp, sep=".")]) > input$cdiff.th, na.rm=T)
      diff$summary["over", gpsamp] <- sum(diff$my.res[, paste("qval", gpsamp, sep=".")] <= input$qval.th & diff$my.res[, paste("cdiff", gpsamp, sep=".")] > input$cdiff.th, na.rm=T)
      diff$summary["under", gpsamp] <- sum(diff$my.res[, paste("qval", gpsamp, sep=".")] <= input$qval.th & diff$my.res[, paste("cdiff", gpsamp, sep=".")] < -input$cdiff.th, na.rm=T)

  }
  diff$groups <- groups
  diff$refs <- refs
  my.res_save <- diff$my.res
  summary_save <- diff$summary
  groups_save <- diff$groups
  refs_save <- diff$refs
  save(my.res_save, summary_save, groups_save, refs_save, file=file.path(init$data_folder, "datasets", input$name, "supervised", paste0(input$name, "_", input$nclust, "_", input$qval.th, "_", input$cdiff.th, "_", input$de_type, ".RData")))

## 4.3 Plotting differential analysis barplot & volcano plots ##

  #Barplot
  png(file.path(init$data_folder,"datasets" , input$name, 'supervised','diff_barplot.png'))
  myylim <- range(c(diff$summary["over",], -diff$summary["under", ]))
  barplot(diff$summary["over",], col="red", las=1, ylim=myylim, main="Differentially bound regions",
          ylab="Number of regions", axes=F)
  barplot(-diff$summary["under", ], col="forestgreen", ylim=myylim, add=T, axes=F, names.arg="")
  z <- axis(2, pos=-10)
  axis(2, at=z, labels=abs(z), las=1)
  dev.off()

  #Volcano plots
  for(k in 1:length(groups)){

    gpsamp <- groups[k]
    png(file.path(init$data_folder,"datasets" , input$name, 'supervised',paste0('volcano_plot_',gpsamp)))
      mycol <- rep("black", nrow(diff$my.res))
      mycol[which(diff$my.res[, paste("qval", gpsamp, sep=".")] <= input$qval.th & diff$my.res[, paste("cdiff", gpsamp, sep=".")] > input$cdiff.th)] <- "red"
      mycol[which(diff$my.res[, paste("qval", gpsamp, sep=".")] <= input$qval.th & diff$my.res[, paste("cdiff", gpsamp, sep=".")] < -input$cdiff.th)] <- "forestgreen"
      plot(diff$my.res[, paste("cdiff", gpsamp, sep=".")], -log10(diff$my.res[, paste("qval", gpsamp, sep=".")]),
           col=mycol, cex=0.7, pch=16, xlab="count difference", ylab="-log10(adjusted p-value)", las=1,
           main=paste(input$gpsamp, "vs the rest","\n", diff$summary["over",input$gpsamp], "enriched,", diff$summary["under", input$gpsamp], "depleted"))
      abline(v=input$cdiff.th, lty=2)
      abline(h=-log10(input$qval.th), lty=2)
      abline(v=-input$cdiff.th, lty=2)
    dev.off()
  }

print("Differential analysis done...")

#If user inputed bam files -> continue towards peak calling and gene set enrichment :
if(file.exists(input$bam1) & file.exists(input$bam2)){

  print("Bam files found, continuing analysis...")

  ##############################################################################################################################
  ### 5. Peak calling
  ##############################################################################################################################

  print("Running peak calling...")


  ## 5.1 Creating peak calling folder ##

  dir.create(file.path(init$data_folder, "datasets",  input$name, "peaks"), showWarnings=FALSE)
  dir.create(file.path(init$data_folder, "datasets",  input$name, "peaks", paste0( input$name, "_k", input$nclust)), showWarnings=FALSE)
  odir = file.path(init$data_folder, "datasets", input$name, "peaks",paste0(input$name,"_k",input$nclust))
  sample_ids <- unique(annot_sel$sample_id)
  input$pc_stat="p.value"
  stat.value <- if(input$pc_stat=="p.value") paste("-p", input$pc_stat_value) else paste("-q", input$pc_stat_value)
  inputBams = as.vector(unlist(inputBams))

  ## 5.2 Merging bam files together ##

  if(length(inputBams) > 1) {
    write(inputBams, file=file.path(odir, "bam_list.txt"))
    system(paste0('samtools merge -@ 4 -f -h ', inputBams[1],' -b ', file.path(odir, "bam_list.txt"), ' ', file.path(odir, 'merged.bam')))
    merged_bam=file.path(odir, 'merged.bam')
  }

  ## 5.3 Writing affected clusters ##

  for(class in levels(factor(affectation$ChromatinGroup))){
    write(as.vector(affectation$barcode[which(affectation$ChromatinGroup == as.character(class))]), file=file.path(odir, paste0(class, ".barcode_class")))
  }
  write(levels(factor(affectation$ChromatinGroup)),file = file.path(odir, "barcodes.barcode_class"))

  ## 5.4 Splitting bam files into clusters ##

  system(paste0('samtools view -H ', merged_bam, ' > ', file.path(odir, 'header.sam')))
  system(paste0('for i in $(cat ', file.path(odir, 'barcodes.barcode_class'), '); do samtools view -h ',merged_bam,' | fgrep -w -f ', file.path(odir,'/$i.barcode_class'), ' > ', file.path(odir,'$i.sam'), ';done'))

  #Reconvert to bam
  system(paste0('for i in $(cat ', file.path(odir,'barcodes.barcode_class'), '); do cat ', file.path(odir, 'header.sam'), ' ', file.path(odir, '$i.sam'), ' | samtools view -b - > ', file.path(odir, '$i.bam'), ' ; done'))

  #BamCoverage
  #system(paste0('for i in $(cat ', file.path(odir,'barcodes.barcode_class'), '); do samtools index ', file.path(odir,'$i.bam'), '; done'))
  #system(paste0('for i in $(cat ', file.path(odir,'barcodes.barcode_class'), '); do bamCoverage --bam ', file.path(odir,'$i.bam'), ' --outFileName ', file.path(odir,'$i.bw'), ' --binSize 50 --smoothLength 500 --extendReads 150 --ignoreForNormalization chrX --numberOfProcessors 4 --normalizeUsing RPKM; done'))

  system(paste0('rm ', file.path(odir,'*.barcode_class'), ' ', file.path(odir,'*.sam')))

  ## 5.5 Call peak on each cluster ##

  #Peak calling with macs2
  stat.value = '-p 0.05'
  for(cluster in levels(factor(affectation$ChromatinGroup))){
    system(paste0('macs2 callpeak ', stat.value, ' --broad -t ', file.path(odir, paste0(cluster,".bam")), " --outdir ", odir," --name ",cluster))
    system(paste0('bedtools merge -delim "\t" -d 20000 -i ', file.path(odir, paste0(cluster, '_peaks.broadPeak')), ' > ', file.path(odir, paste0(cluster, '_merged.bed'))))
  }


  #Clean up files
  unlink(file.path(odir, "bam_list.txt"))
  unlink(file.path(odir, "*.bam"))
  unlink(file.path(odir, "*.bam.bai"))
  unlink(file.path(odir, "*.xls"))
  unlink(file.path(odir, "*.gappedPeak"))
  unlink(file.path(odir, "*_model.r"))

  ## 5.6 Make annotation files ##

  mergeBams <- paste(sapply(levels(factor(affectation$ChromatinGroup)), function(x){ file.path(odir, paste0(x, "_merged.bed")) }), collapse = ';')
  mergeBams <- paste0('"', mergeBams, '"')
  system(paste("bash", file.path("Modules", "makePeakAnnot.sh"), mergeBams, file.path("annotation", input$annotation_id, "chrom.sizes.bed"), file.path("annotation", input$annotation_id, "50k.bed"), file.path("annotation", input$annotation_id, "Gencode_TSS_pc_lincRNA_antisense.bed"), paste0('"', odir, .Platform$file.sep, '"')))

  print("Peak calling done...")

  ##############################################################################################################################
  ### 6. Enrichment Analysis
  ##############################################################################################################################
  
  print("Running enrichment analysis...")
  
  ## 6.1 Loading annotation files ##
  
    myData = new.env()
    load(file.path(init$data_folder,"annotation" , input$annotation_id, 'MSigDB.RData'), envir=myData)
    MSIG.ls <-  myData$MSIG.ls
    MSIG.gs <-  myData$MSIG.gs
    annotFeat_long <-  as.data.frame(cSplit(annotFeat, splitCols="Gene", sep=", ", direction="long"))
  
    pm.annot.window.file <- file.path(init$data_folder, 'datasets', input$name, 'peaks', paste0(input$name, '_k', input$nclust), 'pm.annot.window.bed')
    peak_window <- read.table(pm.annot.window.file, sep="\t", header=F)
    colnames(peak_window) <- c("chr", "start", "end", "w_chr", "w_start", "w_end")
    peak_window$peak_ID <- paste(peak_window$chr, peak_window$start, peak_window$end, sep="_")
    peak_window$window_ID <- paste(peak_window$w_chr, peak_window$w_start, peak_window$w_end, sep="_")
  
    pm.annot.gene.file <- file.path(init$data_folder, 'datasets', input$name, 'peaks', paste0(input$name, '_k', input$nclust), 'pm.annot.gene.bed')
    peak_gene <- read.table(pm.annot.gene.file, sep="\t", header=F)
    peak_gene <- peak_gene[, c(1:3, 7:8)]
    colnames(peak_gene) <- c("chr","start", "end", "Gene_name", "dist_TSS")
    peak_gene$peak_ID <- paste(peak_gene$chr, peak_gene$start, peak_gene$end, sep="_")
  
    Gencode <- read.table(file.path('annotation', input$annotation_id, 'Gencode_TSS_pc_lincRNA_antisense.bed'))
    GencodeGenes <- unique(Gencode$V4)
  
  ## 6.2 Running enrichment test for each cluster ##
  
    Both <- list()
    Enriched <- list()
    Depleted <- list()
  
    for(i in 1:length(diff$groups)){
          gp <- diff$groups[i]
          ref <- diff$refs[i]
          signific <- diff$my.res$ID[which(diff$my.res[, paste("qval", gp, sep=".")] <= input$qval.th & abs(diff$my.res[, paste("cdiff", gp, sep=".")]) > input$cdiff.th)]
          significG <- unique(annotFeat_long$Gene[annotFeat_long$distance < 1000 & annotFeat_long$ID %in% signific])
          over <- diff$my.res$ID[which(diff$my.res[, paste("qval", gp, sep=".")] <= input$qval.th & diff$my.res[, paste("cdiff", gp, sep=".")] > input$cdiff.th)]
          overG <- unique(annotFeat_long$Gene[annotFeat_long$distance < 1000 & annotFeat_long$ID %in% over])
          under <- diff$my.res$ID[which(diff$my.res[, paste("qval", gp, sep=".")] <= input$qval.th & diff$my.res[, paste("cdiff", gp, sep=".")] < -input$cdiff.th)]
          underG <- unique(annotFeat_long$Gene[annotFeat_long$distance < 1000 & annotFeat_long$ID %in% under])
  
          signific_associated_peak <- peak_window[peak_window$window_ID %in% signific, 8]
          over_associated_peak <- peak_window[peak_window$window_ID %in% over, 8]
          under_associated_peak <- peak_window[peak_window$window_ID %in% under, 8]
          signific_associated_gene <- peak_gene[peak_gene$peak_ID %in% signific_associated_peak & peak_gene$dist_TSS<1000, ]
          over_associated_gene <- peak_gene[peak_gene$peak_ID %in% over_associated_peak & peak_gene$dist_TSS<1000, ]
          under_associated_gene <- peak_gene[peak_gene$peak_ID %in% under_associated_peak & peak_gene$dist_TSS<1000, ]
          significG <- unique(signific_associated_gene$Gene_name)
          overG <- unique(over_associated_gene$Gene_name)
          underG <- unique(under_associated_gene$Gene_name)
          
          # if(length(significG)){
          #   enrich.test <- geco.enrichmentTest(gene.sets=MSIG.ls, mylist=significG, possibleIds=GencodeGenes)
          #   enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
          #   enrich.test <- merge(subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
          #   enrich.test <- enrich.test[order(enrich.test$`p-value`),]
          #   ind <- which(enrich.test$`q-value`<= 0.1)
          #   if(!length(ind)){ind <- 1:10}
          #   Both[[i]] <- enrich.test[ind,]
          # }
          if(length(overG)){
            enrich.test <- geco.enrichmentTest(gene.sets=MSIG.ls, mylist=overG, possibleIds=GencodeGenes)
            enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
            enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
            enrich.test <- enrich.test[order(enrich.test$`p-value`),]
            ind <- which(enrich.test$`q-value`<= 0.1)
            if(!length(ind)){ind <- 1:10}
            Enriched[[i]]  <- enrich.test[ind,]
          }
          
          if(length(underG)){
            enrich.test <- geco.enrichmentTest(gene.sets=MSIG.ls, mylist=underG, possibleIds=GencodeGenes)
            enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
            enrich.test <- merge( subset(MSIG.gs, select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
            enrich.test <- enrich.test[order(enrich.test$`p-value`),]
            ind <- which(enrich.test$`q-value`<= 0.1)
            if(!length(ind)){ind <- 1:10}
            Depleted[[i]] <- enrich.test[ind,]
          }
          
    }
    
    enr <- list(Both=NULL, Enriched=NULL, Depleted=NULL)
    # enr$Both <- Both
    enr$Enriched <- Enriched
    enr$Depleted <- Depleted
    
  
  ##  6.3 Saving enrichment tables ##
  
    for(i in 1:length(diff$groups)){
          # if(!is.null(enr$Both[[i]])){
          #   filename <- file.path(init$data_folder,"datasets",input$name,paste0(diff$groups[i], "_significant_gene_sets.csv"))
          #   write.table(enr$Both[[i]], file=filename, quote=FALSE, row.names=FALSE, sep=",")
          # }
          if(!is.null(enr$Enriched[[i]])){
            filename <- file.path(init$data_folder,"datasets",input$name,"supervised",paste0(diff$groups[i], "_enriched_gene_sets.csv"))
            write.table(enr$Enriched[[i]], file=filename, quote=FALSE, row.names=FALSE, sep=",")
          }
          if(!is.null(enr$Depleted[[i]])){
            filename <- file.path(init$data_folder,"datasets",input$name,"supervised",paste0(diff$groups[i], "_depleted_gene_sets.csv"))
            write.table(enr$Depleted[[i]], file=filename, quote=FALSE, row.names=FALSE, sep=",")
           }
    }
    
  ## 6.4 Plotting most depleted/enriched regions ##
    classes_MSIG = as.vector(read.table("annotation/MSIGdb_classes")[,1])
    
    pdf(file.path(init$data_folder,"datasets" , input$name, 'supervised',paste0('Gene_Enrichment_Analysis_Over.pdf')))
    for(i in 1:length(diff$groups)){
      if(!is.null(enr$Enriched[[i]])) {
        over = as_tibble(enr$Enriched[[i]])
        over = over %>% dplyr::filter(Class %in% classes_MSIG)
        if(dim(over)[1]>0){
          over_plot = ggplot(data =over[1:min(10,length(over$Gene.Set)),], aes(x=sort(as.factor(Gene.Set)),label = over$Gene.Set[1:min(10,length(over$Gene.Set))],y= -log10(`q-value`))) + geom_col(fill="#B8B8B8") +geom_text(aes(y = 0), angle = 90, hjust = -.05, size = 2.5) + theme(axis.text.x = element_blank()) + xlab("GeneSets") + ggtitle(paste0("Top 10 Enriched Gene Sets - C",i))
          print(over_plot)
        }
      }
    }
    dev.off()
    
    pdf(file.path(init$data_folder,"datasets" , input$name, 'supervised',paste0('Gene_Enrichment_Analysis_Under.pdf')))
    for(i in 1:length(diff$groups)){
      if(!is.null(enr$Depleted[[i]])) {
        under = as_tibble(enr$Depleted[[i]])
        under = under %>% dplyr::filter(Class %in% classes_MSIG)
        if(dim(under)[1]>0){
          under_plot = ggplot(data =under[1:min(10,length(under$Gene.Set)),], aes(x=sort(as.factor(Gene.Set)),label = under$Gene.Set[1:min(10,length(under$Gene.Set))],y= -log10(`q-value`))) + geom_col(fill="#B8B8B8") +geom_text(aes(y = 0), angle = 90, hjust = -.05, size = 2.5) + theme(axis.text.x = element_blank()) + xlab("GeneSets") + ggtitle(paste0("Top 10 Depleted Gene Sets - C",i))
          print(under_plot)
        }
      }
    }
    dev.off()
  
  print("Gene set enrichment done...")
}
##############################################################################################################################
### END
##############################################################################################################################

print("Run successfully finished !")
