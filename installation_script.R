usePackageBioc <- function(p)
{
  if (!is.element(p, installed.packages()[, 1]))
    BiocManager::install(p, dep = TRUE, version = "3.8")
  require(p, character.only = TRUE)
}
usePackage <- function(p)
{
  if (!is.element(p, installed.packages()[, 1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
usePackageGeco <- function(p)
{
	
  if (!is.element(gsub(".tar.gz", "",p), installed.packages()[, 1]))
    install.packages(file.path("packages",p), repos = NULL, type = "source")
  require(p, character.only = TRUE)
}

#Biocmanager
pkgs_bioc = c("scater",
              "scran",
              "ConsensusClusterPlus",
              "GenomicRanges",
              "IRanges")

for (pkg in pkgs_bioc) {
  usePackageBioc(pkg)
}

#R CRAN
pkgs = c(
  "tibble",
  "dplyr",
  "stringr",
  "irlba",
  "reshape2",
  "Rtsne",
  "DT",
  "tidyr",
  "splitstackshape",
  "DT",
  "tidyr",
  "splitstackshape",
  "rlist",
  "plotly",
  "RColorBrewer",
  "colorRamps",
  "colourpicker",
  "kableExtra",
  "knitr",
  "viridis",
  "ggplot2",
  "gplots",
  "png",
  "gridExtra"
)
for (pkg in pkgs) {
  usePackageBioR(pkg)
}

#geco local packages
pkgs_geco = c(
  "geco.utils.tar.gz",
  "geco.visu.tar.gz",
  "geco.unsupervised.tar.gz",
  "geco.supervised.tar.gz"
)
for (pkg in pkgs_geco) {
  usePackageGeco(pkg)
}



#Monocle
if (!is.element("monocle", installed.packages()[, 1])) {
  source("http://bioconductor.org/biocLite.R")
  
  biocLite()
  biocLite("monocle")
  devtools::install_github("cole-trapnell-lab/DDRTree", ref = "simple-ppt-like")
  devtools::install_github("cole-trapnell-lab/L1-graph")
  install.packages("reticulate")
  library(reticulate)
  py_install('umap-learn',
             pip = T,
             pip_ignore_installed = T) # Ensure the latest version of UMAP is installed
  py_install("louvain")
  
  devtools::install_github("cole-trapnell-lab/monocle-release", ref = "monocle3_alpha")
}
