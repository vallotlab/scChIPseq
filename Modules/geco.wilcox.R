
geco.CompareWilcox <- function(dataMat=NULL, annot=NULL, ref=NULL, groups=NULL, featureTab=NULL){
  res <- featureTab
  res = res[ order(res$ID), ]
  dataMat=  dataMat[ order(row.names(dataMat)), ]
  for(k in 1:length(groups))
  {
    if(length(ref)==1){refsamp <- ref[[1]]}else{refsamp <- ref[[k]]}
    gpsamp <- groups[[k]]
    annot. <- annot[c(refsamp, gpsamp), 1:2]
    annot.$Condition <- c(rep("ref", length(refsamp)), rep("gpsamp", length(gpsamp)))
    mat. <- dataMat[, c(as.character(refsamp), as.character(gpsamp))]
    
    testWilc <- apply(dataMat, 1, function(x) wilcox.test(as.numeric(x[as.character(refsamp)]), as.numeric(x[as.character(gpsamp)])))
    pval.gpsamp <- unlist(lapply(testWilc, function(x) x$p.value))
    qval.gpsamp <- p.adjust(pval.gpsamp, method = "BH")
    Count.gpsamp <- apply(dataMat, 1, function(x) mean(x[as.character(gpsamp)]))
    cdiff.gpsamp <- apply(dataMat, 1, function(x) log(mean(x[as.character(gpsamp)])/mean(x[as.character(refsamp)]), 2))
    # cdiff1.gpsamp <- apply(dataMat, 1, function(x) mean(x[as.character(gpsamp)]) - 2*mean(x[as.character(refsamp)]))
    # cdiff2.gpsamp <- apply(dataMat, 1, function(x) mean(x[as.character(gpsamp)]) - 0.5*mean(x[as.character(refsamp)]))
    # 
    Rank.gpsamp <- rank(qval.gpsamp) # This is different from the rank used in the Wilcox.test !! 
    
    res <- data.frame(res, Rank.gpsamp, Count.gpsamp, cdiff.gpsamp, pval.gpsamp, qval.gpsamp)
    # res <- data.frame(res, Rank.gpsamp, Count.gpsamp, cdiff1.gpsamp,cdiff2.gpsamp, pval.gpsamp, qval.gpsamp)
    
    colnames(res) <- sub("ref", names(ref)[min(c(k, length(ref)))], sub("gpsamp", names(groups)[k], colnames(res)))		  
  }
  res
}
