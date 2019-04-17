geco.annotToCol2 <- 
function (annotS = NULL, annotT = NULL, missing = c("", NA), 
          anotype = NULL, maxnumcateg = 2, categCol = NULL, quantitCol = NULL, 
          plotLegend = T, plotLegendFile = NULL) 
{
  if (is.null(ncol(annotS))) {
    annotS <- data.frame(annotS)
    colnames(annotS) = annotCol
    rownames(annotS) = rownames(annotT)
  }
  for (j in 1:ncol(annotS)) annotS[which(annotS[, j] %in% missing), 
                                   j] <- NA
  if (is.null(anotype)) {
    anotype <- rep("categ", ncol(annotS))
    names(anotype) <- colnames(annotS)
    classes <- sapply(1:ncol(annotS), function(j) class(annotS[, 
                                                               j]))
    nmodal <- sapply(1:ncol(annotS), function(j) length(unique(setdiff(annotS[, 
                                                                              j], NA))))
    anotype[which(classes %in% c("integer", "numeric") & 
                    nmodal > maxnumcateg)] <- "quantit"
    anotype[which(nmodal == 2)] <- "binary"
  }
  anocol <- annotS
  if (plotLegend) 
    pdf(plotLegendFile)
  if (is.null(categCol)) 
    categCol <- c("royalblue", "palevioletred1", "red", "palegreen4", 
                  "skyblue", "sienna2", "slateblue3", "pink2", "slategray", 
                  "black", "orange", "turquoise4", "yellow3", "orangered4", 
                  "orchid", "palegreen2", "orchid4", "red4", "peru", 
                  "orangered", "palevioletred4", "purple", "sienna4", 
                  "turquoise1")
  k <- 1
  for (j in which(anotype == "categ")) {
    tmp <- as.factor(anocol[, j])
    classes <- as.character(levels(tmp))
    ncat <- length(levels(tmp))
    if (k + ncat > length(categCol)) 
      categCol <- c(categCol, categCol)
    levels(tmp) <- categCol[k:(k + ncat - 1)]
    fill <- as.character(levels(tmp))
    anocol[, j] <- as.character(tmp)
    k <- k + ncat
    if (plotLegend) {
      par(mar = c(0, 0, 0, 0))
      plot(-10, axes = F, xlim = c(0, 5), ylim = c(0, 5), 
           xlab = "", ylab = "")
      legend(1, 5, legend = classes, fill = fill, title = colnames(anocol)[j], 
             xjust = 0.5, yjust = 1)
    }
  }
  memcol <- c()
  for (j in which(anotype == "binary")) {
    new <- setdiff(anocol[, j], c(NA, memcol))
    if (length(new) == 2) {
      memcol <- c(memcol, c("aquamarine", "plum1"))
      names(memcol)[(length(memcol) - 1):length(memcol)] <- sort(new)
    }
    if (length(new) == 1) {
      memcol <- c(memcol, setdiff(c("dodgerblue4", "firebrick"), 
                                  memcol[setdiff(anocol[, j], c(NA, new))]))
      names(memcol)[length(memcol)] <- new
    }
    anocol[, j] <- as.character(anocol[, j])
    for (z in 1:length(memcol)) {
      anocol[which(anocol[, j] == names(memcol)[z]), j] <- memcol[z]
    }
    if (plotLegend) {
      par(mar = c(0, 0, 0, 0))
      plot(-10, axes = F, xlim = c(0, 5), ylim = c(0, 5), 
           xlab = "", ylab = "")
      classes <- intersect(names(memcol), annotS[, j])
      fill <- memcol[classes]
      legend(1, 5, legend = classes, fill = fill, title = colnames(anocol)[j], 
             xjust = 0.5, yjust = 1)
    }
  }
  if (is.null(quantitCol)) 
    quantitCol <- c("darkgreen", "darkblue", 
                    "darkgoldenrod4", "darkorchid4", "darkolivegreen4", 
                    "darkorange4", "darkslategray")
  k <- 1
  for (j in which(anotype == "quantit")) {
    colrange <- matlab.like(100)
    anocol[, j] <- colrange[round(geco.changeRange(anocol[, 
                                                          j], newmin = 1, newmax = 100))]
    if (k < length(quantitCol)) {
      k <- k + 1
    }
    else {
      k <- 1
    }
    if (plotLegend) {
      par(mar = c(8, 2, 5, 1))
      lims <- seq(-1, 1, length.out = 200)
      image(matrix(lims, nc = 1), col = colrange, axes = F, 
            xlab = colnames(anocol)[j])
    }
  }
  if (plotLegend) 
    dev.off()
  for (j in 1:ncol(anocol)) anocol[which(is.na(anocol[, j])), 
                                   j] <- "white"
  as.matrix(anocol)
}