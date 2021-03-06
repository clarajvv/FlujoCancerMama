##########################################################
########################## DEA ###########################
##########################################################

##########################################################
### Libraries
source("Code/Settings/0_loadLibraries.R")
loadpkg("dplyr")
loadpkg("limma") # Differential gene expression analysis (DEA)
loadpkg("edgeR")
loadpkg("calibrate") # To label the volcano plot

##########################################################
### Loading data for DEA
load("Data/datosParaDEA.RData")

if(is.null(preDEA)){
  quit()
}else{
  design <- preDEA["design"][[1]]
  counts <- preDEA["counts"][[1]]
  label <- preDEA["label"][[1]]  
}



# -------------------------------------------------------
# Le pasamos el design, count y la combinaci�n


dge <- DGEList(counts=counts) #con datos de expresion genica digital
A <- rowSums(dge$counts)
isexpr <- A > 100 # Keeping genes with total counts more than 100.
dge <- calcNormFactors(dge)
v <- voom(dge[isexpr,], design, plot=FALSE)

# find genes differentially expression between the two groups of samples combined above
fit <- lmFit(v, design) #comparas los modelos de mas de 100 y el de antes
fit <- eBayes(fit)

diff.exp.df <- topTable(fit, coef = "status", n = Inf, sort = "p", p = 0.01) # Positive log-fold-changes mean higher expression in d1
#  hace tabla de genes con le ranking mas alto del modelo lineal
diff.exp.df$gene.name <- rownames(diff.exp.df)

# With the code above, perform the DEA. First, we convert the read counts to log2-cpm, with associated weights, ready for linear modelling. 
# As read counts follow a negative binomial distribution, which has a mathematical theory less tractable than that of the normal distribution, 
# RNAseq data was normalised with the voom methodology. The voom method estimates the mean-variance of the log-counts and generates a precision 
# weight for each observation. This way, a comparative analysis can be performed with all bioinformatic workflows originally developed for microarray
# analyses ([see this](https://bioconductor.org/packages/3.7/bioc/vignettes/CVE/inst/doc/WGCNA_from_TCGA_RNAseq.html) and this paper Charity W Law et al.
# voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. In: Genome biology 15.2 (Jan. 2014), R29–R29).

if(dim(diff.exp.df)[1]<1){
  print("Se para la ejecucion del DEA por falta de genes que cumplan las condiciones establecidas")
}else{
  # Output, Volcano plot
  tab = data.frame(logFC = diff.exp.df$logFC, negLogPval = -log10(diff.exp.df$adj.P.Val))
  tab2 = data.frame(logFC = diff.exp.df$logFC, negLogPval = -log10(diff.exp.df$adj.P.Val), Gene=diff.exp.df$gene.name)
  lfc = 2
  pval = 0.01
  
  HOY  <- format(Sys.time(), "%F_%H.%M.%S")
  nombre_carpeta <- paste(HOY, label, sep = "_" )
  directorio_carpeta_resultados <- paste("Code/Results", nombre_carpeta, sep = "/")
  dir.create(directorio_carpeta_resultados)
  
  write.csv(filter(tab2, abs(logFC) > lfc & negLogPval > -log10(pval)), paste(directorio_carpeta_resultados, "/", label, ".csv", sep = "")) # write output
  
  pdf(file = paste(directorio_carpeta_resultados, "/VolcanoPlot_", label, ".pdf", sep = ""),  width = 9, height = 4.5)
  par(mar = c(5, 4, 4, 5))
  plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue))
  #signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
  points(tab[(abs(tab$logFC) > lfc), ], pch = 16, cex = 0.8, col = "orange") 
  points(tab[(tab$negLogPval > -log10(pval)), ], pch = 16, cex = 0.8, col = "green") 
  points(tab[(abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval)), ], pch = 16, cex = 0.8, col = "red") 
  abline(h = -log10(pval), col = "green3", lty = 2) 
  abline(v = c(-lfc, lfc), col = "blue", lty = 2) 
  mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
  mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
  with(subset(tab2, negLogPval > -log10(pval) & abs(logFC)>lfc), textxy(logFC, negLogPval, labs=Gene, cex=.4))
  dev.off()
  
  
  
}

