## DEA


# Pass arguments to the script
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
        stop("At least one breast cancer type must be supplied (one of Basal, TNBC) ", call.=FALSE)
} else if (length(args)==1) {
        
        args[2] = "LumA"
}

# Load the expression data
load("brca_rnaseq-basal.RData")
load("brca_rnaseq-luminal.RData")
load("brca_rnaseq-tnbc.RData")


source("0_loadLibraries.R")
loadpkg("dplyr")
loadpkg("limma") # Differential gene expression analysis (DEA)
loadpkg("edgeR")
loadpkg("calibrate") # To label the volcano plot

# Pick the datasets to analyse

if(args[1] == "Basal"){d1 = brca_rnaseq.basal}
if(args[1] == "TNBC"){d1 = brca_rnaseq.tnbc}
d2 = brca_rnaseq.luminal

cat("Differential Expression analysis of breast cancer subtypes:",args[1], "vs ",args[2], "\n")





# We Combine the two matrices for gene differential expression analysis (DEA). Further preprocessing included the removal of expression estimates 
# with counts in less than 20% of cases.
rnaseq.for.de <- cbind(d1, d2)
counts = rnaseq.for.de[apply(rnaseq.for.de,1,function(x) sum(x==0))<ncol(rnaseq.for.de)*0.8,]

# Create a design matrix thar contains the RNA samples that are applied to each category (TNBC vs luminal)
df.l <- data_frame("sample" = colnames(d1), "status" = rep(0, length(colnames(d1))) )
df.t <- data_frame("sample" = colnames(d2), "status" = rep(1, length(colnames(d2))) )
df <- rbind(df.t,df.l)
design <- model.matrix(~ status, data = df) #se quedan con 0 las que sean de tnbc y con 1s las que son de luminal


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


# Output, Volcano plot
tab = data.frame(logFC = diff.exp.df$logFC, negLogPval = -log10(diff.exp.df$adj.P.Val))
tab2 = data.frame(logFC = diff.exp.df$logFC, negLogPval = -log10(diff.exp.df$adj.P.Val), Gene=diff.exp.df$gene.name)
lfc = 2
pval = 0.01

write.csv(filter(tab2, abs(logFC) > lfc & negLogPval > -log10(pval)), "dea.csv") # write output

pdf(file = "volcano.pdf", width = 9, height = 4.5)
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