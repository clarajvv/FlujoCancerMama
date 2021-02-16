###############################################
#### Get the RNASeq and the clinical data #####
###############################################

# En bash se puede establecer el directorio de trabajo????

##############################################################
## Libraries 
source("Code/Settings/0_loadLibraries.R")

loadpkg("RTCGAToolbox") # To browse and get TCGA data
loadpkg("readxl")


###############################################################
#### Get the breast cancer gene expression data from TCGA
load("Data/brca_rnaseq.RData")


###############################################################
#### Gene expression data cleaning up
## keep only the tumour data. Sample type 01-09 in the samples barcode means tumour
brca_rnaseq.tumour <- brca_rnaseq[, which(as.numeric(substr(colnames(brca_rnaseq), 14,15)) < 10)]

## Convert barcodes from the rnaseq dataset to sample barcodes, which identify the patients. This way the sample ID is the same in the brca rnaseq matrix and in the supplementary table 1
colnames(brca_rnaseq.tumour) <- substr(colnames(brca_rnaseq.tumour), 1,12)

## Make sure we don't have duplicated samples
brca_rnaseq.tumour <- brca_rnaseq.tumour[, !duplicated(colnames(brca_rnaseq.tumour))]



##################################################################
## Download Supplementary data and import Table 1
## [Supplementary table 1 from Comprehensive molecular portraits of human breast tumours](http://www.nature.com/nature/journal/v490/n7418/full/nature11412.html?foxtrotcallback=true#supplementary-information) has TCGA breast cancer samples classified by PAM50 subtypes and ER Status, PR Status and HER2 Final Status.


sample_data <- read_excel("Data/Supplementary Tables 1-4.xls", sheet = 1, skip = 1)


####################################################################
### Save data
save(brca_rnaseq.tumour, file = "Data/brca_rnaseq_tumour.RData")
save(sample_data, file="Data/sample_data.RData")

