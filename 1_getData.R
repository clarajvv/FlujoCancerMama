## Get the RNASeq and the clinical data


# set seed to be reproductible
set.seed(12345)
options(stringsAsFactors = FALSE)


source("0_loadLibraries.R")

loadpkg("RTCGAToolbox") # To browse and get TCGA data
loadpkg("readxl")
library(RTCGAToolbox)
library(readxl)

#### Get the breast cancer gene expression data from TCGA
brcaData <- getFirehoseData(dataset="BRCA", runDate="20160128",gistic2Date="20160128",forceDownload=F, clinical =TRUE, RNASeq2GeneNorm  =TRUE)

brca_rnaseq <- getData(brcaData,type = "RNASeq2GeneNorm")

#### Gene expression data cleaning up
## keep only the tumour data. Sample type 01-09 in the samples barcode means tumour
brca_rnaseq.tumour <- brca_rnaseq[, which(as.numeric(substr(colnames(brca_rnaseq), 14,15)) < 10)]

## Convert barcodes from the rnaseq dataset to sample barcodes, which identify the patients. This way the sample ID is the same in the brca rnaseq matrix and in the supplementary table 1
colnames(brca_rnaseq.tumour) <- substr(colnames(brca_rnaseq.tumour), 1,12)

## Make sure we don't have duplicated samples
brca_rnaseq.tumour <- brca_rnaseq.tumour[, !duplicated(colnames(brca_rnaseq.tumour))]

### write out the brca tumour rnaseq data matrix
save(brca_rnaseq.tumour, file = "brca_rnaseq.RData")


## Download Supplementary data and import Table 1
## [Supplementary table 1 from Comprehensive molecular portraits of human breast tumours](http://www.nature.com/nature/journal/v490/n7418/full/nature11412.html?foxtrotcallback=true#supplementary-information) has TCGA breast cancer samples classified by PAM50 subtypes and ER Status, PR Status and HER2 Final Status.

temp <- tempfile()
download.file("https://media.nature.com/original/nature-assets/nature/journal/v490/n7418/extref/nature11412-s2.zip",temp)
unzip(temp)
sample_data <- read_excel("nature11412-s2/Supplementary Tables 1-4.xls", sheet = 1, skip = 1)
unlink("nature11412-s2",recursive = T,force = T)
rm(temp)

save(sample_data, file="sample_data.RData")

unlink("20160128*",force = T)
