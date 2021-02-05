source("0_loadLibraries.R")
loadpkg("dplyr")

## Load data from module 1

load("brca_rnaseq.RData")
load("sample_data.RData")


#### Get the breast cancer subtypes data


tnbc_samples <- sample_data %>% dplyr::filter(`ER Status` == "Negative" & `PR Status` == "Negative" & `HER2 Final Status` == "Negative" & `PAM50 mRNA` != "Luminal A")
tnbc_barcodes <- tnbc_samples$`Complete TCGA ID`

luminal_samples <- sample_data %>% dplyr::filter(`PAM50 mRNA` == "Luminal A")
luminal_barcodes <- luminal_samples$`Complete TCGA ID`

basal_samples <- sample_data %>% dplyr::filter(`PAM50 mRNA` == "Basal-like")
basal_barcodes <- basal_samples$`Complete TCGA ID`


#### Subset the expression data according to cancer subtypes

brca_rnaseq.tnbc <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% tnbc_barcodes)]

brca_rnaseq.luminal <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% luminal_barcodes)]

brca_rnaseq.basal <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% basal_barcodes)]

save(brca_rnaseq.basal, file = "brca_rnaseq-basal.RData")
save(brca_rnaseq.luminal, file = "brca_rnaseq-luminal.RData")
save(brca_rnaseq.tnbc, file = "brca_rnaseq-tnbc.RData")
