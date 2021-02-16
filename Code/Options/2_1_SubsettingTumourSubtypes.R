
###########################################################
#### Procesado de datos de los subtipos de tumor ##########
###########################################################

###########################################################
### Libraries  

source("Code/Settings/0_loadLibraries.R") 
loadpkg("dplyr")


###########################################################
### Loading data 

load("Data/brca_rnaseq_tumour.RData") 
load("Data/sample_data.RData")


###########################################################
### Functions 


subsettingTNBC <- function(){
  tnbc_samples <- sample_data %>% dplyr::filter(`ER Status` == "Negative" & `PR Status` == "Negative" & `HER2 Final Status` == "Negative" & `PAM50 mRNA` != "Luminal A")
  tnbc_barcodes <- tnbc_samples$`Complete TCGA ID`
  brca_rnaseq.tnbc <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% tnbc_barcodes)]
  return(brca_rnaseq.tnbc)
}

subsettingLumA <- function(){
  luminal_samples <- sample_data %>% dplyr::filter(`PAM50 mRNA` == "Luminal A")
  luminal_barcodes <- luminal_samples$`Complete TCGA ID`
  brca_rnaseq.luminal <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% luminal_barcodes)]
  return(brca_rnaseq.luminal)
}

subsettingBasal <- function(){
  basal_samples <- sample_data %>% dplyr::filter(`PAM50 mRNA` == "Basal-like")
  basal_barcodes <- basal_samples$`Complete TCGA ID`
  brca_rnaseq.basal <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% basal_barcodes)]
  return(brca_rnaseq.basal)
}


# We Combine the two matrices for gene differential expression analysis (DEA). Further preprocessing included the removal of expression estimates 
# with counts in less than 20% of cases.

matricesCombinationDEA <- function(d1, d2, nombreSubtipo1, nombreSubtipo2){
  rnaseq.for.de <- cbind(d1, d2)
  counts <-  rnaseq.for.de[apply(rnaseq.for.de, 1, function(x) sum(x==0)) < ncol(rnaseq.for.de)*0.8, ]
  # Create a design matrix thar contains the RNA samples that are applied to each category (TNBC vs luminal)
  df.l <- data_frame("sample" = colnames(d1), "status" = rep(0, length(colnames(d1))) )
  df.t <- data_frame("sample" = colnames(d2), "status" = rep(1, length(colnames(d2))) )
  df <- rbind(df.t,df.l)
  design <- model.matrix(~ status, data = df) #se quedan con 0 las que sean de tnbc y con 1s las que son de luminal
  
  # We have to return counts, design and etiqueta para el DEA
  label <- paste("DEA_breast_cancer_subtypes_", nombreSubtipo1, "_vs_", nombreSubtipo2, sep = "")

  return(list("counts" = counts, "design" = design, "label" = label))
}



#############################################################
### Checking entry parameter 
### Code: Basal -> 1. LumA -> 2. TNBC -> 3

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Corregir !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Preguntar por consola
#args <- c(2, 3)
args = commandArgs(trailingOnly = TRUE)
                                 
if (length(args) < 2) {
  stop("Two breast cancer type must be supplied (Basal, LumA or TNBC) ", call.=FALSE)
}else if((args[1] == 1 & args[2] == 2) | (args[1] == 2 & args[2] == 1)){
  # Basal - LuminalA
  if(!file.exists("Data/preDeaBasalLumA.RData")){
    subBasal <- subsettingBasal()
    subLumA <-  subsettingLumA()
    preDEA <- matricesCombinationDEA(d1 = subBasal, d2 = subLumA, nombreSubtipo1 = "Basal", nombreSubtipo2 =  "LuminalA")
    save(preDEA, file = "Data/preDeaBasalLumA.RData")
  }
  
  file_datos_DEA <- "Data/preDeaBasalLumA.RData"
  
}else if((args[1] == 1 & args[2] == 3) | (args[1] == 3 & args[2] == 1)){
  # Basal - TNBC -- Problemaaa
  if(!file.exists("Data/preDeaBasalTNBC.RData")){
    subBasal <- subsettingBasal()
    subTNBC <-  subsettingTNBC()
    preDEA <- matricesCombinationDEA(d1 = subBasal, d2 = subTNBC, nombreSubtipo1 = "Basal", nombreSubtipo2 =  "TNBC")
    save(preDEA, file = "Data/preDeaBasalTNBC.RData")
  }
  file_datos_DEA <- "Data/preDeaBasalTNBC.RData"
  
}else if((args[1] == 2 & args[2] == 3) | (args[1] == 3 & args[2] == 2)){
  # LumA - TNBC
  if(!file.exists("Data/preDeaLumATNBC.RData")){
    subLumA <- subsettingLumA()
    subTNBC <-  subsettingTNBC()
    preDEA <- matricesCombinationDEA(d1 = subLumA, d2 = subTNBC, nombreSubtipo1 = "LuminalA", nombreSubtipo2 =  "TNBC")
    save(preDEA, file = "Data/preDeaLumATNBC.RData")
  }
  file_datos_DEA <- "Data/preDeaLumATNBC.RData"
}else{
  stop("It must be two argument with the following code: Basal-> 1, LumA-> 2 and TNBC-> 3")
}

###########################################################
## Pasar al DEA count y design de la combinación solicitada. 
load(file_datos_DEA)
save(preDEA, file = "Data/datosParaDEA.RData")

