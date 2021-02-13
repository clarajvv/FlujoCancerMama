###########################################################
####### Procesado de datos para un rango de edad ##########
###########################################################

###########################################################
### Libraries  

source("Code/Settings/0_loadLibraries.R") # No sÃ© si es necesario llamar al script o con bash se hace solo
loadpkg("dplyr")

###########################################################
### Loadig data 

load("Data/brca_rnaseq_tumour.RData") # No sÃ© si es la mejor forma de traerme los datos 
load("Data/sample_data.RData")

###########################################################
### Catching arguments and creating subsets

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Quitar el comentario de los args !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#args = commandArgs(trailingOnly = TRUE, asVaule = False)
#Separo el conjunto por el valor introducido. El minimo de edad es 26. El maximo es 90
args <- c(55, 66) #Con este el fisher da menor a 0.05!!!!

onset_age <- function(age){
  if(as.numeric(age) <= args[1]) {return("0")}
  if(as.numeric(age) > args[1] & as.numeric(age) < args[2]) {return(NA)}
  if(as.numeric(age) >= args[2]) {return("1")}
}

df_age_surv <- data.frame(age = as.numeric(sample_data$`Age at Initial Pathologic Diagnosis`), 
                         survival = as.factor(sample_data$`Vital Status`), 
                         barcodes = sample_data$`Complete TCGA ID`)

df_age_surv <- na.omit(df_age_surv)

df_age_surv$age <- sapply(df_age_surv$age, onset_age)

df_age_surv <- na.omit(df_age_surv)

levels(df_age_surv$survival) <- c("0", "1", NA)



###########################################################
### Fisher test

df_fisher <- table(subset(df_age_surv, select = c("age", "survival")))

fisher_result <-  fisher.test(df_fisher)


if(fisher_result$p.value<0.05){
  #####################################################
  ### Preparing data for DEA
      
  first_sample <- sample_data$`Complete TCGA ID`[sample_data$`Age at Initial Pathologic Diagnosis`<args[1]]
  second_sample <- sample_data$`Complete TCGA ID`[sample_data$`Age at Initial Pathologic Diagnosis`>args[1]]
  d1 <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% first_sample)]
  d2 <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% second_sample)]
      
  rnaseq.for.de <- cbind(d1, d2)
  counts <-  rnaseq.for.de[apply(rnaseq.for.de, 1, function(x) sum(x==0)) < ncol(rnaseq.for.de)*0.8, ]
  # Create a design matrix thar contains the RNA samples that are applied to each category (TNBC vs luminal)
  df.l <- data_frame("sample" = colnames(d1), "status" = rep(0, length(colnames(d1))) )
  df.t <- data_frame("sample" = colnames(d2), "status" = rep(1, length(colnames(d2))) )
  df <- rbind(df.t,df.l)
  design <- model.matrix(~ status, data = df) #se quedan con 0 las que sean de tnbc y con 1s las que son de luminal
      
  # We have to return counts, design and etiqueta para el DEA
  label <- paste("DEA_breast_cancer_age_menor", args[1], "_vs_mayor_", args[2], sep = "")
  
  preDEA <- list("counts" = counts, "design" = design, "label" = label)
    
  save(preDEA, file = "Data/datosParaDEA.RData")
    
  print("Se procede a hacer DEA del estudio edad-supervivencia")
}else{
  preDEA <- NULL
  save(preDEA, file = "Data/datosParaDEA.RData")
  print("Segun el test de fisher, no hay correlacion entre los intervalos de edad introducidos y la supervivencia")
}


