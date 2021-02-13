
###########################################################
#### Procesado de datos de los subtipos de tumor ##########
###########################################################

###########################################################
### Libraries  

source("Code/Settings/0_loadLibraries.R") # No s√© si es necesario llamar al script o con bash se hace solo
loadpkg("dplyr")


###########################################################
### Loadig data 

load("Data/brca_rnaseq_tumour.RData") # No s√© si es la mejor forma de traerme los datos 
load("Data/sample_data.RData")


###########################################################
### Functions 


#### Get the breast cancer subtypes data


###########################################################
# Prepocesado de datos

onset_tumor <- function(grade){
  if(grade == "T1") {return("1")}
  if(grade == "T2") {return("2")}
  if(grade == "T3") {return("3")}
  if(grade == "TX" || grade == "NA" || grade == "T4") {return(NA)}
}

onset_vital <- function(vital){
  if(vital == "DECEASED") {return("0")}
  if(vital == "LIVING") {return("1")}
}


tumor_predata <- data.frame(ID = sample_data$`Complete TCGA ID`, tumor = as.factor(sample_data$Tumor), 
                          vital = as.factor(sample_data$`Vital Status`))

tumor_predata <- na.omit(tumor_predata)

tumor_predata$tumor <- sapply(tumor_predata$tumor, onset_tumor)
levels(tumor_predata$vital) <- c("0", "1", NA)

total <- na.omit(tumor_predata)

#####################################################
###  fisher test

tb <- table(total[,-1])

fis<- fisher.test(tb, workspace = 2e8)

if(fis$p.value<0.05){
#####################################################
### Preparing data for DEA
  
  first_sample <- sample_data %>% dplyr::filter(Tumor=="T1")
  second_sample <- sample_data %>% dplyr::filter(Tumor=="T2")
  third_sample <- sample_data %>% dplyr::filter(Tumor=="T3")
  
  first_barcodes <- first_sample$`Complete TCGA ID`
  second_barcodes <- second_sample$`Complete TCGA ID`
  third_barcodes <- third_sample$`Complete TCGA ID`
  
  d1 <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% first_barcodes)]
  d2 <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% second_barcodes)]
  d3 <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% third_barcodes)]
  
  rnaseq.for.de <- cbind(d1, d2, d3)
  counts <-  rnaseq.for.de[apply(rnaseq.for.de, 1, function(x) sum(x==0)) < ncol(rnaseq.for.de)*0.8, ]
  
  # Create a design matrix thar contains the RNA samples that are applied to each category (TNBC vs luminal)
  df.f <- data_frame("sample" = colnames(d1), "status" = rep(0, length(colnames(d1))) )
  df.s <- data_frame("sample" = colnames(d2), "status" = rep(1, length(colnames(d2))) )
  df.t <- data_frame("sample" = colnames(d3), "status" = rep(3, length(colnames(d3))) )
  df <- rbind(df.f,df.s, df.t)
  design <- model.matrix(~ status, data = df) #se quedan con 0 las que sean de tnbc y con 1s las que son de luminal
  
  # We have to return counts, design and etiqueta para el DEA
  
  preDEA <- list("counts" = counts, "design" = design, "label" = "DEA_breast_cancer_grade")
  
  save(preDEA, file = "Data/datosParaDEA.RData")
  
  print("A continuaciÛn se realizar· el DEA del estudio grado-supervivencia")
}else{
  preDEA <- NULL
  save(preDEA, file = "Data/datosParaDEA.RData")
  print("Segun el test de fisher, no hay correlacion entre los intervalos de grade introducidos y la supervivencia")
}

