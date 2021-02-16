###########################################################
#### Procesado de datos de los grados de tumor ##########
###########################################################

###########################################################
### Libraries  

source("Code/Settings/0_loadLibraries.R") # No se si es necesario llamar al script o con bash se hace solo
loadpkg("dplyr")


###########################################################
### Loadig data 

load("Data/brca_rnaseq_tumour.RData") # No se si es la mejor forma de traerme los datos 
load("Data/sample_data.RData")


###########################################################
### Functions 

onset_tumor <- function(grade){
  if(grade == "T1") {return("T1")}
  if(grade == "T2") {return("T2")}
  if(grade == "T3"|| grade == "T4") {return("T3")}
  if(grade == "TX" || grade == "NA" ) {return(NA)}
}

###########################################################
# Prepocesado de datos


tumor_predata <- data.frame(ID = sample_data$`Complete TCGA ID`, tumor = as.factor(sample_data$Tumor), 
                            vital = as.factor(sample_data$`Vital Status`))

tumor_predata$tumor <- sapply(tumor_predata$tumor, onset_tumor)

tumor_predata <- na.omit(tumor_predata)

print("Lo primero que haremos sera comprobar que la cantidad de grados de tumor tiene relacion con la tasa de supervivencia")
tablaSupervivencia <- xtabs(~tumor_predata$vital+ tumor_predata$tumor)
pValor <- unlist(unname(fisher.test(tablaSupervivencia, simulate.p.value = TRUE)[1]))

print(paste("El pvalor es ", round(pValor, 4), ", por lo que se ha encontrado relacion entre ambas variables. Procederemos a efectuar un DEA."))


if(pValor<0.05){
  #####################################################
  ### Preparing data for DEA
  
  first_sample <- sample_data %>% dplyr::filter(Tumor=="T1")
  third_sample <- sample_data %>% dplyr::filter(Tumor=="T3")
  
  first_barcodes <- first_sample$`Complete TCGA ID`
  third_barcodes <- third_sample$`Complete TCGA ID`
  
  d1 <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% first_barcodes)]
  d3 <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% third_barcodes)]
  
  rnaseq.for.de <- cbind(d1, d3)
  counts <-  rnaseq.for.de[apply(rnaseq.for.de, 1, function(x) sum(x==0)) < ncol(rnaseq.for.de)*0.8, ]
  
  # Create a design matrix thar contains the RNA samples that are applied to each category (TNBC vs luminal)
  df.f <- tibble("sample" = colnames(d1), "status" = rep(0, length(colnames(d1))) )
  df.t <- tibble("sample" = colnames(d3), "status" = rep(1, length(colnames(d3))) )
  df <- rbind(df.f, df.t)
  design <- model.matrix(~ status, data = df) #se quedan con 0 las que sean de tnbc y con 1s las que son de luminal
  
  # We have to return counts, design and etiqueta para el DEA
  
  preDEA <- list("counts" = counts, "design" = design, "label" = "DEA_breast_cancer_grade")
  
  save(preDEA, file = "Data/datosParaDEA.RData")
  
  print("A continuación se realizara el DEA del estudio grado-supervivencia")
}else{
  preDEA <- NULL
  save(preDEA, file = "Data/datosParaDEA.RData")
  print("Segun el test de fisher, no hay correlacion entre los intervalos de grade introducidos y la supervivencia")
}

