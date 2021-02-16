#cargar lirberias
source("Code/Settings/0_loadLibraries.R") # no se si hace falta con el bash, creo que no
loadpkg("dplyr")

# Cargar datos
load("Data/brca_rnaseq_tumour.RData")
load("Data/sample_data.RData")

print("Lo primero que haremos será comprobar si el que el cancer tenga receptores de progesterona afecta a la tasa de supervivencia")
d1 <- sample_data[sample_data$`PR Status`== "Positive",]
d2 <- sample_data[sample_data$`PR Status`!= "Positive",]
d2$`PR Status` = "Negative"
datos_limpios <- rbind(d1, d2)
tablaSupervivencia <- xtabs(~datos_limpios$`Vital Status`+ datos_limpios$`PR Status`)
pValor <- unlist(unname(fisher.test(tablaSupervivencia)[1]))

if(pValor<0.05){
  print(paste("El pvalor es ", round(pValor, 5), ", por lo que se ha encontrado relacion entre ambas variables. Procederemos a efectuar un DEA."))
  
  subsettingPR <- function(numNodos){
    samples <- sample_data %>% dplyr::filter(`PR Status` == numNodos)
    barcodes <- samples$`Complete TCGA ID`
    data <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% barcodes)]
    return(data)
  }
  
  brca_rnaseq.PRpos <- subsettingPR("Positive")
  brca_rnaseq.PRneg <- subsettingPR("Negative")
  
  rnaseq.for.de <- cbind(brca_rnaseq.PRpos, brca_rnaseq.PRneg)
  counts <- rnaseq.for.de[apply(rnaseq.for.de,1,function(x) sum(x==0))<ncol(rnaseq.for.de)*0.8,]
  
  df.c1 <- data.frame("sample" = colnames(brca_rnaseq.PRpos), "status" = rep(0, length(colnames(brca_rnaseq.PRpos))) )
  df.c2 <- data.frame("sample" = colnames(brca_rnaseq.PRneg), "status" = rep(1, length(colnames(brca_rnaseq.PRneg))) )
  df <- rbind(df.c1, df.c2)
  
  design <- model.matrix(~ status, data = df)
  label <- "DEA_breast_cancer_PR_status"
  
  preDEA <- list("counts" = counts, "design" = design, "label" = label)
  save(preDEA, file = "Data/datosParaDEA.RData")
} else {
  print(paste("El pvalor es ", round(pValor, 4), ", por lo que no se ha encontrado relacion entre ambas variables. El programa se detendra."))
  quit()
}
