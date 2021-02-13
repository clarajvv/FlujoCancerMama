##borrar
setwd("D:/uma/cuarto/Herramientas y Algoritmos/functional_analysis/flujo cancer")

#cargar lirberias
source("Code/Settings/0_loadLibraries.R") # no se si hace falta con el bash, creo que no
loadpkg("dplyr")

# Cargar datos
load("Data/brca_rnaseq_tumour.RData")
load("Data/sample_data.RData")

print("Lo primero que haremos será comprobar que el que haya nodos afectados tiene relación con la tasa de supervivencia")
tablaSupervivencia <- xtabs(~sample_data$`Vital Status`+ sample_data$`Node-Coded`)
pValor <- unlist(unname(fisher.test(tablaSupervivencia, simulate.p.value = TRUE)[1]))
# si la tiene, no hago un if porque este dato va a ser siempre igual

print(paste("El pvalor es ", round(pValor, 4), ", por lo que se ha encontrado relación entre ambas variables. Procederemos a efectuar un DEA."))

subsettingNodes <- function(numNodos){
  samples <- sample_data %>% dplyr::filter(`Node-Coded` == numNodos)
  barcodes <- samples$`Complete TCGA ID`
  data <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% barcodes)]
  return(data)
}

brca_rnaseq.cant_1 <- subsettingNodes("Positive")
brca_rnaseq.cant_2 <- subsettingNodes("Negative")

rnaseq.for.de <- cbind(brca_rnaseq.cant_1, brca_rnaseq.cant_2)
counts <- rnaseq.for.de[apply(rnaseq.for.de,1,function(x) sum(x==0))<ncol(rnaseq.for.de)*0.8,]

df.c1 <- data.frame("sample" = colnames(brca_rnaseq.cant_1), "status" = rep(0, length(colnames(brca_rnaseq.cant_1))) )
df.c2 <- data.frame("sample" = colnames(brca_rnaseq.cant_2), "status" = rep(1, length(colnames(brca_rnaseq.cant_2))) )
df <- rbind(df.c1, df.c2)

design <- model.matrix(~ status, data = df)
label <- "DEA_breast_cancer_nodes_affected"

preDEA <- list("counts" = counts, "design" = design, "label" = label)
save(preDEA, file = "Data/datosParaDEA.RData")
