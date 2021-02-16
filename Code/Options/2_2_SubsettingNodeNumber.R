#cargar lirberias
source("Code/Settings/0_loadLibraries.R") # no se si hace falta con el bash, creo que no
loadpkg("dplyr")

# Cargar datos
load("Data/brca_rnaseq_tumour.RData")
load("Data/sample_data.RData")

print("Lo primero que haremos sera comprobar que la cantidad de nodos afectados tiene relacion con la tasa de supervivencia")
tablaSupervivencia <- xtabs(~sample_data$`Vital Status`+ sample_data$Node)
pValor <- unlist(unname(fisher.test(tablaSupervivencia, simulate.p.value = TRUE)[1]))
# si la tiene, no hago un if porque este dato va a ser siempre igual

print(paste("El pvalor es ", round(pValor, 4), ", por lo que se ha encontrado relacion entre ambas variables. Procederemos a efectuar un DEA."))

args<-c(1, 0, 3)
par_nodos <- args[c(2,3)] #irene esto lo cambias como veas pq no se como se van a pasar los args, asumo que son los de las posiciones segunda y tercera (hablame si eso)

cant_1 <- paste("N", par_nodos[1], sep = "")
cant_2 <- paste("N", par_nodos[2], sep = "")

# samples_cant_1 <- sample_data %>% dplyr::filter(Node == cant_1)
# samples_cant_2 <- sample_data %>% dplyr::filter(Node == cant_2)

# barcodes_nodos_elegidos <- samples_nodos_elegidos$`Complete TCGA ID`
# brca_rnaseq.nodos <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% barcodes_nodos_elegidos)]
# samples_nodos_elegidos <- sample_data %>% dplyr::filter(Node == cant_1 | Node == cant_2)

subsettingNodes <- function(numNodos){
  samples <- sample_data %>% dplyr::filter(Node == numNodos)
  barcodes <- samples$`Complete TCGA ID`
  data <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% barcodes)]
  return(data)
}

brca_rnaseq.cant_1 <- subsettingNodes(cant_1)
brca_rnaseq.cant_2 <- subsettingNodes(cant_2)

rnaseq.for.de <- cbind(brca_rnaseq.cant_1, brca_rnaseq.cant_2)
counts <- rnaseq.for.de[apply(rnaseq.for.de,1,function(x) sum(x==0))<ncol(rnaseq.for.de)*0.8,]

df.c1 <- data.frame("sample" = colnames(brca_rnaseq.cant_1), "status" = rep(0, length(colnames(brca_rnaseq.cant_1))) )
df.c2 <- data.frame("sample" = colnames(brca_rnaseq.cant_2), "status" = rep(1, length(colnames(brca_rnaseq.cant_2))) )
df <- rbind(df.c1, df.c2)

design <- model.matrix(~ status, data = df)
label <- paste("DEA_breast_cancer_affected_nodes", par_nodos[1], "vs", par_nodos[2], sep = "_")

preDEA <- list("counts" = counts, "design" = design, "label" = label)
save(preDEA, file = "Data/datosParaDEA.RData")
