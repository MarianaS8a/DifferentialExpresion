if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GO.db")
BiocManager::install("annotate")
#Human libraries
BiocManager::install("org.Hs.eg.db")
BiocManager::install("hugene20stprobeset.db")
BiocManager::install("hugene20sttranscriptcluster.db")
BiocManager::install("affy")
BiocManager::install("limma")


# 1) Leemos los datos en R
# 2) Vemos las dimensiones de la matriz de datos
# 3) Damos un vistazo al tipo de datos que tenemos
getwd()
#setwd("C:/Users/feliy/Desktop/Tec/2do Semestre/Análisis de biología computacional/Analisis de datos expresion diferencial")
CoCaData = read.table(file="ColonData.txt", header=TRUE)
dim(CoCaData)
CoCaData[1:10,1:5]

# Annotation de genes
library(GO.db)
install.packages("XML")
library(annotate)
# Human libraries
library(org.Hs.eg.db)
library(hugene20stprobeset.db)
library(hugene20sttranscriptcluster.db)
annDB = "hugene20sttranscriptcluster.db"
NoTranscripts = 53617
# Librer\’ias de an\’alisis
library(affy)
install.packages("gplots")
library(gplots)
# Librer\’ias de expresi\’on diferencial
library(limma)

# Identificando los genes a traves de los IDs
# ---------------------------------------------
probenames <- as.character(CoCaData$ID_REF)
entID <- getEG(probenames,annDB)
sym <- getSYMBOL(probenames,annDB)
#Crea variable de datos de expresion de genes
#concatena simbolo con matriz de datos
geneExpData <- data.frame(sym, CoCaData)
indx.NAs = !is.na(geneExpData$sym)
#Genera matriz filtrada
geneExpData = geneExpData[indx.NAs,]
dim(geneExpData)
geneExpData[1:10,1:5]

#Generar matriz de datos ordenada con nombre 
#en renglones
datosColon = as.matrix(geneExpData[,3:86])
rownames(datosColon) = geneExpData$sym
datosColon[1:5,1:5]

#Obtener etiquetas del csv
etiquetas = read.csv(file="DetalleMuestrasCoCa.csv", header=TRUE)
#Cambiar nombre de columnas por etiquetas
colnames(datosColon) = as.character(etiquetas$Etiqueta)
#Ordenar columnas de matriz por Normal, tumor, etc...
indx= c(which(substring(etiquetas$Etiqueta,1,6) == "Normal"),
        which(substring(etiquetas$Etiqueta,1,6) == "TumorP"),
        which(substring(etiquetas$Etiqueta,1,3) == "T2D"),
        which(substring(etiquetas$Etiqueta,1,8) == "TumorT2D"))
datosColon = datosColon[,indx]
datosColon[1:5,1:5]
colnames(datosColon)

#Realizar graficas para conocer distribucion y dispersion
colores = c(rep("lightgray", 19), rep("steelblue",19),
            rep("orange",23), rep("indianred",23))
#Compara la densidad para cada gen de cada uno de los pacientes. 
#Cada una de las 84 lineas representa cada paciente y muestra la densidad de cada gen
#Es evidente que la densidad varía entre pacientes
plotDensity(datosColon, col=colores, main = "Datos Normalizados",
            xlab = "Expression Levels")
par(mar=c(7,4,4,4))


#Diagrama caja y bigotes
#Grises normales, azul colon cancer, amarillo diabetes, rojo ambas


#se rompe el programa
#boxplot(datosColon, col=colores, main = "Datos Normalizados",
#        las =2, ylab = "Expression", cex.axis=.7, pch=".")




# Identificando los genes diferencialmente expresados 
# ----------------------------------------------------------- 
# initialize variables 
pv = 2 
M = 0.3 
# design matrix 
smpls = 84 #num de muestras total
gps = 4 #numero de grupos para los datos
design = matrix(rep(0,gps*smpls), nrow=smpls) 
#Nombres de cada uno de los 4 casos
colnames(design) = c("Normal", "CoCa", "T2D", "CoCaT2D") 
design[1:19, 1]=1 
design[20:38,2]=1 
design[39:61,3]=1 
design[62:84,4]=1 
# contrast matrix 
#Genera matriz de contraste entre los casos
cont.matrix = makeContrasts("CoCa - Normal", "CoCaT2D - T2D", "CoCaT2D - CoCa", "T2D - Normal", levels=design) 
# linear model fit 
#Ingresa datos filtrados por orden segun caso a matriz
fit = lmFit(datosColon,design) 
fitC = contrasts.fit(fit, cont.matrix) 
fitCB = eBayes(fitC) 

# resultados comparacion colon cancer vs normal
TT = topTable(fitCB, coef=1, adjust = "fdr", sort.by = "logFC", 
              number=NoTranscripts, genelist=fit$genes) 
genesDifExp = TT[(-log10(TT$P.Value)>pv & abs(TT$logFC)>M),]
genesDifExp
Up.regulated = sum(genesDifExp$logFC>M) 
Up.regulated
Dw.regulated = sum(genesDifExp$logFC<=M)
Dw.regulated
dim(genesDifExp) 
head(genesDifExp) 

#install.packages("xlsx") 
library(xlsx) 
write.xlsx (genesDifExp, file="CoCa vs Normal Genes diferencialmente expresados.xlsx", sheetName="Hoja 1", 
            col.names=TRUE, row.names=TRUE, append=FALSE) 
write.table(genesDifExp, file="CoCa vs Normal Genes diferencialmente expresados.txt", sep = "\t",
            col.names=TRUE, row.names=TRUE, append=FALSE)

#Generar archivos de contraste entre CoCaT2D vs T2D
TT2 = topTable(fitCB, coef=2, adjust = "fdr", sort.by = "logFC", 
               number=NoTranscripts, genelist=fit$genes)
genesDifExp2 = TT2[(-log10(TT2$P.Value)>pv & abs(TT2$logFC)>M),]
write.xlsx (genesDifExp2, file="CoCaT2D vs T2D Genes diferencialmente expresados.xlsx", sheetName="Hoja 1", 
            col.names=TRUE, row.names=TRUE, append=FALSE)
write.table(genesDifExp2, file="CoCaT2D vs T2D Genes diferencialmente expresados.txt", sep = "\t",
            col.names=TRUE, row.names=TRUE, append=FALSE)

#Generar archivos de contraste entre CoCaT2D vs CoCa
TT3 = topTable(fitCB, coef=3, adjust = "fdr", sort.by = "logFC", 
               number=NoTranscripts, genelist=fit$genes)
genesDifExp3 = TT3[(-log10(TT3$P.Value)>pv & abs(TT3$logFC)>M),]
write.xlsx (genesDifExp3, file="CoCaT2D vs CoCa Genes diferencialmente expresados.xlsx", sheetName="Hoja 1", 
            col.names=TRUE, row.names=TRUE, append=FALSE)
write.table(genesDifExp3, file="CoCaT2D vs CoCa Genes diferencialmente expresados.txt", sep = "\t",
            col.names=TRUE, row.names=TRUE, append=FALSE)


#Grafica de volcan
par(mar=c(5,4,4,2)+0.1)
volcanoplot(fitCB, coef=4)
volcanoplot(fitCB, coef="CoCa - Normal", main = "Contraste CoCa vs
Normal")
volcanoplot(fitCB, col = "green", coef="CoCaT2D - T2D", main = "Contraste CoCaT2D vs
T2D")
volcanoplot(fitCB, coef="CoCaT2D - CoCa", main = "Contraste CoCaT2D vs
CoCa")


volcanoplot(fitCB, col="blue", coef=1)
genelabels = genesDifExp$ID
datos.cluster = datosColon[as.numeric(rownames(genesDifExp)),1:38]
par(oma = c(3,1,3,4),mar=c(14,5,2,2)+0.1)
ind.hmap = heatmap.2(datos.cluster, col=greenred(75), scale="row",
                     key=TRUE, symkey=FALSE, density.info="none", trace="none",
                     cexRow=.55, cexCol=.7, main = "CoCa vs Normal", labRow =
                       genelabels, ColSideColors=colores[1:38])

datos.cluster2 = datosColon[as.numeric(rownames(genesDifExp)),39:84]
ind.hmap = heatmap.2(datos.cluster2, col=greenred(75), scale="row",
                     key=TRUE, symkey=FALSE, density.info="none", trace="none",
                     cexRow=.55, cexCol=.7, main = "T2D vs CoCaT2D", labRow =
                       genelabels, ColSideColors=colores[39:84])

coef3 = c(20:38,62:84)
datos.cluster3 = datosColon[as.numeric(rownames(genesDifExp)),coef3]
ind.hmap = heatmap.2(datos.cluster3, col=greenred(75), scale="row",
                     key=TRUE, symkey=FALSE, density.info="none", trace="none",
                     cexRow=.55, cexCol=.7, main = "CoCa vs. CoCaT2D", labRow =
                       genelabels, ColSideColors=colores[coef3])