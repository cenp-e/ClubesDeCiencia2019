
install.packages("devtools")
library(devtools)
install_github("satijalab/seurat", ref = REF)
library(Seurat)
install.packages("devtools")
library(devtools)
library(useful)

library(dplyr)
#mover datos dentro Rstudio
human.data=read.table("HumanDataSet/HumanSingleCellData.txt",sep="\t",header=TRUE,row.names=1)

#para ver una esquina de los datos
corner(human.data)


#crear seurat object para seguir analizando
# tomar todos los genes en > 2 cells, todas las celulas con > 50 genes
human <- CreateSeuratObject(raw.data = human.data, min.cells = 2, min.genes = 50, project = "Human")
human <- NormalizeData(object = human, normalization.method = "LogNormalize", scale.factor = 10000)

#encontrar genes variables
human <- FindVariableGenes(object = human, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
human <- ScaleData(object = human)
human <- RunPCA(object = human, pc.genes = human@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
human<- JackStraw(object = human, num.replicate = 100)


#agrupar las celulas y crear TSNEPlot
human<- FindClusters(object = human, reduction.type = "pca", dims.use = 1:12, 
                     resolution = 0.9, print.output = 0, save.SNN = TRUE)
human<- RunTSNE(object = human, dims.use = 1:12, do.fast = TRUE)
TSNEPlot(object = human)
#guarda la gráfica y datos
saveRDS(human, file = "~/human_tutorial.rds")

#Encuentra marcadores para cada cluster.
cluster1.markers <- FindMarkers(object = human, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))
human.markers <- FindAllMarkers(object = human, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)
human.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

#visualizar
VlnPlot(object = human, features.plot = c("SOX2", "PAX6", "HOPX"))
FeaturePlot(object = human, features.plot = c("SOX2", "SATB2", "PAX6", "HOPX"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

#renombrar después de la identificación
current.cluster.ids <- c(0, 1, 2)
new.cluster.ids <- c("X0", "X1", "X2")
human@ident <- plyr::mapvalues(x = human@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = human, do.label = TRUE, pt.size = 0.5)









