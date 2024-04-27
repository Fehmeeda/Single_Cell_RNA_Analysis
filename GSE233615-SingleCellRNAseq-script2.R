# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(SeuratObject)
library(patchwork)
library(celldex)
library(openxlsx)
library(clusterProfiler)
library(enrichplot)
library(pathview)


# first load normal samples from a dataset and create Seurat Object as well.

options(Seurat.object.assay.version = "v3")

#sample1
sample1 <- Read10X_h5(filename = "GSE118127-GSE233615_RAW-dataset/GSE118127-Normalsamplesdataset/GSM3319043_sample_3-16_filtered_gene_bc_matrices_h5.h5")
# Initialize the Seurat object with the raw (non-normalized data).
sample1 <- CreateSeuratObject(counts = sample1, project = "Normal", min.cells = 3, min.features = 200)

#sample2
sample2 <- Read10X_h5(filename = "GSE118127-GSE233615_RAW-dataset/GSE118127-Normalsamplesdataset/GSM3319044_sample_3-17_filtered_gene_bc_matrices_h5.h5")
# Initialize the Seurat object with the raw (non-normalized data).
sample2 <- CreateSeuratObject(counts = sample2, project = "Normal", min.cells = 3, min.features = 200)

#sample3
sample3 <- Read10X_h5(filename = "GSE118127-GSE233615_RAW-dataset/GSE118127-Normalsamplesdataset/GSM3557969_sample6a_B1_i12H_filtered_gene_bc_matrices_h5.h5")
# Initialize the Seurat object with the raw (non-normalized data).
sample3 <- CreateSeuratObject(counts = sample3, project = "Normal", min.cells = 3, min.features = 200)

#sample4
sample4 <- Read10X_h5(filename = "GSE118127-GSE233615_RAW-dataset/GSE118127-Normalsamplesdataset/GSM3557961_sample11_B2_i10F_filtered_gene_bc_matrices_h5.h5")
# Initialize the Seurat object with the raw (non-normalized data).
sample4 <- CreateSeuratObject(counts = sample4, project = "Normal", min.cells = 3, min.features = 200)

# Now load ovarian samples and also create seurat objects:

#sample1
sample5<- Read10X_h5(filename = "GSE118127-GSE233615_RAW-dataset/OvarianCancerdataset-GSE233615_RAW/GSM7431434_E1_filtered_feature_bc_matrix.h5")
# Initialize the Seurat object with the raw (non-normalized data).
sample5 <- CreateSeuratObject(counts = sample5, project = "Ovarian", min.cells = 3, min.features = 200)

#sample2
sample6<- Read10X_h5(filename = "GSE118127-GSE233615_RAW-dataset/OvarianCancerdataset-GSE233615_RAW/GSM7431435_E2_filtered_feature_bc_matrix.h5")
# Initialize the Seurat object with the raw (non-normalized data).
sample6 <- CreateSeuratObject(counts = sample6, project = "Ovarian", min.cells = 3, min.features = 200)

#sample3
sample7 <- Read10X_h5(filename = "GSE118127-GSE233615_RAW-dataset/OvarianCancerdataset-GSE233615_RAW/GSM7431436_E3_filtered_feature_bc_matrix.h5")
# Initialize the Seurat object with the raw (non-normalized data).
sample7 <- CreateSeuratObject(counts = sample7, project = "Ovarian", min.cells = 3, min.features = 200)

#sample4
sample8 <- Read10X_h5(filename = "GSE118127-GSE233615_RAW-dataset/OvarianCancerdataset-GSE233615_RAW/GSM7431437_E4_filtered_feature_bc_matrix.h5")
# Initialize the Seurat object with the raw (non-normalized data).
sample8 <- CreateSeuratObject(counts = sample8, project = "Ovarian", min.cells = 3, min.features = 200)

# Create a vector of cell IDs
cell_ids <- paste0("sample", 1:8)  
# Merge Seurat Objects Normal
Mergedsamples<- merge(sample1, y = c(sample2, sample3,sample4,sample5,sample6,sample7,sample8),          
                      add.cell.ids = cell_ids,
                      project = 'MergedSamplesNormalandOvarian')
#view(Mergedsamples@meta.data)


# Quality control and selecting cells for further analysis
# Data normalization 
# Identification of highly variable features (feature selection)
# Data scaling

# view QC metrics

range(Mergedsamples$nFeature_RNA)

range(Mergedsamples$nCount_RNA)

# store mitochondrial percentage in object meta data
Mergedsamples <- PercentageFeatureSet(Mergedsamples, pattern = "^MT-", col.name = "percent.mt")

#view(Mergedsamples@meta.data)

range(Mergedsamples$percent.mt)

##### Selecting cells for further analysis
# Visualize QC metrics as a violin plot
print("Selecting cells for further analysis")
VlnPlot(Mergedsamples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Mergedsamples <- subset(Mergedsamples, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)

VlnPlot(Mergedsamples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Data normalization (NormalizeData)
# Feature selection/Identification of highly variable features (FindVariableFeatures)
# Data scaling (ScaleData)

#View(Mergedsamples)

# Data normalization
print("Data normalization")
#Mergedsamples@assays[["RNA"]]@data@x

Mergedsamples <- NormalizeData(Mergedsamples, 
                               normalization.method ="LogNormalize", scale.factor = 10000)

#Mergedsamples@assays[["RNA"]]@data@x

# Identification of highly variable features (feature selection)
print(" feature selection")
Mergedsamples

Mergedsamples <- FindVariableFeatures(Mergedsamples, 
                                      selection.method = "vst", nfeatures = 2000)

# Top features and featurePlot
print("Top features and featurePlot")

variable_genes<-Mergedsamples@assays[["RNA"]]@var.features
write.table(variable_genes, file = "variable_genes.txt", sep = "\t", row.names = FALSE)

VariableFeaturePlot(Mergedsamples)

Mergedsamples

# Data scaling
print("Data scaling")
Mergedsamples <- ScaleData(Mergedsamples) #2000 identified variable features 

all.genes <- rownames(Mergedsamples)
Mergedsamples <- ScaleData(Mergedsamples, features = all.genes) 

#Mergedsamples@assays[["RNA"]]@data@x

#View(Mergedsamples@commands)

#Perform PCA on the scaled data (linear dimensional reduction)
print("linear dimensional reduction")
Mergedsamples <- RunPCA(Mergedsamples, features = VariableFeatures(object = Mergedsamples))

# Examine and visualize PCA results a few different ways
print("PCA")
# DimPlot(), VizDimReduction() and DimHeatmap ()
#DimPlot(Mergedsamples, reduction = "pca", dims = c(1,2))
DimPlot(Mergedsamples, reduction = "pca", dims = c(1, 10))
#DimPlot(Mergedsamples, reduction = "pca", dims = c(1, 50))

# Determine the ‘dimensional’ of the dataset
print("dimensional")
Mergedsamples <- JackStraw(Mergedsamples, num.replicate = 100)
Mergedsamples <- ScoreJackStraw(Mergedsamples, dims = 1:20)
JackStrawPlot(Mergedsamples, dims = 1:20)

#ElbowPlot(Mergedsamples)
ElbowPlot(Mergedsamples, ndims = 50, reduction = "pca")

# Cluster the cells
print("Cluster the cells")
Mergedsamples <- FindNeighbors(Mergedsamples, dims = 1:20)
Mergedsamples <- FindClusters(Mergedsamples, resolution=c(0.1,0.3,0.5,0.7,1))

# Run non-linear dimensional reduction (UMAP/tSNE)
print("non-linear dimensional reduction ")
Mergedsamples <- RunUMAP(Mergedsamples, dims = 1:20)
#DimPlot(Mergedsamples, reduction = "umap", label = TRUE, repel = TRUE)

Mergedsamples <- RunTSNE(object = Mergedsamples)
#DimPlot(object = Mergedsamples, reduction = "tsne")

Mergedsamples <- RunTSNE(object = Mergedsamples)
#DimPlot(object = Mergedsamples, reduction = "tsne")
DimPlot(object = Mergedsamples, reduction = "tsne", group.by = 'orig.ident')

######Setup the Seurat objects
# split the object into a list of two repeats
print("Setup the Seurat objects")
Mergedsamples.list<-SplitObject(Mergedsamples, split.by = 'orig.ident')

#list()

# normalize and identify variable features for each dataset independently 
print(" normalize and identify variable features for each dataset independently ")
Mergedsamples.list <- lapply(Mergedsamples.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  
})

#select features that are repeatedly variable across dataset for integration
print("select features that are repeatedly variable across dataset for integration")

features<-SelectIntegrationFeatures(object.list = Mergedsamples.list)

#Find integration Anchors
print("Find integration Anchors")
Mergedsamples.anchors<-FindIntegrationAnchors(object.list = Mergedsamples.list, anchor.features = features)

#perform integration to create an 'integrated' data assay
print("integrated")
Mergedsamples.integrated<-IntegrateData(anchorset = Mergedsamples.anchors)

###Perfrom an integrated analysis
print("Perfrom an integrated analysis")
#specify that we will perform downstream analysis on the integrated data
DefaultAssay(Mergedsamples.integrated)<-"integrated"

#run the standard workflow for visulaization and clustering
print("run the standard workflow for visulaization and clustering")

Mergedsamples.integrated<-ScaleData(Mergedsamples.integrated,verbose = FALSE)
Mergedsamples.integrated<-RunPCA(Mergedsamples.integrated, npcs = 50, verbose = FALSE)
Mergedsamples.integrated<-FindNeighbors(Mergedsamples.integrated, 
                                        reduction="pca", dims=1:30)
Mergedsamples.integrated<-FindClusters(Mergedsamples.integrated,resolution=0.3)
Mergedsamples.integrated<-RunUMAP(Mergedsamples.integrated, reduction="pca",
                                  dims=1:30)

#Visualization
print("Visualization")

#DimPlot(Mergedsamples.integrated,reduction="umap",label = TRUE)

#Compare
print("Compare")
plot1<-DimPlot(Mergedsamples,reduction="umap", group.by = 'orig.ident')
plot2<-DimPlot(Mergedsamples.integrated,reduction="umap", group.by = 'orig.ident')
plot1+plot2

ref<-celldex::HumanPrimaryCellAtlasData()
#View(as.data.frame(colData(ref)))

pbmc_counts<-GetAssayData(Mergedsamples.integrated)
library(SingleR)
pred<-SingleR(test=pbmc_counts,
              ref=ref,
              labels=ref$label.main)

pred

Mergedsamples.integrated$singleR.labels<-pred$labels[match(rownames(Mergedsamples.integrated@meta.data),rownames(pred))]
#DimPlot(Mergedsamples.integrated, reduction = "umap", group.by = "singleR.labels")+NoLegend()


plotScoreHeatmap(pred)

plotDeltaDistribution(pred)
library(pheatmap)
tab <- table(Assigned=pred$labels, Clusters=Mergedsamples.integrated$seurat_clusters)
pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))

# setting Idents as Seurat annotations provided (also a sanity check!)
print("setting Idents as Seurat annotations provided")

Idents(Mergedsamples.integrated) <- Mergedsamples.integrated@meta.data$singleR.labels
Idents(Mergedsamples.integrated)

#DimPlot(Mergedsamples.integrated, reduction = 'umap', label = TRUE)


# findMarkers between conditions ---------------------
print("findMarkers between conditions")

Mergedsamples.integrated$celltype.cnd <- paste0(Mergedsamples.integrated$singleR.labels,'_', Mergedsamples.integrated$orig.ident)
#View(Mergedsamples.integrated@meta.data)
Idents(Mergedsamples.integrated) <- Mergedsamples.integrated$celltype.cnd

#DimPlot(Mergedsamples.integrated, reduction = 'umap', label = TRUE)


all_Marker<-FindAllMarkers(Mergedsamples.integrated,
                           logfc.threshold = 0.25,
                           min.pct = 0.1,
                           only.pos = TRUE)

# Add a column indicating gene status
print("Add a column indicating gene status")

all_Marker$status <- ifelse(all_Marker$avg_log2FC >= 2 & all_Marker$p_val_adj < 0.0001, "Upregulated",
                            ifelse(all_Marker$avg_log2FC < 2 & all_Marker$p_val_adj < 0.0001, "Downregulated",
                                   "Not significant"))
library(openxlsx)
write.xlsx(all_Marker,"singlecellDEGs.xlsx",rowNames=FALSE)
print("DONE")
