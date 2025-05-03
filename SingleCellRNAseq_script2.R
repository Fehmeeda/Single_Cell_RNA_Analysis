# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(SeuratObject)
library(patchwork)
library(celldex)
library(openxlsx)
library(SingleR)


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
sample5<- Read10X_h5(filename = "GSE118127-GSE233615_RAW-dataset/GSE118127-Normalsamplesdataset/GSM7431434_E1_filtered_feature_bc_matrix.h5")
# Initialize the Seurat object with the raw (non-normalized data).
sample5 <- CreateSeuratObject(counts = sample5, project = "Ovarian", min.cells = 3, min.features = 200)

#sample2
sample6<- Read10X_h5(filename = "GSE118127-GSE233615_RAW-dataset/GSE118127-Normalsamplesdataset/GSM7431435_E2_filtered_feature_bc_matrix.h5")
# Initialize the Seurat object with the raw (non-normalized data).
sample6 <- CreateSeuratObject(counts = sample6, project = "Ovarian", min.cells = 3, min.features = 200)

#sample3
sample7 <- Read10X_h5(filename = "GSE118127-GSE233615_RAW-dataset/GSE118127-Normalsamplesdataset/GSM7431436_E3_filtered_feature_bc_matrix.h5")
# Initialize the Seurat object with the raw (non-normalized data).
sample7 <- CreateSeuratObject(counts = sample7, project = "Ovarian", min.cells = 3, min.features = 200)

#sample4
sample8 <- Read10X_h5(filename = "GSE118127-GSE233615_RAW-dataset/GSE118127-Normalsamplesdataset/GSM7431437_E4_filtered_feature_bc_matrix.h5")
# Initialize the Seurat object with the raw (non-normalized data).
sample8 <- CreateSeuratObject(counts = sample8, project = "Ovarian", min.cells = 3, min.features = 200)

# Create a vector of cell IDs
cell_ids <- paste0("sample", 1:8)  
# Merge Seurat Objects Normal
merged_samples<- merge(sample1, y = c(sample2, sample3,sample4,sample5,sample6,sample7,sample8),          
                      add.cell.ids = cell_ids,
                      project = 'MergedSamplesNormalandOvarian')
view(merged_samples@meta.data)

total_cells <- ncol(merged_samples@assays$RNA@data)
print(paste("Total number of cells NormaL AND Ovarian:", total_cells))

# Quality control and selecting cells for further analysis
# Data normalization 
# Identification of highly variable features (feature selection)
# Data scaling

# view QC metrics

# Perform quality control
merged_samples <- PercentageFeatureSet(merged_samples, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(merged_samples, features = c("nFeature_RNA", "nCount_RNA" , "percent.mt") ,  ncol = 3 )
merged_samples <- subset(merged_samples, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 &  nCount_RNA < 25000 & percent.mt <= 25)
VlnPlot(merged_samples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


total_cells <- ncol(merged_samples@assays$RNA@data)
print(paste("Total number of cells NormaL AND Ovarian:", total_cells))

# Normalize data and find variable features
merged_samples <- NormalizeData(merged_samples, normalization.method = "LogNormalize", scale.factor = 10000)
merged_samples <- FindVariableFeatures(merged_samples, selection.method = "vst", nfeatures = 2000)

# Top features and featurePlot

variable_genes<-merged_samples@assays[["RNA"]]@var.features
write.table(variable_genes, file = "variable_genes.txt", sep = "\t", row.names = FALSE)

#VariableFeaturePlot(Mergedsamples)
#Identify the 10 most highly variable genes:

top10 <- head(VariableFeatures(merged_samples), 10)
plot1 <- VariableFeaturePlot(merged_samples)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

# Set up for integration

###Setup the Seurat objects
# split the object into a list of two repeats
samples_list<-SplitObject(merged_samples, split.by = 'orig.ident')

# Select integration features and find anchors
features <- SelectIntegrationFeatures(object.list = samples_list)
anchors <- FindIntegrationAnchors(object.list = samples_list, anchor.features = features)
integrated <- IntegrateData(anchorset = anchors)

# Perform integrated analysis
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)

# Visualize PCA results
DimPlot(integrated, reduction = "pca", dims = c(1, 10))
DimPlot(integrated, reduction = "pca", dims = c(1, 30))
DimPlot(integrated, reduction = "pca", dims = c(1, 50))

# Determine the ‘dimensional’ of the dataset
Mergedsamples <- JackStraw(integrated, num.replicate = 100)
Mergedsamples <- ScoreJackStraw(Mergedsamples, dims = 1:20)
JackStrawPlot(Mergedsamples, dims = 1:20)

ElbowPlot(integrated)
ElbowPlot(integrated, ndims = 50, reduction = "pca")


# Perform clustering and UMAP
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.3)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)

# Annotation using SingleR
ref <- celldex::HumanPrimaryCellAtlasData()
pbmc_counts <- GetAssayData(integrated)
pred <- SingleR(test = pbmc_counts, ref = ref, labels = ref$label.main)
integrated$singleR.labels <- pred$labels[match(rownames(integrated@meta.data), rownames(pred))]
DimPlot(integrated, reduction = "umap", group.by = "singleR.labels") + NoLegend()

# Plot SingleR results
plotScoreHeatmap(pred)
plotDeltaDistribution(pred)
tab <- table(Assigned = pred$labels, Clusters = integrated$seurat_clusters)
pheatmap(log10(tab + 10), color = colorRampPalette(c("darkblue", 'lavender'))(10))

# Set Idents as Seurat annotations
Idents(integrated) <- integrated@meta.data$singleR.labels
DimPlot(integrated, reduction = 'umap', label = TRUE)

# Step 1: Identify the first 10 unique labels
first_10_labels <- unique(integrated@meta.data$singleR.labels)[1:10]

# Step 2: Set identities based on the singleR.labels
Idents(integrated) <- integrated@meta.data$singleR.labels

# Step 3: Extract cells that match the first 10 labels
cells_to_keep <- WhichCells(integrated, idents = first_10_labels)
# Step 4: Subset the Seurat object with these cells
subset_integrated <- subset(integrated, cells = cells_to_keep)

# Step 5: Create the DimPlot
DimPlot(subset_integrated, reduction = 'umap', label = TRUE)

# Find markers between conditions
integrated$celltype.cnd <- paste0(integrated$singleR.labels, '_', integrated$orig.ident)
Idents(integrated) <- integrated$celltype.cnd
DimPlot(integrated, reduction = 'umap', label = TRUE)

all_Marker <- FindAllMarkers(integrated, logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE)

# Add a column indicating gene status
all_Marker$status <- ifelse(all_Marker$avg_log2FC >= 2 & all_Marker$p_val_adj < 0.0001, "Upregulated",
                            ifelse(all_Marker$avg_log2FC < 2 & all_Marker$p_val_adj < 0.0001, "Downregulated",
                                   "Not significant"))


# Write markers to Excel
write.xlsx(all_Marker, "singlecellDEGs_integrated-2.xlsx", rowNames = FALSE)
