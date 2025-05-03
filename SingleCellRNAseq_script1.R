# Load libraries
library(Seurat)
library(celldex)
library(SingleR)
library(openxlsx)
library(pheatmap)

# Set Seurat object assay version
options(Seurat.object.assay.version = "v3")

#sample1
sample1 <- Read10X_h5(filename = "C:/Users/Muzi/Downloads/GSE118127_RAW/GSM3557963_sample13_B1_i12C_filtered_gene_bc_matrices_h5.h5")
# Initialize the Seurat object with the raw (non-normalized data).
sample1 <- CreateSeuratObject(counts = sample1, project = "Normal", min.cells = 3, min.features = 200)

#sample2
sample2 <- Read10X_h5(filename = "C:/Users/Muzi/Downloads/GSE118127_RAW/GSM3319039_sample_1-8_filtered_gene_bc_matrices_h5.h5")
# Initialize the Seurat object with the raw (non-normalized data).
sample2 <- CreateSeuratObject(counts = sample2, project = "Normal", min.cells = 3, min.features = 200)

#sample3
sample3 <- Read10X_h5(filename = "C:/Users/Muzi/Downloads/GSE118127_RAW/GSM3319047_sample_3-6_filtered_gene_bc_matrices_h5.h5")
# Initialize the Seurat object with the raw (non-normalized data).
sample3 <- CreateSeuratObject(counts = sample3, project = "Normal", min.cells = 3, min.features = 200)

# Now load ovarian samples and also create seurat objects:

#sample1
sample4<- Read10X_h5(filename = "C:/Users/Muzi/Downloads/GSE233615_RAW/GSM7431438_E7_filtered_feature_bc_matrix.h5")
# Initialize the Seurat object with the raw (non-normalized data).
sample4 <- CreateSeuratObject(counts = sample4, project = "Ovarian", min.cells = 3, min.features = 200)

#sample2
sample5<- Read10X_h5(filename = "C:/Users/Muzi/Downloads/GSE233615_RAW/GSM7431439_E9_filtered_feature_bc_matrix.h5")
# Initialize the Seurat object with the raw (non-normalized data).
sample5 <- CreateSeuratObject(counts = sample5, project = "Ovarian", min.cells = 3, min.features = 200)

#sample3
sample6 <- Read10X_h5(filename = "C:/Users/Muzi/Downloads/GSE233615_RAW/GSM7431440_E10_filtered_feature_bc_matrix.h5")
# Initialize the Seurat object with the raw (non-normalized data).
sample6 <- CreateSeuratObject(counts = sample6, project = "Ovarian", min.cells = 3, min.features = 200)

#sample4
sample7 <- Read10X_h5(filename = "C:/Users/Muzi/Downloads/GSE233615_RAW/GSM7431441_E11_filtered_feature_bc_matrix.h5")
# Initialize the Seurat object with the raw (non-normalized data).
sample7 <- CreateSeuratObject(counts = sample7, project = "Ovarian", min.cells = 3, min.features = 200)

# Create a vector of cell IDs
cell_ids <- paste0("sample", 1:7)  
# Merge Seurat Objects Normal
merged_samples<- merge(sample1, y = c(sample2, sample3,sample4,sample5,sample6,sample7),          
                       add.cell.ids = cell_ids,
                       project = 'merged_samplesNormalandOvarian')
#view(merged_samples@meta.data)

# Perform quality control
merged_samples <- PercentageFeatureSet(merged_samples, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(merged_samples, features = c("nFeature_RNA", "nCount_RNA" , "percent.mt") ,  ncol = 3 )

merged_samples <- subset(merged_samples, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 25000 & percent.mt <= 25)
VlnPlot(merged_samples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Normalize data and find variable features
merged_samples <- NormalizeData(merged_samples, normalization.method = "LogNormalize", scale.factor = 10000)
merged_samples <- FindVariableFeatures(merged_samples, selection.method = "vst", nfeatures = 2000)

# Top features and featurePlot

variable_genes<-merged_samples@assays[["RNA"]]@var.features
write.table(variable_genes, file = "variable_genes.txt", sep = "\t", row.names = FALSE)

#VariableFeaturePlot(merged_samples)
#Identify the 10 most highly variable genes:

top10 <- head(VariableFeatures(merged_samples), 10)
plot1 <- VariableFeaturePlot(merged_samples)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

# Data scaling
merged_samples <- ScaleData(merged_samples) #2000 identified variable features 

all.genes <- rownames(merged_samples)
merged_samples <- ScaleData(merged_samples, features = all.genes) 

merged_samples@assays[["RNA"]]@data@x

#View(merged_samples@commands)

#Perform PCA on the scaled data (linear dimensional reduction)

merged_samples <- RunPCA(merged_samples, features = VariableFeatures(object = merged_samples))

# Examine and visualize PCA results a few different ways
# DimPlot(), VizDimReduction() and DimHeatmap ()
DimPlot(merged_samples, reduction = "pca", dims = c(1,2))
DimPlot(merged_samples, reduction = "pca", dims = c(1, 10))
DimPlot(merged_samples, reduction = "pca", dims = c(1, 50))

# Determine the ‘dimensional’ of the dataset
merged_samples <- JackStraw(merged_samples, num.replicate = 100)
merged_samples <- ScoreJackStraw(merged_samples, dims = 1:20)
JackStrawPlot(merged_samples, dims = 1:20)

ElbowPlot(merged_samples)
ElbowPlot(merged_samples, ndims = 50, reduction = "pca")

# Run non-linear dimensional reduction (UMAP/tSNE)
merged_samples <- RunUMAP(merged_samples, dims = 1:20)
DimPlot(merged_samples, reduction = "umap", label = TRUE, repel = TRUE)

merged_samples <- RunTSNE(object = merged_samples)
DimPlot(object = merged_samples, reduction = "tsne")

merged_samples <- RunTSNE(object = merged_samples)
DimPlot(object = merged_samples, reduction = "tsne")
DimPlot(object = merged_samples, reduction = "tsne", group.by = 'orig.ident')+labs(title = "", x = "tSNE_1", y = "tSNE_1", color = "Sample Type")

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
DimPlot(integrated, reduction = "pca", dims = c(1, 2))
DimPlot(integrated, reduction = "pca", dims = c(1, 10))
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

#Visualization
DimPlot(integrated,reduction="umap",label = TRUE)

# Generate the UMAP plots
plot1 <- DimPlot(merged_samples, reduction = "umap", group.by = 'orig.ident')+labs(title = "Before\n Batch Effect Correction", x = "UMAP 1", y = "UMAP 2", color = "Sample Type")

plot2 <- DimPlot(integrated, reduction = "umap", group.by = 'orig.ident') +
  labs(title = "After\n Batch Effect Correction", x = "UMAP 1", y = "UMAP 2", color = "Sample Type")

# Combine the plots
combined_plot <- plot1 + plot2

# Display the combined plot
print(combined_plot)

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
pheatmap(log10(tab + 10), color = colorRampPalette(c('white', 'blue'))(10))

# Set Idents as Seurat annotations
Idents(integrated) <- integrated@meta.data$singleR.labels
DimPlot(integrated, reduction = 'umap', label = TRUE)

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
write.xlsx(all_Marker, "singlecellDEGs_integrated-1.xlsx", rowNames = FALSE)