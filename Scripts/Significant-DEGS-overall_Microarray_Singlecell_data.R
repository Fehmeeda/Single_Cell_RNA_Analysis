library(readxl)
library(VennDiagram)
library(openxlsx)

# File paths for the single-cell RNA-Seq data (Identification of common DEGs of single cell Data)
file1 <- "C:/Users/PMLS/Documents/singlecellDEGs1.xlsx"              
file2 <- "C:/Users/PMLS/Documents/singlecellDEGsdataset2.xlsx"     

# Read data from Excel files
df1 <- read_excel(file1)
df2 <- read_excel(file2)

# Filter genes based on status (Upregulated and Downregulated)
upregulated_genes_df1 <- df1$Gene_Symbol[df1$status == "Upregulated"]
downregulated_genes_df1 <- df1$Gene_Symbol[df1$status == "Downregulated"]
upregulated_genes_df2 <- df2$Gene_Symbol[df2$status == "Upregulated"]
downregulated_genes_df2 <- df2$Gene_Symbol[df2$status == "Downregulated"]

# Find common upregulated and downregulated genes
common_upregulated_genes <- intersect(upregulated_genes_df1, upregulated_genes_df2)
common_downregulated_genes <- intersect(downregulated_genes_df1, downregulated_genes_df2)

# Create data frames for common genes with their status
common_upregulated_df <- data.frame(Gene_Symbol = common_upregulated_genes, Status = "Upregulated")
common_downregulated_df <- data.frame(Gene_Symbol = common_downregulated_genes, Status = "Downregulated")

# Write common upregulated genes and their status to an Excel file
write.xlsx(common_upregulated_df, "common_upregulated_genes.xlsx")

# Write common downregulated genes and their status to an Excel file
write.xlsx(common_downregulated_df, "common_downregulated_genes.xlsx")


# Create Venn diagram for common upregulated genes
venn_upregulated <- venn.diagram(
  x = list(Dataset1 = upregulated_genes_df1, Dataset2 = upregulated_genes_df2),
  category.names = c("file 1", "file 2"),
  main = "Venn Diagram of Common Upregulated Marker Genes across Single Cell RNA-Seq Data",
  filename = NULL,
  fill = c("lightgreen", "orange"), # Fill colors for the circles
  main.fontface = "bold" # Font style for the main title
  
)

# Create Venn diagram for common downregulated genes
venn_downregulated <- venn.diagram(
  x = list(Dataset1 = downregulated_genes_df1, Dataset2 = downregulated_genes_df2),
  category.names = c("file 1", " file 2"),
  main = "Venn Diagram of Common Downregulated Marker Genes across Single Cell RNA-Seq Data",
  filename = NULL,fill = c("skyblue", "pink"), # Fill colors for the circles
  main.fontface = "bold" # Font style for the main title
 
)

# Save Venn diagrams as TIFF
tiff("venn_upregulated-singlecell.tiff", width = 10, height = 10, units = "in", res = 300)
grid.draw(venn_upregulated)
dev.off()

tiff("venn_downregulated-singlecell.tiff", width = 10, height = 10, units = "in", res = 300)
grid.draw(venn_downregulated)
dev.off()

# Create a combined Excel file
combined_df <- rbind(common_upregulated_df, common_downregulated_df)
write.xlsx(combined_df, "combined_common_genes-singlecell.xlsx")

####################################

###################################

#### Identification of Common DEGs overall from microarray data and single cell RNA sequencing data ########

library(readxl)
library(VennDiagram)
library(grid)

# File paths for all Excel sheets
file1 <- "C:/Users/PMLS/Documents/GSE18520_DEGs.xlsx"
file2 <- "C:/Users/PMLS/Documents/GSE26712_DEGS.xlsx"
file3 <- "C:/Users/PMLS/Documents/GSE40595_DEGs.xlsx"
file4 <- "C:/Users/PMLS/Documents/GSE54388_DEG.xlsx"
file5 <- "C:/Users/PMLS/Documents/updated_file1-singlecell.xlsx"

file_paths <- c(file1, file2, file3, file4, file5)

# Lists to store gene sets for upregulated and downregulated genes from each sheet
upregulated_sets <- list()
downregulated_sets <- list()

# Read genes and their status from each Excel sheet
for (file_path in file_paths) {
  data <- na.omit(read_excel(file_path))
  upregulated_sets[[file_path]] <- data$Gene_Symbol[data$status == "Upregulated"]
  downregulated_sets[[file_path]] <- data$Gene_Symbol[data$status == "Downregulated"]
}

# Find common upregulated genes
common_upregulated <- Reduce(intersect, upregulated_sets)
# Find common downregulated genes
common_downregulated <- Reduce(intersect, downregulated_sets)

# Write common genes to text files
write.table(common_upregulated, "common_upregulated_genes-all.txt", row.names = FALSE, col.names = FALSE)
write.table(common_downregulated, "common_downregulated_genes-all.txt", row.names = FALSE, col.names = FALSE)

# Create a Venn diagram for upregulated genes
venn.plot.up <- venn.diagram(
  x = list(
    A = upregulated_sets[[1]],
    B = upregulated_sets[[2]],
    C = upregulated_sets[[3]],
    D = upregulated_sets[[4]],
    E = upregulated_sets[[5]]
  ),
  category.names = c("GSE18520", "GSE26712", "GSE40595", "GSE54388", "SingleCellData"),
  filename = NULL,
  fill = rainbow(5),
  main = list(label = "Venn Diagram of Common Upregulated DEGs across Microarray and Single cell RNA-Seq Data", cex = 1.5, fontface = "bold"),
  main.cex = 1.5,
  cex = 1.2,
  fontface = "bold",
  cat.fontface = "bold",
  cat.cex = 1.2,
  cat.default.pos = "outer",
  cat.dist = c(0.05, 0.05, 0.05, 0.05, 0.03),
  cat.pos = 0,
  scaled = TRUE,
  rotation.degree = 0,
  margin = 0.1
)

# Save the Venn diagram for upregulated genes in TIFF format
tiff("venn_diagram_upregulated_genes-all.tiff", width = 900, height = 800)
grid.newpage()
grid.draw(venn.plot.up)
dev.off()

# Create a Venn diagram for downregulated genes
venn.plot.down <- venn.diagram(
  x = list(
    A = downregulated_sets[[1]],
    B = downregulated_sets[[2]],
    C = downregulated_sets[[3]],
    D = downregulated_sets[[4]],
    E = downregulated_sets[[5]]
  ),
  category.names = c("GSE18520", " GSE26712", "GSE40595", "GSE54388", "SingleCellData"),
  filename = NULL,
  fill = rainbow(5),
  main = list(label = "Venn Diagram of Common Downregulated DEGs across microarray and Single Cell RNA-Seq Data", cex = 1.5, fontface = "bold"),
  main.cex = 1.5,
  cex = 1.2,
  fontface = "bold",
  cat.fontface = "bold",
  cat.cex = 1.2,
  cat.default.pos = "outer",
  cat.dist = c(0.05, 0.05, 0.05, 0.05, 0.03),
  cat.pos = 0,
  scaled = TRUE,
  rotation.degree = 0,
  margin = 0.1
)

# Save the Venn diagram for downregulated genes in TIFF format
tiff("venn_diagram_downregulated_genes-all.tiff", width = 900, height = 800)
grid.newpage()
grid.draw(venn.plot.down)
dev.off()
