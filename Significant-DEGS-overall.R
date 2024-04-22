library(readxl)

#Reading single cell RNA-Seq data scripts
file1<-   "C:/Users/PMLS/Documents/singlecellDEGs1.xlsx"              
file2<-  "C:/Users/PMLS/Documents/singlecellDEGsdataset2.xlsx"     

df1 <- read_excel(file1)
df2 <- read_excel(file2)

genes_df1 <- df1$Gene_Symbol
genes_df2 <- df2$gene

common_genes <- intersect(genes_df1, genes_df2)

# Replace "common_genes.txt" with the desired path and file name
write.table(common_genes, "common_genes_single-cell.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


# Now we find Common DEGs between `microarray data scripts and Single cell data analysis scripts

file3<-   "C:/Users/PMLS/Documents/common_genes-DEGs.txt"              
file4<-  "C:/Users/PMLS/Documents/common_genes_single-cell.txt"  

# Read the first text file
genes_file1 <- readLines(file3)

#Remove double quotes:

genes_file1 <- gsub('"', '', genes_file1)

genes_file1
# Read the second text file
genes_file2 <- readLines(file4)

# Find common genes
common_genes <- intersect(genes_file1, genes_file2)

# Save common genes to a text file
writeLines(common_genes, "common_DEGs_both.txt")

`

