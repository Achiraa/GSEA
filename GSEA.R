setwd("C:/geo_data/GFP")

library(enrichplot)
library(ggplot2)
library("org.Mm.eg.db")
library(clusterProfiler)
library(gridExtra)

file1=read.csv("forpathwayanalysisGFP.csv")
head(file1)

# we want the log2 fold change 
ensembl_gene_list1 = file1$log2FoldChange
ensembl_gene_list1

# name the vector
names(ensembl_gene_list1) <- file1$geneids

# sort the list in decreasing order (required for clusterProfiler)
gene_list1 = sort(ensembl_gene_list1, decreasing = TRUE)
View(gene_list1)

resultgseaBP= gseGO(geneList= gene_list1, 
           ont ="BP", 
           keyType = "ENSEMBL",
           minGSSize = 3, 
           maxGSSize = 800, 
           pvalueCutoff = 0.05, 
           verbose = TRUE, 
           OrgDb = org.Mm.eg.db, 
           pAdjustMethod = "none")

resultgseaCC= gseGO(geneList= gene_list1, 
                   ont ="CC", 
                   keyType = "ENSEMBL",
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   verbose = TRUE, 
                   OrgDb = org.Mm.eg.db, 
                   pAdjustMethod = "none")
resultgseaMF= gseGO(geneList= gene_list1, 
                    ont ="MF", 
                    keyType = "ENSEMBL",
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = TRUE, 
                    OrgDb = org.Mm.eg.db, 
                    pAdjustMethod = "none")

head(resultgseaBP)
head(resultgseaCC)
head(resultgseaMF)

BP1=dotplot(resultgseaBP, showCategory=20) + ggtitle("Dotplot for GSEA:Gene Ontology:BP")
CC1=dotplot(resultgseaCC, showCategory=20) + ggtitle("Dotplot for GSEA:Gene Ontology:CC")
MF1=dotplot(resultgseaMF, showCategory=20) + ggtitle("Dotplot for GSEA:Gene Ontology:MF")
  
BP1
CC1
MF1

grid.arrange(BP1, CC1, MF1, nrow =1)

#BP=barplot(resultgseaBP, showCategory=4)+ ggtitle("Dotplot for GSEA:Gene Ontology:Biological Processes")
#CC=barplot(resultgseaCC, showCategory=7)+ ggtitle("Dotplot for GSEA:Gene Ontology:Cellular Components")
#MF=barplot(resultgseaMF, showCategory=10)+ ggtitle("Dotplot for GSEA:Gene Ontology:Molecular Functions")



