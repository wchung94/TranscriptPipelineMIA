gene_expression <- read.table("DiffParserOut.txt", header = TRUE)
Genes_FPKM <- read.table("genes.fpkm_tracking", header = TRUE)
row.names(Genes_FPKM) <- Genes_FPKM$gene_id
Useful_Genes_FPKM <- Genes_FPKM[Genes_FPKM$gene_id %in% gene_expression$gene_id, ]
Useful_FPKM_data <- Useful_Genes_FPKM[c(10,14,18)]
data.cordist <- as.dist(1-cor(Useful_FPKM_data))
hc.corout <- hclust(data.cordist, method = "complete")
plot(hc.corout, cex = 1)