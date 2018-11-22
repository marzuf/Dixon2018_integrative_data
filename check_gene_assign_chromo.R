
setDir <- "/media/electron"

ds1 <- "gene_data_final/GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460/genes2tad/all_genes_positions.txt"
ds2 <- "gene_data_final/GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP/genes2tad/all_genes_positions.txt"
ds3 <- "gene_data_final/GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP_vs_GSE105697_ENCFF364CWZ_T47D/genes2tad/all_genes_positions.txt"

ds4 <- file.path(setDir,
                "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt")


d1 <- read.delim(ds1, header=F, stringsAsFactors = F)
d1$V1 <- as.character(d1$V1)

d2 <- read.delim(ds2, header=F, stringsAsFactors = F)
d2$V1 <- as.character(d2$V1)

d3 <- read.delim(ds3, header=F, stringsAsFactors = F)
d3$V1 <- as.character(d3$V1)

d4 <- read.delim(ds4, header=F, stringsAsFactors = F)
d4$V1 <- as.character(d4$V1)

commonGenes <- Reduce(intersect, list(d1$V1, d2$V1, d3$V1, d4$V1))

d1 <- d1[d1$V1 %in% commonGenes,c("V1", "V2")]
d2 <- d2[d2$V1 %in% commonGenes,c("V1", "V2")]
d3 <- d3[d3$V1 %in% commonGenes,c("V1", "V2")]
d4 <- d4[d4$V1 %in% commonGenes,c("V1", "V2")]



all.equal(d1,d2)
all.equal(d1,d3)
all.equal(d1,d4)