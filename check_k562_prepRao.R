library(data.table)
setwd("/media/electron/mnt/etemp/marie/Dixon2018_integrative_data")

raoFile <- "leukemia/K562_hic/GSE63525_K562_chr9_10kb_hic_observedNONE.txt"
raoDT <- fread(raoFile)
colnames(raoDT) <- c("binA", "binB", "count")
head(raoDT)

prepFile <- "leukemia/K562_hic/GSE63525_K562_chr9_10kb_hic_observedNONE_ICE_noNormTopDom.matrix"
prepDT <- fread(prepFile)
prepDT[1:5,1:5]