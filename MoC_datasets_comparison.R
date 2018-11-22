library(foreach)
library(ggplot2)

outFolder <- "MOC_DATASETS_COMPARISON_CONSENSUS"
system(paste0("mkdir -p ", outFolder))

plotType <- "svg"
myHeightGG <- 7 
myWidthGG <- 10

#MOC_BETWEEN_DATASETS/GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460/chr_MoC_table.txt
# setwd("/media/electron/mnt/etemp/marie/Dixon2018_integrative_data")

# chrOrder <- paste0("chr", c(1:22, "X"))

mainFold <- "MOC_BETWEEN_DATASETS"
mainFold <- "MOC_BETWEEN_DATASETS_CONSENSUS"

all_files <- list.files(mainFold, full.names = TRUE, pattern="chr_MoC_table.txt$", recursive = TRUE)

moc_file <- "MOC_BETWEEN_DATASETS/GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460/chr_MoC_table.txt"
all_MoC_dt <- foreach(moc_file = all_files, .combine='rbind') %do% {
  dt <- read.delim(moc_file, header=T, stringsAsFactors = F)
  dt$comparison <- basename(dirname(moc_file))  
  dt
}

all_MoC_dt$comparison_label <- gsub("_vs_", " vs.\n", all_MoC_dt$comparison)


lungCells <- c("A549", "NCIH460")
breastCells <- c("T47D", "MCF7") 

all_MoC_dt$comparison_label <- gsub("_vs_", " vs.\n", all_MoC_dt$comparison)

all_MoC_dt$cell_type <- sapply(all_MoC_dt$comparison, function(x) {
  ifelse(any(sapply(lungCells, function(i) grepl(i,x))), "lung", 
         ifelse(any(sapply(breastCells, function(i) grepl(i,x))), "breast", NA))})

if(any(is.na(all_MoC_dt$cell_type))) {

all_MoC_dt$cell_type <- ifelse(grepl("LUNG", all_MoC_dt$comparison), "lung", "breast")

}

stopifnot(!is.na(all_MoC_dt$cell_type))



boxCol <- c(lung = "dodgerblue4", breast = "darkorange2")
  

p_all <- ggplot(all_MoC_dt, aes(x = comparison_label, y = MoC, fill = cell_type)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter() +
  scale_x_discrete(name="")+
  scale_y_continuous(name=paste0("MoC values"),
                     breaks = scales::pretty_breaks(n = 5))+
  scale_fill_manual(values = boxCol)+
  labs(fill="")+
  ggtitle(label = paste0("MoC between datasets"))+
  theme( # Increase size of axis lines
    # top, right, bottom and left
    # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size=10),
    panel.grid = element_blank(),
    # panel.grid.major = element_line(colour = "lightpink"),
    # strip.text.x = element_text(),
    axis.text.x = element_text( hjust=1,vjust = 0.5, angle=90, size=10),
    axis.line.x = element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .3, color = "black"),
    #    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
    axis.title.y = element_text(color="black", size=12),
    axis.title.x = element_text(color="black", size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.key = element_blank()
    # axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1)
  ) #+
# geom_hline(yintercept = 1, linetype = 2)

outFile <- file.path(outFolder, paste0("MoC_comparisons_boxplot.", plotType))
ggsave(p_all, filename = outFile, height = myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


meanMoC_dt <- aggregate(MoC ~ comparison, FUN=mean, data = all_MoC_dt)
meanMoC_dt$MoC <- round(meanMoC_dt$MoC, 4)
meanMoC_dt <- meanMoC_dt[order(meanMoC_dt$MoC, decreasing=TRUE),]
colnames(meanMoC_dt)[colnames(meanMoC_dt) == "MoC"] <- "meanMoC"
outFile <- file.path(outFolder, paste0("meanMoC_dt.txt"))
write.table(meanMoC_dt, file = outFile, sep="\t", quote=F, col.names=T, row.names=F)
cat(paste0("... written: ", outFile, "\n"))
write.table(meanMoC_dt, file = "", sep="\t", quote=F, col.names=T, row.names=F)


