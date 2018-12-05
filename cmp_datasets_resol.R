startTime <- Sys.time()

library(foreach)
library(doMC)
library(ggplot2)
library(stringi)

options(scipen=100)


cat("> START: cmp_datasets_resol.R\n")
# Rscript cmp_datasets_resol.R

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

source(file.path(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data", "utils_fct.R"))

if(SSHFS) setwd("/media/electron/mnt/etemp/marie/Dixon2018_integrative_data")

registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path("CMP_DATASETS_RESOL")
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "cmp_datasets_resol_logFile.txt")
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""

plotType <- "svg"

widthBoxplot <- 10
heightBoxplot <- 7

strWidthSplit <- 35

binSize <- 40000

txt <- paste0("!! hard-coded bin size:\t", binSize, "\n")
printAndLog(txt, logFile)


mainFold <- "CHECK_MATRESOL"

consensusPattern=""
source("datasets_settings.R")

all_files <- list.files(mainFold, full.names=T, pattern="Rdata")

curr_file <- all_files[1]
all_resol_DT <- foreach(curr_file = all_files, .combine='rbind') %do% {
  eval(parse(text=load(curr_file)))
}

all_resol_DT$countSum_log10 <- log10(all_resol_DT$countSum)

all_vars <- colnames(all_resol_DT)[!colnames(all_resol_DT) %in% c("dataset", "chromo")]

#all_resol_DT$datasetLabel <- unlist(sapply(all_resol_DT$dataset, function(x) paste0(stri_wrap(str = x, width = strWidthSplit), collapse="\n")))

all_resol_DT$datasetLabel <- unlist(sapply(all_resol_DT$dataset, function(x) names(ds_mapping[ds_mapping == x])))

var_to_plot = "rowAbove1000"
for(var_to_plot in all_vars) {
  
  mean_all_resol_DT <- aggregate(as.formula(paste0(var_to_plot, " ~ dataset")), FUN=mean, data = all_resol_DT)
  mean_all_resol_DT <- mean_all_resol_DT[order(mean_all_resol_DT[, var_to_plot], decreasing = TRUE),]
  all_resol_DT$dataset <- factor(as.character(all_resol_DT$dataset), levels = mean_all_resol_DT$dataset)
  
  mean_all_resol_DT <- aggregate(as.formula(paste0(var_to_plot, " ~ datasetLabel")), FUN=mean, data = all_resol_DT)
  mean_all_resol_DT <- mean_all_resol_DT[order(mean_all_resol_DT[, var_to_plot], decreasing = TRUE),]
  all_resol_DT$datasetLabel <- factor(as.character(all_resol_DT$datasetLabel), levels = mean_all_resol_DT$datasetLabel)
  
  # p_common <- ggplot(all_resol_DT, aes_string(x = "dataset", y = var_to_plot)) + 
  p_common <- ggplot(all_resol_DT, aes_string(x = "datasetLabel", y = var_to_plot)) + 
    geom_boxplot(outlier.shape=NA) +
    # geom_jitter(aes(colour = chromo)) +
    scale_x_discrete(name="")+
    # scale_y_continuous(name=paste0("-log10(", padjVarGO, ")"),
    scale_y_continuous(name=paste0(var_to_plot),
                       breaks = scales::pretty_breaks(n = 10))+ #, limits = c(0, max(auc_DT_m$value)+0.05))+
    # coord_cartesian(expand = FALSE) +
    # scale_fill_manual(values = c(selectGenes = "dodgerblue4", selectTADs_genes = "darkorange2"),
    #                   labels = c(selectGenes = "selectGenes", selectTADs_genes = "selectTADs_genes"))+
    # scale_colour_manual(values = c(selectGenes = "dodgerblue4", selectTADs_genes = "darkorange2"),
    #                     labels = c(selectGenes = "selectGenes", selectTADs_genes = "selectTADs_genes"), guide = F)+
    labs(colour  = "") +
    ggtitle(label = paste0(var_to_plot))+
    theme( # Increase size of axis lines
      # top, right, bottom and left
      # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size=10),
      panel.grid = element_blank(),
      # panel.grid.major = element_line(colour = "lightpink"),
      # strip.text.x = element_text(),
      axis.text.x = element_text( hjust=1,vjust = 0.5, size=8, angle = 90),
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
  
  if(SSHFS) p_all
  
  p_dot <- p_common +  geom_jitter(aes(colour = chromo))
  if(SSHFS) p_dot
  
  p_txt <- p_common + geom_text(aes(label=chromo, colour=chromo, fontface="bold"),size=2.5, position = position_jitter(w = 0.3)) + guides(colour = "none")
  if(SSHFS) p_txt
  
  
  outFile <- file.path(outFold, paste0(var_to_plot, "_boxplot_chromoDots.", plotType))
  ggsave(plot=p_dot, file = outFile, width = widthBoxplot, height = heightBoxplot)
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFold, paste0(var_to_plot, "_boxplot_chromoLabs.", plotType))
  ggsave(plot=p_txt, file = outFile, width = widthBoxplot, height = heightBoxplot)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
}



######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))






