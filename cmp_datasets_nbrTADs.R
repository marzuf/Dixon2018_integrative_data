startTime <- Sys.time()

library(foreach)
library(doMC)
library(ggplot2)
library(stringi)

options(scipen=100)

# done differently previously in: /mnt/ed4/marie/TAD_call_pipeline_DI/compare_TopDom_tissue_consensus

cat("> START: cmp_datasets_nbrTADs.R\n")
# Rscript cmp_datasets_nbrTADs.R

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

source(file.path(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data", "utils_fct.R"))

if(SSHFS) setwd("/media/electron/mnt/etemp/marie/Dixon2018_integrative_data")

registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path("CMP_DATASETS_NBRTADS")
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "cmp_datasets_nbrTADs_logFile.txt")
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""

plotType <- "svg"

widthBoxplot <- 10
heightBoxplot <- 7

strWidthSplit <- 35

binSize <- 40000

txt <- paste0("!! hard-coded bin size:\t", binSize, "\n")
printAndLog(txt, logFile)

consensusPattern <- "_conservedTADs.txt$"

source("datasets_settings.R")

i=1
ds1="breastCL1"
all_nbr_dt <- foreach(ds1 = all_ds, .combine="rbind") %dopar% {
  
  all_files1 <- eval(parse(text = paste0(ds1, "Files")))
  
  name1 <- eval(parse(text = paste0(ds1, "name")))
  
  cat(paste0("*** START: ", ds1,  "\n"))
    
  chromo = "chr1"
  chr_nbr_dt <- foreach(chromo = intersectChromos, .combine='rbind') %do% {
    
    cat(paste0("> ", ds1, "  - ", chromo, "\n"))
    
    file1 <- all_files1[grepl(paste0(chromo, "_"), basename(all_files1)) & grepl(paste0(name1), all_files1)]
    stopifnot(length(file1) == 1)
    
    
    dt1 <- read.delim(file1, stringsAsFactors = FALSE, header=F, col.names = c("chromo", "start", "end"))
    if(is.character(dt1[1,2]) & is.character(dt1[2,2])) {
      dt1 <- read.delim(file1, stringsAsFactors = FALSE, header=T, col.names = c("chromo", "start", "end"))
    }
    
    
    stopifnot(ncol(dt1) == 3)
    stopifnot(is.numeric(dt1[,2]))
    stopifnot(is.numeric(dt1[,3]))
    
    head(dt1, 2)
    
    chrEnd <- dt1[nrow(dt1), 3]
    chromoCover <- sum(dt1[,3] - dt1[,2] + 1)/chrEnd
    stopifnot(chromoCover >= 0 & chromoCover <= 1)
    
    data.frame(
      ds1 = ds1,
      chromo = chromo,
      chromoCover = chromoCover,
      nTADs = nrow(dt1),
      meanSizeTADs = mean(dt1$end - dt1$start + 1),
      stringsAsFactors = FALSE
    )
  }
  chr_nbr_dt
}

outFile <- file.path(outFold, "all_nbr_dt.Rdata")
save(all_nbr_dt, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

all_nbr_dt <- eval(parse(text = load(outFile)))
stopifnot(nrow(all_nbr_dt) == (length(intersectChromos) * length(all_ds)))

all_nbr_dt$chromo <- factor(as.character(all_nbr_dt$chromo), levels = paste0("chr", c(1:22, "X")))

stopifnot(!is.na(all_nbr_dt))

all_vars <- colnames(all_nbr_dt)[!colnames(all_nbr_dt) %in% c("ds1", "chromo")]

plot_tit <- c(
meanSizeTADs = "TAD mean size (bp)",
chromoCover = "Chromosome coverage",
nTADs = "Nbr of TADs"
)

stopifnot(all_vars %in% names(plot_tit))


for(var_to_plot in all_vars) {
  
  mytit <- plot_tit[var_to_plot]
  
  mean_all_nbr_dt <- aggregate(as.formula(paste0(var_to_plot, " ~ ds1")), FUN=mean, data = all_nbr_dt)
  mean_all_nbr_dt <- mean_all_nbr_dt[order(mean_all_nbr_dt[, var_to_plot], decreasing = TRUE),]
  all_nbr_dt$ds1 <- factor(as.character(all_nbr_dt$ds1), levels = mean_all_nbr_dt$ds1)
  
  p_common <- ggplot(all_nbr_dt, aes_string(x = "ds1", y = var_to_plot)) + 
    geom_boxplot(outlier.shape=NA) +
    # geom_jitter(aes(colour = chromo)) +
    scale_x_discrete(name="")+
    # scale_y_continuous(name=paste0("-log10(", padjVarGO, ")"),
#    scale_y_continuous(name=paste0(var_to_plot),
    scale_y_continuous(name=paste0(mytit),
                       breaks = scales::pretty_breaks(n = 10))+ #, limits = c(0, max(auc_DT_m$value)+0.05))+
    # coord_cartesian(expand = FALSE) +
    # scale_fill_manual(values = c(selectGenes = "dodgerblue4", selectTADs_genes = "darkorange2"),
    #                   labels = c(selectGenes = "selectGenes", selectTADs_genes = "selectTADs_genes"))+
    # scale_colour_manual(values = c(selectGenes = "dodgerblue4", selectTADs_genes = "darkorange2"),
    #                     labels = c(selectGenes = "selectGenes", selectTADs_genes = "selectTADs_genes"), guide = F)+
    labs(colour  = "") +
#    ggtitle(label = paste0(var_to_plot))+
    ggtitle(label = paste0(mytit))+
    theme( # Increase size of axis lines
      # top, right, bottom and left
      # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size=10),
      panel.grid = element_blank(),
      # panel.grid.major = element_line(colour = "lightpink"),
      # strip.text.x = element_text(),
      axis.text.x = element_text( hjust=1,vjust = 0.5, size=12, angle = 90),
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






