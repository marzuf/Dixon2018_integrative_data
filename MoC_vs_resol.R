startTime <- Sys.time()

library(foreach)
library(doMC)
library(ggplot2)
library(stringi)

options(scipen=100)


cat("> START: MoC_vs_resol.R\n")
# Rscript MoC_vs_resol.R

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

source(file.path(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data", "utils_fct.R"))

if(SSHFS) setwd("/media/electron/mnt/etemp/marie/Dixon2018_integrative_data")

registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path("MOC_VS_RESOL")
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, "moc_vs_resol_logFile.txt")
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""

plotType <- "svg"

myHeight <- 7
myWidth <- 7

strWidthSplit <- 35

binSize <- 40000

consensusPattern <- ""

source("datasets_settings.R")

txt <- paste0("!! hard-coded bin size:\t", binSize, "\n")
printAndLog(txt, logFile)

### prepare resolution data
resolInFold <- file.path("CHECK_MATRESOL")
all_resol_files <- list.files(resolInFold, full.names=T, pattern="Rdata")
curr_file <- all_resol_files[1]
all_resol_DT <- foreach(curr_file = all_resol_files, .combine='rbind') %do% {
  eval(parse(text=load(curr_file)))
}
all_resol_DT$countSum_log10 <- log10(all_resol_DT$countSum)

all_resol_DT$ds <- unlist(sapply(all_resol_DT$dataset, function(x) names(ds_mapping[ds_mapping == x])))
stopifnot(!is.na(all_resol_DT$ds))

### prepare MoC data
mocInFold <- file.path("CMP_DATASETS_MOC")
mocFile <- file.path(mocInFold, "all_MoC_dt.Rdata")
stopifnot(file.exists(mocFile))
moc_DT <- eval(parse(text = load(mocFile)))

moc_DT <- moc_DT[moc_DT$ds1 %in% all_resol_DT$ds & moc_DT$ds2 %in% all_resol_DT$ds ,]

moc_DT$diffRowAbove1000 <- unlist(foreach(i = 1:nrow(moc_DT), .combine='c') %dopar% {
  ds1 <- moc_DT$ds1[i]
  ds2 <- moc_DT$ds2[i]
  curr_chromo <- moc_DT$chromo[i]
  stopifnot(ds1 %in% all_resol_DT$ds)
  stopifnot(ds2 %in% all_resol_DT$ds)
  abs(all_resol_DT$rowAbove1000[all_resol_DT$ds == ds1 & all_resol_DT$chromo == curr_chromo] - all_resol_DT$rowAbove1000[all_resol_DT$ds == ds2 & all_resol_DT$chromo == curr_chromo])
})


moc_DT$diffCountSumLog10 <- unlist(foreach(i = 1:nrow(moc_DT), .combine='c') %dopar% {
  ds1 <- moc_DT$ds1[i]
  ds2 <- moc_DT$ds2[i]
  curr_chromo <- moc_DT$chromo[i]
  stopifnot(ds1 %in% all_resol_DT$ds)
  stopifnot(ds2 %in% all_resol_DT$ds)
  abs(all_resol_DT$countSum_log10[all_resol_DT$ds == ds1 & all_resol_DT$chromo == curr_chromo] - all_resol_DT$countSum_log10[all_resol_DT$ds == ds2 & all_resol_DT$chromo == curr_chromo])
})


ref_var <- "MoC"

vars_to_plot <- c("diffRowAbove1000", "diffCountSumLog10")
plot_var <- vars_to_plot[1]

nDatasets <- length(unique(c(moc_DT$ds1, moc_DT$ds2)))

for(plot_var in vars_to_plot) {
  
  myx <- moc_DT[,plot_var]
  myy <- moc_DT[,ref_var]
  
  outFile <- file.path(outFold, paste0(ref_var, "_vs_", plot_var, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width = myWidth))
  plot(x = myx,
       y = myy,
       xlab = plot_var,
       ylab = ref_var,
       pch = 16, cex = 0.7,
       main = paste0(ref_var, " vs. ", plot_var))
  
  mtext(side=3, text = paste0("(nDatasets = ", nDatasets, ")"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}


######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




