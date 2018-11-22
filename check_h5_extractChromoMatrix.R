library(data.table)

options(scipen = 100)

startTime <- Sys.time()

printAndLog <- function(txt, logFile){
  cat(txt)
  cat(txt, file = logFile, append=T)
}

cat("> START check_h5_extractChromoMatrix.R\n")

# Rscript check_h5_extractChromoMatrix.R breast/T47D/ENCSR549MGQ_GSE105697/GSE105697_ENCFF364CWZ_chromatin_interactions_hg19_

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

# inFile <- file.path(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data/breast/T47D/ENCSR549MGQ_GSE105697/GSE105697_ENCFF364CWZ_chromatin_interactions_hg19_chr5_TopDom.matrix")
# stopifnot(file.exists(inFile))

matrixPattern <- file.path(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data/breast/T47D/ENCSR549MGQ_GSE105697/GSE105697_ENCFF364CWZ_chromatin_interactions_hg19_")

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1) {
  matrixPattern <- args[1]
}
stopifnot(file.exists(dirname(matrixPattern)))

outFold <- file.path("CHECK_H5_EXTRACTCHROMOMATRIX")
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, paste0(basename(matrixPattern), "_check_h5_extractchromomatrix_logFile.txt"))
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""

all_files <- list.files(dirname(matrixPattern), full.names = TRUE, pattern = basename(matrixPattern))

#all_files <- all_files[1]

for(inFile in all_files) {
  inMat <- fread(inFile)
  
  txt <- paste0("> START file: ", basename(inFile), "\n")
  printAndLog(txt, logFile)
  
  txt <- paste0("... matrix dimension:\t", nrow(inMat), "x", ncol(inMat) -3, "\n")
  printAndLog(txt, logFile)
  
  txt <- paste0("... found chromo:\t", unique(inMat$V1), "\n")
  printAndLog(txt, logFile)
  
  txt <- paste0("... found unique bin size:\t", unique(inMat$V3-inMat$V2)/1000, " kb\n")
  printAndLog(txt, logFile)
  
  start <- inMat$V2
  start <- start[2:length(start)]
  end <- inMat$V3
  end <- end[1:(length(end) - 1)]
  stopifnot(start == end)
  
  txt <- paste0("... bins are contiguous:\t", "OK", "\n")
  printAndLog(txt, logFile)
}

######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
