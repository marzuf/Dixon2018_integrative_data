library(data.table)

options(scipen = 100)

startTime <- Sys.time()

printAndLog <- function(txt, logFile){
  cat(txt)
  cat(txt, file = logFile, append=T)
}

cat("> START check_matResol.R\n")
# Rscript check_matResol.R breast/T47D/ENCSR549MGQ_GSE105697/GSE105697_ENCFF364CWZ_chromatin_interactions_hg19_

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


outFold <- file.path("CHECK_MATRESOL")
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, paste0(basename(matrixPattern), "_check_matresol_logFile.txt"))
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""

all_files <- list.files(dirname(matrixPattern), full.names = TRUE, pattern = basename(matrixPattern))

#all_files <- all_files[1]

for(inFile in all_files) {
  
  txt <- paste0("> START file: ", basename(inFile), "\n")
  printAndLog(txt, logFile)
  
  inMat <- fread(inFile)
  
  coordDT <- inMat[,c(1:3)]
  
  stopifnot(length(unique(coordDT[,1])) == 1)
  stopifnot(length(unique(coordDT[,3]-coordDT[,2])) == 1)
  
  binSize <- unique(coordDT[,3]-coordDT[,2])
  curr_chromo <- unique(coordDT[,1])
  
  txt <- paste0("... found chromo:\t", curr_chromo, "\n")
  printAndLog(txt, logFile)
  
  txt <- paste0("... found bin size:\t", binSize/1000, " kb", "\n")
  printAndLog(txt, logFile)
  
  hic_DT <- inMat[,-c(1:3)]
  stopifnot(nrow(hic_DT) == ncol(hic_DT))
  
  
  txt <- paste0("... ", curr_chromo, " - matrix dim.:\t", paste0(dim(hic_DT), collapse = " x "), "\n")
  printAndLog(txt, logFile)
  
  matrixRowSum <- rowSums(hic_DT, na.rm=T)
  
  txt <- "... summary matrix row sum:\n"
  printAndLog(txt, logFile)
  sink(logFile, append=T)
  print(summary(matrixRowSum))
  sink()
  
  txt <- paste0("... # rows with >= 1000 counts:\t", sum(matrixRowSum >= 1000), "/", length(matrixRowSum), " (", round(sum(matrixRowSum >= 1000)/length(matrixRowSum)*100, 2),"%)\n")
  printAndLog(txt, logFile)
  txt <- "\n"
  printAndLog(txt, logFile)
}
  
  
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




