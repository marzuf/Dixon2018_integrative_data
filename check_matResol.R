library(data.table)
library(foreach)

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

dataset <- gsub("(.+)_chromatin_inter.+", "\\1", basename(matrixPattern))

chromoLevels <- paste0("chr", c(1:22, "X"))

check_resolDT <- foreach(inFile = all_files, .combine='rbind') %do% {
  
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
  
  outFile <- file.path(outFold, paste0(dataset, "_", chromo, "_matrixRowSum.Rdata"))
  save(matrixRowSum, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  countSum <- sum(matrixRowSum)
  
  txt <- paste0("... matrix count sum:\t", round(countSum, 4), "\n")
  printAndLog(txt, logFile)
  
  txt <- "... summary matrix row sum:\n"
  printAndLog(txt, logFile)
  sink(logFile, append=T)
  print(summary(matrixRowSum))
  sink()
  
  rowAbove1000 <- sum(matrixRowSum >= 1000)/length(matrixRowSum)
  
  txt <- paste0("... # rows with >= 1000 counts:\t", sum(matrixRowSum >= 1000), "/", length(matrixRowSum), " (", round(rowAbove1000*100, 2),"%)\n")
  printAndLog(txt, logFile)
  txt <- "\n"
  printAndLog(txt, logFile)
  
  tmpDT <- data.frame(
    dataset = dataset,
    chromo = curr_chromo,
    countSum = countSum,
    rowAbove1000 = rowAbove1000,
    stringsAsFactors = FALSE
  )
  colnames(tmpDT)[2] <- "chromo"
  
  tmpDT
}

check_resolDT$countSum <- round(check_resolDT$countSum, 4)
check_resolDT$rowAbove1000 <- round(check_resolDT$rowAbove1000, 4)

stopifnot(check_resolDT$chromo %in% chromoLevels)
check_resolDT$chromo <- factor(check_resolDT$chromo, levels = chromoLevels)
check_resolDT <- check_resolDT[order(as.numeric(check_resolDT$chromo), decreasing=F),]

outFile <- file.path(outFold, paste0(dataset, "_check_resolDT.Rdata"))
save(check_resolDT, file=outFile)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, paste0(dataset, "_check_resolDT.txt"))
write.table(check_resolDT, col.names=T, row.names=F, sep="\t", quote=F, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




