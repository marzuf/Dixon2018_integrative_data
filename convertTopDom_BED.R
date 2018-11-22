library(data.table)

options(scipen = 100)

startTime <- Sys.time()

printAndLog <- function(txt, logFile){
  cat(txt)
  cat(txt, file = logFile, append=T)
}

cat("> START convertTopDom_BED.R\n")
# Rscript convertTopDom_BED.R 40000 breast/MCF7/GSM1631185_GSE66733/MCF7/TopDom/GSM1631185_MCF7_40kb_chr.+.bed

SSHFS <- FALSE

setDir <- ifelse(SSHFS, "/media/electron", "")

# inFile <- file.path(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data/breast/MCF7/GSM1631185_GSE66733/MCF7/TopDom/GSM1631185_MCF7_40kb_chr1.bed")
# stopifnot(file.exists(inFile))

binSize <- 40000
matrixPattern <- file.path(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data/breast/MCF7/GSM1631185_GSE66733/MCF7/TopDom/GSM1631185_MCF7_40kb_chr.+.bed$")

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 2) {
  binSize <- as.numeric(args[1])
  matrixPattern <- args[2]
} else {
  stop("erorr")
}
stopifnot(!is.na(binSize))
stopifnot(file.exists(dirname(matrixPattern)))

outFold <- dirname(matrixPattern)

logFile <- file.path(outFold, gsub("\\.bed\\$", "_convert_topdom_logFile.txt", paste0(basename(matrixPattern))))
# logFile <- file.path(outFold, paste0(basename(matrixPattern), "_check_matresol_logFile.txt"))
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""

all_files <- list.files(dirname(matrixPattern), full.names = TRUE, pattern = basename(matrixPattern))

txt <- paste0("... set bin size:\t", binSize/1000, " kb", "\n")
printAndLog(txt, logFile)

stopifnot(length(all_files) <= 23)

for(inFile in all_files) {
  
  txt <- paste0("> START file: ", basename(inFile), "\n")
  printAndLog(txt, logFile)
  
  stopifnot(grepl("\\.bed$", inFile))
  outFile <- gsub("\\.bed$", "_final_domains.txt", inFile)
  
  tadDT <- read.delim(inFile, header=F, stringsAsFactors = FALSE)
  
  # !!! SO WEIRD !!!
  # IN FILE /mnt/etemp/marie/Dixon2018_integrative_data/breast/MCF7/GSM1631185_GSE66733/MCF7/TopDom/GSM1631185_MCF7_40kb_chr4.bed
  ## chr4    190680000       190680000       domain
  ##chr4    190680000       190960000       domain
  
  tadDT <- tadDT[tadDT[,3] > tadDT[,2],]
  
  coordDT <- tadDT[,c(1:3)]
  
  stopifnot(length(unique(coordDT[,1])) == 1)
  
  curr_chromo <- unique(coordDT[,1])
  
  start <- coordDT[,2]
  start <- start[2:length(start)]
  end <- coordDT[,3]
  end <- end[1:(length(end) - 1)]
  stopifnot(start == end)
  
  txt <- paste0("... found chromo:\t", curr_chromo,  "\n")
  printAndLog(txt, logFile)
  
  tadDT <- tadDT[ tadDT[,4] == "domain", ]
  stopifnot(nrow(tadDT) > 0)
  
  tadDT[,2] <- tadDT[,2] + 1
  
  outDT <- tadDT[,1:3]
  
  ### ROUND THE LAST DOMAIN !!!
  
  last_before <- outDT[nrow(outDT), 3]
  last_after <- ceiling(last_before/binSize)*binSize
  stopifnot(last_after >= last_before)
  outDT[nrow(outDT), 3] <- last_after
  
  stopifnot(outDT[,3] %% binSize == 0)
  stopifnot( (outDT[,2]-1) %% binSize == 0)
  
  txt <- paste0("... change last end from\t", last_before, "\tto\t", last_after, "\n" )
  printAndLog(txt, logFile)
  
  txt <- "\n"
  printAndLog(txt, logFile)
  
  write.table(outDT, file = outFile, sep="\t", quote=F, col.names=F, row.names=F)
  cat(paste0("... written: ", outFile, "\n"))
}
  
  
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




