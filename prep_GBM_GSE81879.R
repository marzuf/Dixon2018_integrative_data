library(data.table)
library(foreach)
library(doMC)

options(scipen=100)

source("utils_fct.R")

startTime <- Sys.time()

SSHFS <- T
setDir <- ifelse(SSHFS, "/media/electron", "")

registerDoMC(ifelse(SSHFS, 2, 40))

inFile <- file.path(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data/glioblastoma/GSE81879/GSM2176966_GB176_norm_int.csv")
# <outPrefix>_chr#_<outSuffix>
outPrefix <- file.path(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data/glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176966_GB176")
outSuffix <- paste0("_norm_int_count.txt")
binSize <- 40000

# Rscript prep_GBM_GSE81879.R /mnt/etemp/marie/Dixon2018_integrative_data/glioblastoma/GSE81879/GSM2176966_GB176_norm_int.csv glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176966_GB176 norm_int_count.txt 40000

# Rscript prep_GBM_GSE81879.R <inFile> <outPrefix> <outSuffix> <binSize>

inFile <- "/mnt/etemp/marie/Dixon2018_integrative_data/glioblastoma/GSE81879/GSM2176966_GB176_norm_int.csv"
outPrefix <- "glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176966_GB176" 
outSuffix <- "norm_int_count.txt" 
binSize <- 40000

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)
inFile <- args[1]
outPrefix <- args[2]
outSuffix <- args[3]
binSize <- args[4]
binSize <- as.numeric(binSize)
stopifnot(!is.na(binSize))

outFold <- dirname(outPrefix)
system(paste0("mkdir -p ", outFold))

logFile <- file.path(outFold, paste0("prep_GBM_", basename(outPrefix), "_logFile.txt"))
if(!SSHFS) system(paste0("rm -f ", logFile))

cat(paste0("... read input file\n"))

interactionsDT <- read.delim(inFile, stringsAsFactors = FALSE, sep = ",", check.names = FALSE)
dim(interactionsDT)

cat(paste0("... melt dataframe\n"))
stopifnot(colnames(interactionsDT)[1] == "")
countDT <- melt(interactionsDT, id="")
head(countDT)
colnames(countDT)[colnames(countDT) == ""] <- "binA_frag"
colnames(countDT)[colnames(countDT) == "variable"] <- "binB_frag"
colnames(countDT)[colnames(countDT) == "value"] <- "count"

head(countDT)

# remove the 0 count
stopifnot(is.numeric(countDT$count))

txt <- paste0("... remove zero counts, # rows before = ", nrow(countDT), "\n")
printAndLog(txt, logFile)
countDT <- countDT[countDT$count > 0,]
txt <- paste0("... remove zero counts, # rows before = ", nrow(countDT), "\n")
printAndLog(txt, logFile)

head(countDT)

countDT$chrA <- gsub("(chr.+):.+", "\\1", countDT$binA_frag)
countDT$chrB <- gsub("(chr.+):.+", "\\1", countDT$binB_frag)

countDT$chrA <- as.character(countDT$chrA)
countDT$chrB <- as.character(countDT$chrB)

countDT$binA_fragStart <-gsub("chr.+:(.+)-.+", "\\1", countDT$binA_frag)
countDT$binA_fragEnd <-gsub("chr.+:.+-(.+)", "\\1", countDT$binA_frag)

countDT$binB_fragStart <-gsub("chr.+:(.+)-.+", "\\1", countDT$binB_frag)
countDT$binB_fragEnd <-gsub("chr.+:.+-(.+)", "\\1", countDT$binB_frag)

countDT$binA_fragStart <- as.numeric(as.character(countDT$binA_fragStart))
countDT$binA_fragEnd <- as.numeric(as.character(countDT$binA_fragEnd))
countDT$binB_fragStart <- as.numeric(as.character(countDT$binB_fragStart))
countDT$binB_fragEnd <- as.numeric(as.character(countDT$binB_fragEnd))

head(countDT)

stopifnot(!is.na(countDT$binA_fragStart))
stopifnot(!is.na(countDT$binA_fragEnd))
stopifnot(!is.na(countDT$binB_fragStart))
stopifnot(!is.na(countDT$binB_fragEnd))

# remove already the interchromo counts for speed:

txt <- paste0("... remove interchromo counts, # rows before = ", nrow(countDT), "\n")
printAndLog(txt, logFile)

countDT <- countDT[countDT$chrA == countDT$chrB,]
stopifnot(countDT$chrA == countDT$chrB)

txt <- paste0("... remove interchromo counts, # rows after = ", nrow(countDT), "\n")
printAndLog(txt, logFile)

head(countDT)

countDT$chromo <- countDT$chrA
countDT$chrA <- NULL
countDT$chrB <- NULL

# assign the end and start of the fragments to the bins
# multiply by binSize to have it in the same format as Rao data
countDT$binA_start <- (countDT$binA_fragStart %/% binSize) * binSize
countDT$binA_end <- (countDT$binA_fragEnd %/% binSize)  * binSize
countDT$binB_start <- (countDT$binB_fragStart %/% binSize) * binSize
countDT$binB_end <- (countDT$binB_fragEnd %/% binSize) * binSize


rangeBinA <- range(countDT$binA_fragEnd - countDT$binA_fragStart)
txt <- paste0("rangeBinA = ", paste0(rangeBinA, collapse=" - "), "\n")
printAndLog(txt, logFile)
rangeBinB <- range(countDT$binB_fragEnd - countDT$binB_fragStart)
txt <- paste0("rangeBinB = ", paste0(rangeBinB, collapse=" - "), "\n")
printAndLog(txt, logFile)

head(countDT)

# remove if the fragment span many bins:

#                txt <- paste0("... remove fragments spanning many bins, # rows before = ", nrow(countDT), "\n")
#                printAndLog(txt, logFile)

#                countDT <- countDT[(countDT$binA_start == countDT$binA_end) &
#                                     (countDT$binB_start == countDT$binB_end),]

#                txt <- paste0("... remove fragments spanning many bins, # rows after = ", nrow(countDT), "\n")
#                printAndLog(txt, logFile)

#                stopifnot(countDT$binA_start == countDT$binA_end)
#                stopifnot(countDT$binB_start == countDT$binB_end)

                #countDT$binA <- countDT$binA_start
                #countDT$binB <- countDT$binB_start
                #countDT$binA_start <- NULL
                #countDT$binA_end <- NULL
                #countDT$binB_start <- NULL
                #countDT$binB_end <- NULL

# I loose all ??? !!!


#                txt <- paste0("... remove fragments spanning many bins, # rows before = ", nrow(countDT), "\n")
#                printAndLog(txt, logFile)

#                countDT <- countDT[(countDT$binA_start == countDT$binA_end) & 
#                                     (countDT$binB_start == countDT$binB_end),]

#                txt <- paste0("... remove fragments spanning many bins, # rows after = ", nrow(countDT), "\n")
#                printAndLog(txt, logFile)

#                stopifnot(countDT$binA_start == countDT$binA_end)
#                stopifnot(countDT$binB_start == countDT$binB_end)

                countDT$binA <- countDT$binA_start
                countDT$binB <- countDT$binB_start
                countDT$binA_start <- NULL
                countDT$binA_end <- NULL
                countDT$binB_start <- NULL
                countDT$binB_end <- NULL


stopifnot(is.numeric(countDT$count))

# iterate over the chromosomes and aggregate the counts

all_chromo <- unique(countDT$chromo)
stopifnot(length(all_chromo) > 0)

cat("... found following chromo: ", paste0(all_chromo, collapse=","), "\n")

#all_chromo <- all_chromo[1]

foo <- foreach(chromo = all_chromo) %dopar% {
  
  cat("... start ", chromo, "\n")
  
  chromo_countDT <- countDT[countDT$chromo == chromo,]
  
  dt <- as.data.table(chromo_countDT[,c("binA", "binB", "count")])
  agg_dt <- dt[,  list(count=sum(count)), by=c('binA', 'binB')]

  agg_dt$binA2 <- pmin(agg_dt$binA, agg_dt$binB)
  agg_dt$binB2 <- pmax(agg_dt$binA, agg_dt$binB)

  agg_dt$binA <- agg_dt$binA2
  agg_dt$binB <- agg_dt$binB2
  
  agg_dt$binA2 <- NULL
  agg_dt$binB2 <- NULL

  outFile <- paste0(outPrefix, "_", chromo, "_", binSize/1000, "kb_", outSuffix )
  write.table(agg_dt, file = outFile, row.names = FALSE, col.names=F, sep="\t", quote=F)
  cat(paste0("... written: ", outFile, "\n"))
}



######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



