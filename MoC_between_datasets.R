options(scipen = 100)

startTime <- Sys.time()

printAndLog <- function(txt, logFile){
  cat(txt)
  cat(txt, file = logFile, append=T)
}

cat("> START MoC_between_datasets.R\n")
# Rscript MoC_between_datasets.R  40000 \
# breast/MCF7/GSM1631185_GSE66733/MCF7/TopDom GSM1631185_MCF7_40kb_chr.+final_domains.txt GSM1631185_MCF7 \
# /mnt/ed4/marie/TAD_call_pipeline_TopDom/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final chr.+_conservedTADs.txt PIPELINE_CONSENSUS

SSHFS <- FALSE

setDir <- ifelse(SSHFS, "/media/electron", "")

source("utils_fct.R")

#*****************************
binSize <- 40000

inFolder1 <- file.path(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data/breast/MCF7/GSM1631185_GSE66733/MCF7/TopDom")
inPattern1 <- "GSM1631185_MCF7_40kb_chr.+final_domains.txt"
name1 <- "GSM1631185_MCF7"
  
inFolder2 <- file.path(setDir, "/mnt/ed4/marie/TAD_call_pipeline_TopDom/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final")
inPattern2 <- "chr.+_conservedTADs.txt"
name2 <- "PIPELINE_CONSENSUS"
#*****************************

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 7)

binSize <- as.numeric(as.character(args[1]))

inFolder1 <- args[2]
inPattern1 <- args[3]
name1 <- args[4]

inFolder2 <- args[5]
inPattern2 <- args[6]
name2 <- args[7]

cat(paste0("... binSize\t=\t", binSize, "\n"))
cat(paste0("... inFolder1\t=\t", inFolder1, "\n"))
cat(paste0("... inPattern1\t=\t", inPattern1, "\n"))
cat(paste0("... name1\t=\t", name1, "\n"))
cat(paste0("... inFolder2\t=\t", inFolder2, "\n"))
cat(paste0("... inPattern2\t=\t", inPattern2, "\n"))
cat(paste0("... name2\t=\t", name2, "\n"))

###############################################################################################################
###############################################################################################################
###############################################################################################################

stopifnot(file.exists(inFolder1))
stopifnot(file.exists(inFolder2))
stopifnot(!is.na(binSize))

outFolder <- file.path("MOC_BETWEEN_DATASETS", paste0(name1, "_vs_", name2))
system(paste0("mkdir -p ", outFolder))

logFile <- file.path(outFolder, paste0(name1, "_vs_", name2, "_MoC_between_datasets.logFile.txt"))
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""

inFiles1 <- list.files(inFolder1, full.names = TRUE, pattern=inPattern1)
inFiles2 <- list.files(inFolder2, full.names = TRUE, pattern=inPattern2)

stopifnot(length(inFiles1) <= 23)
stopifnot(length(inFiles2) <= 23)

all_chromo1 <- gsub(".*(chr[0-9]+|chrX).+", "\\1", basename(inFiles1))
stopifnot(all_chromo1 %in% paste0("chr", c(1:22, "X")))

all_chromo2 <- gsub(".*(chr[0-9]+|chrX).+", "\\1", basename(inFiles2))
stopifnot(all_chromo2 %in% paste0("chr", c(1:22, "X")))

intersectChromo <- intersect(all_chromo1, all_chromo2)

txt <- paste0("... found common chromo:\t", length(intersectChromo), "\n")
printAndLog(txt, logFile)

# intersectChromo <- intersectChromo[1]

all_mocs <- c()

for(chromo in intersectChromo) {
  
  
  txt <- paste0("> START chromo:\t", chromo, "\n")
  printAndLog(txt, logFile)
  
  file1 <- inFiles1[grep(paste0(chromo, "_"), basename(inFiles1))]
  stopifnot(length(file1) == 1)
  stopifnot(file.exists(file1))
  
  file2 <- inFiles2[grep(paste0(chromo, "_"), basename(inFiles2))]
  stopifnot(length(file2) == 1)
  stopifnot(file.exists(file2))

  txt <- paste0("... file1=", file1, "\n")
  printAndLog(txt, logFile)
  txt <- paste0("... file2=", file2, "\n")
  printAndLog(txt, logFile)
    
  dt1 <- read.delim(file1, stringsAsFactors = F, header=F, col.names=c("chromo", "start", "end"))
  dt2 <- read.delim(file2, stringsAsFactors = F, header=F, col.names=c("chromo", "start", "end"))
  
  if(dt1$chromo[1] == "chromo") {
    dt1 <- read.delim(file1, stringsAsFactors = F, header=T, col.names=c("chromo", "start", "end"))
  }
  if(dt2$chromo[1] == "chromo") {
    dt2 <- read.delim(file2, stringsAsFactors = F, header=T, col.names=c("chromo", "start", "end"))
  }
  
  head(dt1,3)
  head(dt2,3)

  stopifnot(is.numeric(dt1$start))
  stopifnot(is.numeric(dt1$end))
  stopifnot(is.numeric(dt2$start))
  stopifnot(is.numeric(dt2$end))
  
  last_before1 <- dt1$end[nrow(dt1)]
  last_after1 <- ceiling(last_before1/binSize)*binSize
  stopifnot(last_after1 >= last_before1)
  dt1$end[nrow(dt1)] <- last_after1
  stopifnot(dt1$end %% binSize == 0)
  stopifnot( (dt1$start-1) %% binSize == 0)
  txt <- paste0("... change last end from\t", last_before1, "\tto\t", last_after1, "\n" )
  printAndLog(txt, logFile)
  
  last_before2 <- dt2$end[nrow(dt2)]
  last_after2 <- ceiling(last_before2/binSize)*binSize
  stopifnot(last_after2 >= last_before2)
  dt2$end[nrow(dt2)] <- last_after2
  stopifnot(dt2$end %% binSize == 0)
  stopifnot( (dt2$start-1) %% binSize == 0)
  txt <- paste0("... change last end from\t", last_before2, "\tto\t", last_after2, "\n" )
  printAndLog(txt, logFile)
  
  minEnd <- min(c(dt1$end[nrow(dt1)], dt2$end[nrow(dt2)]))
  
  txt <- paste0("... nrow before cropping to same minEnd\t", "dt1=", nrow(dt1), "\t", "dt2=", nrow(dt2), "\n")
  printAndLog(txt, logFile)
  
  dt1 <- dt1[dt1$end <= minEnd,]
  dt2 <- dt2[dt2$end <= minEnd,]
  
  txt <- paste0("... nrow after cropping to same minEnd\t", "dt1=", nrow(dt1), "\t", "dt2=", nrow(dt2), "\n")
  printAndLog(txt, logFile)
  
  moc <- calculate_MoC_with_domainTypes(set1=dt1, set2=dt2, chr_len=minEnd,
                                        gaps_as_clusters = TRUE, file_as_input = FALSE)
  
  txt <- paste0("*** MoC\t=\t", round(moc, 4), "\n")
  printAndLog(txt, logFile)
  
  txt <- paste0("\n")
  printAndLog(txt, logFile)
  
  all_mocs <- c(all_mocs, moc)
  
}


mocDT <- data.frame(chromo = intersectChromo,
                    MoC = round(all_mocs,4),
                    stringsAsFactors = FALSE)

write.table(mocDT, col.names=T, row.names=F, sep="\t", quote=F, file="")

write.table(mocDT, col.names=T, row.names=F, sep="\t", quote=F, file=logFile, append=T)

outFile <- file.path(outFolder, "chr_MoC_table.txt")
write.table(mocDT, col.names=T, row.names=F, sep="\t", quote=F, file=outFile)
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




