options(scipen = 100)

library(doMC)
library(foreach)

startTime <- Sys.time()

cat("> START find_consensusTADs.R\n")
# Rscript find_consensusTADs.R  \
# breast/MCF7/GSM1631185_GSE66733/MCF7/TopDom GSM1631185_MCF7_40kb_chr.+final_domains.txt GSM1631185_MCF7 \
# breast/MCF7/GSM1942100_GSM1942101_GSE75070/shGFP/TopDom GSE75070_MCF7_shGFP_40kb_chr.+final_domains.txt GSE75070_MCF7_shGFP \
# breast/T47D/ENCSR549MGQ_GSE105697/TopDom GSE105697_ENCFF364CWZ_T47D_40kb_chr.+final_domains.txt GSE105697_ENCFF364CWZ_T47D


SSHFS <- FALSE

registerDoMC(ifelse(SSHFS, 2, 40))

setDir <- ifelse(SSHFS, "/media/electron", "")

if(SSHFS) setwd("/media/electron/mnt/etemp/marie/Dixon2018_integrative_data")

source("utils_fct.R")

source("wsm_TADconsensus_withCoverage_version2_fct_withAllTxt.R")


#*****************************
inFolder1 <- file.path(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data/breast/MCF7/GSM1631185_GSE66733/MCF7/TopDom")
inPattern1 <- "GSM1631185_MCF7_40kb_chr.+final_domains.txt"
name1 <- "GSM1631185_MCF7"
inFolder2 <- file.path(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data/breast/MCF7/GSM1942100_GSM1942101_GSE75070/shGFP/TopDom")
inPattern2 <- "GSE75070_MCF7_shGFP_40kb_chr.+final_domains.txt"
name2 <- "GSE75070_MCF7_shGFP"
inFolder3 <- file.path(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data/breast/T47D/ENCSR549MGQ_GSE105697/TopDom")
inPattern3 <- "GSE105697_ENCFF364CWZ_T47D_40kb_chr.+final_domains.txt"
name3 <- "GSE105697_ENCFF364CWZ_T47D"
args <- c(inFolder1, inPattern1, name1,
          inFolder2, inPattern2, name2,
          inFolder3, inPattern3, name3)

#*****************************
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %% 3 == 0)

all_folders <- args[seq(from=1, to=length(args), by=3)]
all_patterns <- args[seq(from=2, to=length(args), by=3)]
all_names <- args[seq(from=3, to=length(args), by=3)]

stopifnot(file.exists(all_folders))

outFolder <- file.path("FIND_CONSENSUS_TADS", paste0(all_names, collapse="_vs_"))
system(paste0("mkdir -p ", outFolder))

logFile <- file.path(outFolder, paste0(paste0(all_names, collapse="_vs_"), "_consensusTADs.logFile.txt"))
if(!SSHFS) system(paste0("rm -f ", logFile))
if(SSHFS) logFile <- ""




# hard-coded settings for the consensus
txt  <- paste0("*** HARD-CODED SETTINGS FOR TAD CALLING ALGORITHM:\n")
printAndLog(txt, logFile)

set_tolRad <- 80000
set_conservThresh <- ifelse(length(args)/3 == 3, 2/3, 
                            ifelse(length(args)/3 == 2, 1, NA))
if(is.na(set_conservThresh)) stop("error: check set_conservThresh\n")
set_weightValue <- NULL
set_coverageThresh <-  1
set_tadfileHeader <-  FALSE
set_chrSize <- NULL
set_ncpu <- ifelse(SSHFS, 2, 40)
set_fileAsInput <- TRUE
set_logFile <- NULL

txt  <- paste0("... set_tolRad\t=\t", set_tolRad, "\n")
printAndLog(txt, logFile)
txt  <- paste0("... set_conservThresh\t=\t", set_conservThresh, "\n")
printAndLog(txt, logFile)
txt  <- paste0("... set_weightValue\t=\t", set_weightValue, "\n")
printAndLog(txt, logFile)
txt  <- paste0("... set_coverageThresh\t=\t", set_coverageThresh, "\n")
printAndLog(txt, logFile)
txt  <- paste0("... set_tadfileHeader\t=\t", as.character(set_tadfileHeader), "\n")
printAndLog(txt, logFile)
txt  <- paste0("... set_chrSize\t=\t", as.character(set_chrSize), "\n")
printAndLog(txt, logFile)
txt  <- paste0("... set_ncpu\t=\t", set_ncpu, "\n")
printAndLog(txt, logFile)
txt  <- paste0("... set_fileAsInput\t=\t", set_fileAsInput, "\n")
printAndLog(txt, logFile)
txt  <- paste0("... set_logFile\t=\t", as.character(set_logFile), "\n")
printAndLog(txt, logFile)

txt  <- paste0("*** SETTINGS RETRIEVED FROM COMMAND LINE:\n")
printAndLog(txt, logFile)

txt  <- paste0("... all_folders\t=\t", paste0(all_folders, collapse=", "), "\n")
printAndLog(txt, logFile)
txt  <- paste0("... all_patterns\t=\t", paste0(all_patterns, collapse=", "), "\n")
printAndLog(txt, logFile)
txt  <- paste0("... all_names\t=\t", paste0(all_names, collapse=", "), "\n")
printAndLog(txt, logFile)


###############################################################################################################
###############################################################################################################
###############################################################################################################

all_files <- lapply(seq_along(all_folders), function(i) {
  list.files(all_folders[i], full.names = TRUE, pattern=all_patterns[i])
})

stopifnot( unlist(lapply(all_files, length)) <= 23)

all_chromos <- lapply(seq_along(all_files), function(i) {
  chrs <- gsub(".*(chr[0-9]+|chrX).+", "\\1", basename(all_files[[i]]))
  stopifnot(chrs %in% paste0("chr", c(1:22, "X")))
  chrs
})

intersectChromo <- Reduce(intersect, all_chromos)

txt <- paste0("... found common chromo:\t", length(intersectChromo), "\n")
printAndLog(txt, logFile)

# intersectChromo <- intersectChromo[1]
# intersectChromo <- "chrX"

for(chromo in intersectChromo) {
  
  cat(paste0("> START chromo = ", chromo, "\n"))
  
  chromo_files <- lapply(seq_along(all_files), function(i) {
    curr_files <- all_files[[i]]
    curr_chr_file <- curr_files[grep(paste0(chromo, "_"), basename(curr_files))]
    stopifnot(file.exists(curr_chr_file))
    curr_chr_file
  })
  
    domainsDT <- get_consensus_TADs(chromo_files,
                                    tolRad=set_tolRad,
                                    conservThresh=set_conservThresh,
                                    weightValue = set_weightValue,
                                    coverageThresh = set_coverageThresh,
                                  tadfileHeader = set_tadfileHeader,
                                  chrSize=set_chrSize,
                                  ncpu=set_ncpu,
                                  fileAsInput=set_fileAsInput,
                                  logFile=set_logFile)
  
    
  cat("nrow consensus:", nrow(domainsDT), "\n")
  
  if(nrow(domainsDT) > 1 ){
    for(i in 2:nrow(domainsDT)) {
      # cat(paste0(i, "\t"))
      stopifnot( domainsDT[i, 2] > domainsDT[i-1, 3]  )
    }
  }
  outFile <- file.path(outFolder, paste0(chromo, "_conservedTADs.txt"))
  write.table(domainsDT, file = outFile, sep="\t", quote=F, col.names=F, row.names=F)
  cat(paste0("... written:\t", outFile, "\n"))  
}
 

######################################################################################
cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




