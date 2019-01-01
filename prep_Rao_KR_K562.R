startTime <- Sys.time()

options(scipen=100)

### :
# Rscript prep_Rao_data.R <inFile> <chromo> <binSize> <newBinSize> <outFile>
# Rscript prep_Rao_data.R leukemia/K562/GSE63525/K562/10kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_10kb.RAWobserved chr1 10000 40000 leukemia/K562/GSE63525/GSE63525_K562_40kb_ICE_chr1_TopDom.matrix
# Rscript prep_Rao_data.R leukemia/K562/GSE63525/K562/10kb_resolution_intrachromosomal/chr21/MAPQGE30/chr21_10kb.RAWobserved chr21 10000 40000 leukemia/K562/GSE63525/GSE63525_K562_40kb_ICE_chr21_TopDom.matrix

# Rscript prep_Rao_data.R leukemia/K562/GSE63525/K562/10kb_resolution_intrachromosomal/chr9/MAPQGE30/chr9_10kb.RAWobserved chr9 10000 10000 leukemia/K562/GSE63525/GSE63525_K562_10kb_ICE_chr9_TopDom.matrix

# For example, here is a line from the GM12878_combined 5kb chr1 MAPQGE30 raw observed contact matrix 
# (GM12878_combined/5kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_5kb.RAWobserved):
#   40000000	40100000	59.0
# To normalize this entry using the KR normalization vector, one would divide 59.0 by the 8001st line ((40000000/5000)+1=8001) 
# and the 8021st line ((40100000/5000)+1=8021) of GM12878_combined/5kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_5kb.KRnorm.
# The 8001st line of the KR norm file is 1.2988778370674694;The 8021st line of the KR norm file is 1.6080499717941548. 
# So the corresponding KR normalized entry for the entry above is 59.0/(1.2988778370674694*1.6080499717941548) or 28.24776973966101.


args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 5)

inFile <- "leukemia/K562/GSE63525/K562/10kb_resolution_intrachromosomal/chr9/MAPQGE30/chr9_10kb.RAWobserved"
chromo <- "chr9"
binSize <- 10000
newBinSize <- 40000
outFile <- "foo.txt"

inFile <- "leukemia/K562/GSE63525/K562/10kb_resolution_intrachromosomal/chr10/MAPQGE30/chr10_10kb.RAWobserved"
chromo <- "chr10"
binSize <- 10000
newBinSize <- 40000
outFile <- "foo.txt"

inFile <- args[1]
chromo <- args[2]
binSize <- args[3]
newBinSize <- args[4]
outFile <- args[5]

biasFile <- gsub("RAWobserved$", "KRnorm", inFile)

binSize <- as.numeric(binSize)
newBinSize <- as.numeric(newBinSize)
stopifnot(!is.na(binSize))
stopifnot(!is.na(newBinSize))

stopifnot(file.exists(inFile))
stopifnot(file.exists(biasFile))

dir.create(dirname(outFile), recursive=TRUE)

SSHFS=F
setDir <- ifelse(SSHFS, "/media/electron", "")

suppressPackageStartupMessages(library(HiTC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(data.table, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # for preparing data for LGF
suppressPackageStartupMessages(library(Matrix, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # for preparing data for LGF
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # for preparing data for LGF
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # for preparing data for LGF

ifelse(SSHFS, registerDoMC(2,40))

source(file.path(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/ezh2_utils_fct.R"))
#source("createHTC.R")
##########################################################################################################################

cat(paste0("init bin size=", binSize, "\n"))
cat(paste0("new bin size=", newBinSize, "\n"))

### CHANGE HERE FOR GSE63525 !!!
binInCoordinates <- TRUE  
binfileHeader <- FALSE

ICE_maxiter <- 1000

countDT <- read.delim(inFile, stringsAsFactors = FALSE, header=FALSE, col.names=c("binA", "binB", "count"))
biasDT <- read.delim(biasFile, stringsAsFactors = FALSE, header=FALSE, col.names=c("bias"))

if(all(is.na(biasDT))) {
  
  stop("-- STOP: all biases are NA !!! \n")
  
}

#   40000000	40100000	59.0
# To normalize this entry using the KR normalization vector,
# one would divide 59.0 by the 8001st line ((40000000/5000)+1=8001) 
# and the 8021st line ((40100000/5000)+1=8021) of GM12878_combined/5kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_5kb.KRnorm.
# The 8001st line of the KR norm file is 1.2988778370674694;The 8021st line of the KR norm file is 1.6080499717941548. 
# So the corresponding KR normalized entry for the entry above is 59.0/(1.2988778370674694*1.6080499717941548) or 28.24776973966101.

newCount <- foreach(i = seq_len(nrow(countDT)), .combine='c') %dopar% {
  
  idx_line1 <- (countDT$binA[i]/binSize)+1
  idx_line2 <- (countDT$binB[i]/binSize)+1
  
  bias1 <- biasDT$bias[idx_line1]
  bias2 <- biasDT$bias[idx_line2]
  
  rawCount <- countDT$count[i]
  
  normCount <- rawCount(bias1*bias2)
  normCount
}
stopifnot(length(newCount) == nrow(countDT))
countDT$count <- newCount

######################## with HTC
htc_object <- createHTC(file=inFile, 
                        bin.size = binSize, 
                        chr=chromo, 
                        dim = -1, 
                        reindex = binInCoordinates, 
                        header=binfileHeader,
                        inputIsFile = FALSE
                        )
hicMat_v0 <- as.matrix(intdata(htc_object))
dim(hicMat_v0)
if(newBinSize > binSize) {
  htc_object <- binningC(htc_object, binsize=newBinSize, bin.adjust=FALSE, upa=TRUE, method="sum", optimize.by = "speed")
} 
hicMat_v0b <- as.matrix(intdata(htc_object))
rownames(hicMat_v0b) <- colnames(hicMat_v0b) <- NULL

######################## with sparseMatrix
countDT <- fread(inFile, header=F, stringsAsFactors = F, col.names = c("bin1", "bin2", "count"))
### CHANGE HERE FOR GSE63525 !!!
countDT$bin1 <- countDT$bin1/binSize
countDT$bin2 <- countDT$bin2/binSize
hicMat_v1 <- symMatrix_from_sparseListDT(countDT)
dim(hicMat_v1)

if(newBinSize > binSize) {
  rebinDT <- rebin_sparseMatrix(sparseCountDT=countDT, initBinSize = binSize, newBinSize=newBinSize, filled=FALSE)
  stopifnot(all(colnames(rebinDT) == c("bin1", "bin2", "count")))
  hicMat_v1b <- symMatrix_from_sparseListDT(rebinDT) 
} else {
  hicMat_v1b <- hicMat_v1
}

dim(hicMat_v1b)

######################## COMPARISON HERE:
stopifnot(all.equal(as.matrix(hicMat_v0), as.matrix(hicMat_v1), check.attributes=F))
stopifnot(all.equal(as.matrix(hicMat_v0b), as.matrix(hicMat_v1b), check.attributes=F))

all.equal(as.matrix(hicMat_v0), as.matrix(hicMat_v1), check.attributes=F)
all.equal(as.matrix(hicMat_v0b), as.matrix(hicMat_v1b), check.attributes=F)

stopifnot(all.equal(as.matrix(hicMat_v0), as.matrix(hicMat_v1), check.attributes=F))
stopifnot(all.equal(as.matrix(hicMat_v0b), as.matrix(hicMat_v1b), check.attributes=F))

######################################################################### ICE NORMALIZATION

cat("... start ICE normalization\n")

# htc_object_ICE <- normICE(htc_object, max_iter=ICE_maxiter)
htc_object_ICE <- htc_object

ice_matrix <- as.data.frame(as.matrix(intdata(htc_object_ICE)))

cat("... ICE matrix dim: ", paste0(dim(ice_matrix), collapse=" x "), "\n")

coordDT <- data.frame(chromo = chromo,
                      start = seq(from=0, by = newBinSize, length.out = nrow(ice_matrix)),
                      end = seq(from=newBinSize, by = newBinSize, length.out = nrow(ice_matrix)),
                      stringsAsFactors=F
)

outDT <- cbind(coordDT, ice_matrix)

cat("... write ICE normalized matrix to file\n")
write.table(outDT, file = outFile, col.names=F, row.names=F, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))

###
cat("*** done ***\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))