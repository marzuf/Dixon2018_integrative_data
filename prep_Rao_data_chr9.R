startTime <- Sys.time()

options(scipen=100)


# !!! BEFORE normICE
# -> source("chr9_normICE.R") to avoid error at iteration 172 for chr9

### :
# Rscript prep_Rao_data.R <inFile> <chromo> <binSize> <newBinSize> <outFile>
# Rscript prep_Rao_data.R leukemia/K562/GSE63525/K562/10kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_10kb.RAWobserved chr1 10000 40000 leukemia/K562/GSE63525/GSE63525_K562_40kb_ICE_chr1_TopDom.matrix
# Rscript prep_Rao_data.R leukemia/K562/GSE63525/K562/10kb_resolution_intrachromosomal/chr21/MAPQGE30/chr21_10kb.RAWobserved chr21 10000 40000 leukemia/K562/GSE63525/GSE63525_K562_40kb_ICE_chr21_TopDom.matrix

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

binSize <- as.numeric(binSize)
newBinSize <- as.numeric(newBinSize)
stopifnot(!is.na(binSize))
stopifnot(!is.na(newBinSize))

stopifnot(file.exists(inFile))

system(paste0("mkdir -p ", dirname(outFile)))

SSHFS=F
setDir <- ifelse(SSHFS, "/media/electron", "")

suppressPackageStartupMessages(library(HiTC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(data.table, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # for preparing data for LGF
suppressPackageStartupMessages(library(Matrix, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # for preparing data for LGF

source(file.path(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/ezh2_utils_fct.R"))
#source("createHTC.R")
##########################################################################################################################

cat(paste0("init bin size=", binSize, "\n"))
cat(paste0("new bin size=", newBinSize, "\n"))

### CHANGE HERE FOR GSE63525 !!!
binInCoordinates <- TRUE  
binfileHeader <- FALSE

ICE_maxiter <- 1000

######################## with HTC
htc_object <- createHTC(file=inFile, 
                        bin.size = binSize, 
                        chr=chromo, 
                        dim = -1, 
                        reindex = binInCoordinates, 
                        header=binfileHeader)
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

# !!! avoid error at iteration 172
source("chr9_normICE.R")

save(htc_object,file="chr9_htc_object.Rdata")

htc_object_ICE <- normICE(htc_object, max_iter=ICE_maxiter)

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
