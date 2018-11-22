startTime <- Sys.time()
cat(paste0("> Rscript convert_h5_TopDom_matrix.R\n"))

# Rscript convert_h5_TopDom_matrix.R <inFile> <chr> 
# Rscript convert_h5_TopDom_matrix.R GSE105600_ENCFF852YOE_chromatin_interactions_hg19.h5 chr1


options(scipen=100)

inFile = "/mnt/etemp/marie/Dixon2018_integrative_data/lung/A549/ENCSR444WCZ/GSE105600_ENCFF852YOE_chromatin_interactions_hg19.h5"
chromo = "chr1"

library(data.table)
library(rhdf5)

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 2)

inFile <- args[1]
chromo <- args[2]

outFile <- gsub("\\.h5$", paste0("_", chromo, "_TopDom.matrix"), inFile)

cat("... list data:\n")

dataList <- h5ls(inFile)
stopifnot("chrs" %in% dataList$name)

dataList

cat("... load chrs\n")

all_chromos <- h5read(inFile, "chrs")
stopifnot(chromo %in% all_chromos)

stopifnot("bin_positions" %in% dataList$name)

cat("... load bin_positions\n")
hicmat_binPos <- h5read(inFile, "bin_positions")

cat("... load interaction matrix\n")
stopifnot("interactions" %in% dataList$name)
hicMat <- h5read(inFile, "interactions")

chrIdx <- which(all_chromos == chromo)
stopifnot(length(chrIdx) == 1)

stopifnot(ncol(hicmat_binPos) == ncol(hicMat))

matToKeep <- which(hicmat_binPos[1,] == chrIdx)
length(matToKeep)

hicMat_chromo <- as.data.frame(hicMat[matToKeep, matToKeep])

gc()
gc()

# stopifnot(nrow(hicMat_chromo) == ncol(hicMat_chromo))
# stopifnot(nrow(hicMat_chromo) == length(matToKeep))

cat(paste0("... chromo ", chromo, " - matrix dim:\t", nrow(hicMat_chromo), "x", ncol(hicMat_chromo), "\n"))

bin_start <- hicmat_binPos[2,matToKeep] - 1
stopifnot(bin_start %% 10 == 0)
stopifnot(length(unique(diff(bin_start))) == 1)

bin_end <- hicmat_binPos[3,matToKeep]
stopifnot(length(bin_end) == length(bin_start))
stopifnot(ncol(hicMat_chromo) == length(bin_start))

# round the last bin if needed
stopifnot(length(unique(diff(bin_end))) == 1 |length(unique(diff(bin_end))) ==  2)
if(length(unique(diff(bin_end))) == 2){
  bin_end[length(bin_end)] <-   bin_end[length(bin_end)-1] + unique(diff(bin_start))
}
stopifnot(length(unique(diff(bin_end))) == 1)
stopifnot( unique(diff(bin_end)) == unique(diff(bin_start)))
stopifnot(bin_end %% 10 == 0)
stopifnot(bin_end > bin_start)

cat(paste0("... bin size:\t", unique(diff(bin_start)), "\n"))

coord_DT <- data.frame(
  chromo = chromo,
  start = bin_start,
  end = bin_end,
stringsAsFactors = FALSE
)

stopifnot(nrow(coord_DT) == nrow(hicMat_chromo))

retrievedDim <- nrow(coord_DT)

inMat <- hicMat_chromo

# replace NA

cat("!!! Warning: replace NA and NaN with 0 !!!\n")

cat(paste0("... # of NA: ", sum(is.na(inMat)),  "/", retrievedDim*retrievedDim, " = ",
           round( sum(is.na(inMat))/(retrievedDim*retrievedDim) * 100, 2), " %\n"))
inMat[is.na(inMat)] <- 0

topdomDT <- cbind(coord_DT, inMat)

stopifnot(nrow(topdomDT) + 3 == ncol(topdomDT))

cat("... write matrix to file ...\n")
write.table(topdomDT, file = outFile, sep="\t", quote=F, row.names = F, col.names = F)
stopifnot(file.exists(outFile))
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
