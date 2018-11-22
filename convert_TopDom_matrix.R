startTime <- Sys.time()
cat(paste0("> Rscript convert_TopDom_matrix.R\n"))

# Rscript convert_TopDom_matrix.R <inFile> <outFile>
# Rscript convert_TopDom_matrix.R HiCStein-MCF7-WT__hg19__chr1__C-40000-iced.matrix HiCStein-MCF7-WT__hg19__chr1__C-40000-iced_TopDom.matrix


options(scipen=100)


inFile = "/mnt/etemp/marie/Dixon2018_integrative_data/breast/MCF7/GSM1631185/Hi-C_MCF7_MCF10A_processed_HiCfiles/Heatmaps/chrxchr/40kb/MCF7/HiCStein-MCF7-WT__hg19__chr1__C-40000-iced.matrix"
outFile = "/mnt/etemp/marie/Dixon2018_integrative_data/breast/MCF7/GSM1631185/Hi-C_MCF7_MCF10A_processed_HiCfiles/Heatmaps/chrxchr/40kb/MCF7/HiCStein-MCF7-WT__hg19__chr1__C-40000-iced_TopDom.matrix"

library(data.table)

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 2)

inFile <- args[1]
outFile <- args[2]

inMat <- fread(inFile)
stopifnot(ncol(inMat) == nrow(inMat) + 1)

retrievedDim <- as.numeric(as.character(gsub("(.+)x(.+)", "\\1", colnames(inMat)[1])))
stopifnot(!is.na(retrievedDim))

coord <- inMat[,1]
head(coord)
coord_chromo <- gsub(".+\\|.+\\|(.+):(.+)-(.+)", "\\1", unlist(coord))
stopifnot(length(unique(coord_chromo)) == 1)
stopifnot(grepl("chr", coord_chromo))
coord_start <- as.numeric(as.character( gsub(".+\\|.+\\|(.+):(.+)-(.+)", "\\2", unlist(coord))))
stopifnot(!is.na(coord_start))
# for TopDom: -1
stopifnot(coord_start[1] == 1)
coord_start <- coord_start - 1
coord_end <- as.numeric(as.character( gsub(".+\\|.+\\|(.+):(.+)-(.+)", "\\3", unlist(coord))))
head(coord_end)
stopifnot(!is.na(coord_end))

inMat[,1] <- NULL
stopifnot(ncol(inMat) == nrow(inMat))
stopifnot(ncol(inMat) == retrievedDim)

# replace NA

cat("!!! Warning: replace NA and NaN with 0 !!!\n")

cat(paste0("... # of NA: ", sum(is.na(inMat)),  "/", retrievedDim*retrievedDim, " = ",
           round( sum(is.na(inMat))/(retrievedDim*retrievedDim) * 100, 2), " %\n"))
inMat[is.na(inMat)] <- 0

coordDT <- data.frame(chromo = coord_chromo,
                      start = coord_start,
                      end = coord_end,
                      stringsAsFactors =FALSE)

topdomDT <- cbind(coordDT, inMat)

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

