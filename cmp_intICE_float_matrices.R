SSHFS <- TRUE
setDir <- ifelse(SSHFS, "~/media/electron", "")

intMatrixFile <- file.path(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF714TMN/GSE105318_ENCFF714TMN_chromatin_interactions_hg19_chr1_TopDom.matrix")# → int ???
iceMatrixFile <- file.path(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF714TMN/GSE105318_ENCFF714TMN_chromatin_interactions_hg19_chr1_ICE_TopDom.matrix")
floatMatrixFile <- file.path(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF439QFU/GSE105318_ENCFF439QFU_chromatin_interactions_hg19_chr1_TopDom.matrix")# → float"

intMatrixDT <- read.delim(intMatrixFile, header=F, stringsAsFactors = FALSE)
dim(intMatrixDT)
intMatrixDT[1:5,1:5]
intMatrix_s <- intMatrixDT
intMatrixDT <- intMatrixDT[,-c(1:3)]
dim(intMatrixDT)
intLT <- intMatrixDT[lower.tri(intMatrixDT)]

iceMatrixDT <- read.delim(iceMatrixFile, header=F, stringsAsFactors = FALSE)
dim(iceMatrixDT)
iceMatrixDT[1:5,1:5]
iceMatrix_s <- iceMatrixDT
iceMatrixDT <- iceMatrixDT[,-c(1:3)]
dim(iceMatrixDT)
iceLT <- iceMatrixDT[lower.tri(iceMatrixDT)]


floatMatrixDT <- read.delim(floatMatrixFile, header=F, stringsAsFactors = FALSE)
dim(floatMatrixDT)
floatMatrixDT[1:5,1:5]
floatMatrix_s <- floatMatrixDT
floatMatrixDT <- floatMatrixDT[,-c(1:3)]
dim(floatMatrixDT)
floatLT <- floatMatrixDT[lower.tri(floatMatrixDT)]


stopifnot(length(floatLT) == length(iceLT))
stopifnot(length(intLT) == length(iceLT))

head(floatLT);head(intLT);head(iceLT)

cor.test(floatLT, intLT)
0.8789686 

cor.test(floatLT, iceLT)
0.9842409 

cor.test(iceLT, intLT)
0.8722768 
