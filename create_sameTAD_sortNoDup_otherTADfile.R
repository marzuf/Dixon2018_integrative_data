startTime <- Sys.time()
cat(paste0("> Rscript create_sameTAD_sortNoDup_otherTADfile.R\n"))

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

# Rscript create_sameTAD_sortNoDup_otherTADfile.R GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal
# or
# Rscript create_sameTAD_sortNoDup_otherTADfile.R gene_data_final/GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal/genes2tad/all_genes_positions.txt CREATE_SAME_TAD_SORTNODUP/GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal


### UPDATE sortNoDup 30.06.2018
# -> sort rows of the gene coord. table to ensure alphabetical order of the genes !!
#    so that after melt the 1st gene will always be the 1st in alphabetical order
# -> add as.character() in apply !!! use "gene1" etc. instead of 1 index in apply 

### => this version 27.12.2018: pass the input TAD list and output folder 
### (used to prepare data, based on TAD file other than the one from the pipeline)

# outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "CREATE_SAME_TAD_SORTNODUP")
# gene2tadDT_file <-  paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
# gene2tadDT_file <- file.path("gene_data_final/GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal/genes2tad/all_genes_positions.txt")
# outFold <- file.path("CREATE_SAME_TAD_SORTNODUP", "GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal")

caller <- "TopDom"

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 1) {
  ds <- args[1]
  outFold <- file.path("CREATE_SAME_TAD_SORTNODUP", ds)
  gene2tadDT_file <- file.path("gene_data_final", ds, "genes2tad", "all_genes_positions.txt")
} else {
  stopifnot(length(args) == 2)
  gene2tadDT_file <- args[1]
  outFold <- args[2]
}

dir.create(outFold, recursive=TRUE)
stopifnot(file.exists(gene2tadDT_file))

gene2tadDT <- read.delim(gene2tadDT_file, header = F, stringsAsFactors = F, col.names = c("entrezID", "chromo", "start", "end", "region"))
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
stopifnot(!any(duplicated(gene2tadDT$entrezID)))
gene2tadDT <- gene2tadDT[grep("_TAD", gene2tadDT$region),]


all_tads <- unique(gene2tadDT$region)

all_TAD_pairs <- foreach(tad = all_tads, .combine='rbind') %dopar% {
  tad_g2t_dt <- gene2tadDT[gene2tadDT$region == tad,]
  if(nrow(tad_g2t_dt) == 1) return(NULL)
  # UPDATE 30.06.2018 -> ENSURE AS.CHARACTER + ALPHABETICAL ORDER !!!
  tad_g2t_dt$entrezID <- as.character(tad_g2t_dt$entrezID)
  tad_g2t_dt <- tad_g2t_dt[order(tad_g2t_dt$entrezID),]
  tadDT <- as.data.frame(t(combn(tad_g2t_dt$entrezID, m=2)))
  colnames(tadDT) <- c("gene1", "gene2")
  tadDT$region <- tad
  tadDT$gene1 <- as.character(tadDT$gene1)
  tadDT$gene2 <- as.character(tadDT$gene2)
  stopifnot(tadDT$gene1 < tadDT$gene2)
  tadDT
}

outFile <- file.path(outFold, "all_TAD_pairs.Rdata")
save(all_TAD_pairs, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))