# AUC coexprDist: look at the number of gene pairs as a function of distance  families by TADs for the consensus and the tissue TADs

# Rscript cmp_nbrGenes_dist.R GSE105566_ENCFF358MNA_Panc1 TCGApaad_wt_mutKRAS

SSHFS <- F

setDir <- ifelse(SSHFS, "/media/electron", "")
setDir <- ifelse(SSHFS, "~/media/electron", "")

if(SSHFS) setwd("~/media/electron/mnt/etemp/marie/Dixon2018_integrative_data")

source("utils_fct.R")

plotType <- "png"
myHeight <- ifelse(plotType=="svg", 7, 400) 
myWidth <- myHeight * 1.2 

TAD_ds <- "GSE105566_ENCFF358MNA_Panc1"
expr_ds <- "TCGApaad_wt_mutKRAS"
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
TAD_ds <- args[1]
expr_ds <- args[2]

outFold <- file.path("CMP_NBRGENES_DIST", TAD_ds, expr_ds)
dir.create(outFold, recursive=TRUE)

coexprFile_dataset <- file.path("CREATE_COEXPR_SORTNODUP", TAD_ds, paste0(expr_ds, "_pearson"), "coexprDT.Rdata")
stopifnot(file.exists(coexprFile_dataset))

distFile_dataset <- file.path("CREATE_DIST_SORTNODUP", TAD_ds, "all_dist_pairs.Rdata")
stopifnot(file.exists(distFile_dataset))

coexprFile_consensus <- file.path(setDir,
                                  "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom",
                                  "CREATE_COEXPR_SORTNODUP",
                                  paste0(expr_ds, "_pearson"), "coexprDT.Rdata")
stopifnot(file.exists(coexprFile_consensus))

distFile_consensus <- file.path(setDir,
                                "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom",
                                "CREATE_DIST_SORTNODUP",
                                "all_dist_pairs.Rdata")
stopifnot(file.exists(distFile_consensus))

cat("... load dataset_coexpr_DT \n")
dataset_coexpr_DT <- eval(parse(text = load(coexprFile_dataset)))
cat("... load dataset_dist_DT \n")
dataset_dist_DT <- eval(parse(text = load(distFile_dataset)))

dataset_coexpr_DT$gene1 <- as.character(dataset_coexpr_DT$gene1)
dataset_coexpr_DT$gene2 <- as.character(dataset_coexpr_DT$gene2)
dataset_dist_DT$gene1 <- as.character(dataset_dist_DT$gene1)
dataset_dist_DT$gene2 <- as.character(dataset_dist_DT$gene2)

cat("... merge into dataset_coexprDist_DT \n")
dataset_coexprDist_DT <- merge(dataset_coexpr_DT, dataset_dist_DT, by=c("gene1", "gene2"))

cat("... load consensus_coexpr_DT \n")
consensus_coexpr_DT <- eval(parse(text = load(coexprFile_consensus)))
cat("... load consensus_dist_DT \n")
consensus_dist_DT <- eval(parse(text = load(distFile_consensus)))

consensus_coexpr_DT$gene1 <- as.character(consensus_coexpr_DT$gene1)
consensus_coexpr_DT$gene2 <- as.character(consensus_coexpr_DT$gene2)
consensus_dist_DT$gene1 <- as.character(consensus_dist_DT$gene1)
consensus_dist_DT$gene2 <- as.character(consensus_dist_DT$gene2)

cat("... merge into consensus_coexprDist_DT \n")
consensus_coexprDist_DT <- merge(consensus_coexpr_DT, consensus_dist_DT, by=c("gene1", "gene2"))

write.table(consensus_coexprDist_DT[1:10,], file="", quote=F, sep="\n", col.names=T, row.names=F)

outFile <- file.path(outFold, "consensus_coexprDist_DT.Rdata")
save(consensus_coexprDist_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

# aggregate the number of genes by distance bin (using hist function)


