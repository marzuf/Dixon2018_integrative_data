SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

source("utils_fct.R")

if(SSHFS) setwd("/media/electron/mnt/etemp/marie/Dixon2018_integrative_data")

cat("> START compare_pipeline_results.R", "\n")
# Rscript compare_pipeline_results.R TCGAbrca_lum_bas GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP consensus

# pipeline gene list
# ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 
# "7105"         "64102"          "8813"         "57147"         "55732" 
exprDS <- "TCGAbrca_lum_bas"
ds1 <- "GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP"
ds2 <- "consensus"

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
exprDS <- args[1]
ds1 <- args[2]
ds2 <- args[3]

outFold <- file.path("COMPARE_PIPELINE_SIGNIF", paste0(ds1, "_", ds2), exprDS)
system(paste0("mkdir -p ", outFold))

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 400, 7)

logFile <- file.path(outFold, paste0("compare_pip_results_", ds1, "_", ds2, "_logFile.txt"))
if(SSHFS) logFile <- ""
if(!SSHFS) system(paste0("rm -f ", logFile))

############################## DATASET 1
stopifnot(ds1 != "consensus")
mainDir1 <- file.path(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data/PIPELINE/OUTPUT_FOLDER", ds1, exprDS)

g2t_file1 <- file.path(setDir,
                       "/mnt/etemp/marie/Dixon2018_integrative_data/gene_data_final",
                        ds1, "genes2tad/all_genes_positions.txt")
g2t1 <- read.delim(g2t_file1, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
g2t1$entrezID <- as.character(g2t1$entrezID)

geneList1 <- eval(parse(text = load(file.path(mainDir1, "0_prepGeneData", "pipeline_geneList.Rdata"))))


############################## DATASET 2 OR CONSENSUS

if(ds2 == "consensus") {
  mainDir2 <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER", exprDS)
  g2t_file2 <- file.path(setDir,
                         "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt")

} else {
  mainDir2 <- file.path(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data/PIPELINE/OUTPUT_FOLDER", ds2, exprDS)
  g2t_file2 <- file.path(setDir,
                         "/mnt/etemp/marie/Dixon2018_integrative_data/gene_data_final",
                         ds2, "genes2tad/all_genes_positions.txt")
}
g2t2 <- read.delim(g2t_file2, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
g2t2$entrezID <- as.character(g2t2$entrezID)
geneList2 <- eval(parse(text = load(file.path(mainDir2, "0_prepGeneData", "pipeline_geneList.Rdata"))))


######################################################################################################################## 

commonGenes <- intersect(geneList1, geneList2)

stopifnot(commonGenes %in% g2t1$entrezID)
stopifnot(commonGenes %in% g2t2$entrezID)

txt <- paste0("... genes from ds1:\t", length(geneList1), "\n")
printAndLog(txt, logFile)
txt <- paste0("... genes from ds2:\t", length(geneList2), "\n")
printAndLog(txt, logFile)
txt <- paste0("... genes in common:\t", length(commonGenes), "\n")
printAndLog(txt, logFile)

tadpval1 <- eval(parse(text = load(file.path(mainDir1, "11_runEmpPvalCombined", "emp_pval_combined.Rdata"))))
tadpval1 <- p.adjust(tadpval1, method ="BH")
tadpval1 <- sort(tadpval1, decreasing=FALSE)
tadrank1 <- rank(tadpval1, ties.method = "min")
stopifnot(names(tadrank1) %in% g2t1$region)
tadrank1_dt <- data.frame(region = names(tadrank1), region_rank = tadrank1, stringsAsFactors = FALSE)
rownames(tadrank1_dt) <- NULL

tadpval2 <- eval(parse(text = load(file.path(mainDir2, "11_runEmpPvalCombined", "emp_pval_combined.Rdata"))))
tadpval2 <- p.adjust(tadpval2, method ="BH")
tadpval2 <- sort(tadpval2, decreasing=FALSE)
tadrank2 <- rank(tadpval2, ties.method = "min")
stopifnot(names(tadrank2) %in% g2t2$region)
tadrank2_dt <- data.frame(region = names(tadrank2), region_rank = tadrank2, stringsAsFactors = FALSE)
rownames(tadrank2_dt) <- NULL

g2t1_dt <- g2t1[g2t1$entrezID %in% commonGenes,]
stopifnot(nrow(g2t1_dt) > 0)
g2t2_dt <- g2t2[g2t2$entrezID %in% commonGenes,]
stopifnot(nrow(g2t2_dt) > 0)


g2t1_dt_withRank <- merge(g2t1_dt, tadrank1_dt, by="region", all.x=TRUE, all.y=FALSE)
stopifnot(commonGenes %in% g2t1_dt_withRank$entrezID)
stopifnot(!duplicated(g2t1_dt_withRank$entrezID))
rownames(g2t1_dt_withRank) <- as.character(g2t1_dt_withRank$entrezID)
g2t1_dt_withRank <- g2t1_dt_withRank[as.character(commonGenes),]
stopifnot(g2t1_dt_withRank$entrezID == rownames(g2t1_dt_withRank))
stopifnot(g2t1_dt_withRank$entrezID == commonGenes)

g2t2_dt_withRank <- merge(g2t2_dt, tadrank2_dt, by="region", all.x=TRUE, all.y=FALSE)
stopifnot(commonGenes %in% g2t2_dt_withRank$entrezID)
stopifnot(!duplicated(g2t2_dt_withRank$entrezID))
rownames(g2t2_dt_withRank) <- as.character(g2t2_dt_withRank$entrezID)
g2t2_dt_withRank <- g2t2_dt_withRank[as.character(commonGenes),]
stopifnot(g2t2_dt_withRank$entrezID == rownames(g2t2_dt_withRank))
stopifnot(g2t2_dt_withRank$entrezID == commonGenes)

myTit <- paste0(exprDS, ": gene TAD ranking comparison")
myxlab <- paste0("gene TAD rank ", ds1, " (", length(commonGenes), "/", length(geneList1), ")")
myylab <- paste0("gene TAD rank ", ds2, " (", length(commonGenes), "/", length(geneList2), ")")
mySub <- paste0(ds1, " vs. ", ds2)
# mySubBot <- paste0("(# intersect genes = ", length(commonGenes), ")")
outFile <- file.path(outFold, paste0(ds1, "_", ds2, "_gene_tad_rank.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(x = g2t1_dt_withRank$region_rank,
     y = g2t2_dt_withRank$region_rank,
     xlab=myxlab,
     ylab=myylab,
     pch = 16, cex = 0.7,
     main = myTit
     )
mtext(side=3, text = mySub)

add_curv_fit(x = g2t1_dt_withRank$region_rank, 
             y=g2t2_dt_withRank$region_rank, withR2 = FALSE, lty=2, col="darkgray")
  
addCorr(x=g2t1_dt_withRank$region_rank, 
          y = g2t2_dt_withRank$region_rank,
        bty="n",
  legPos="bottomright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFold, paste0(ds1, "_", ds2, "_gene_tad_rank_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

densplot(x=g2t1_dt_withRank$region_rank,
         y=g2t2_dt_withRank$region_rank,
         xlab=myxlab,
         ylab=myylab,
         pch = 16, cex = 0.7,
         main = myTit
)
mtext(side=3, text = mySub)
add_curv_fit(x = g2t1_dt_withRank$region_rank, 
             y=g2t2_dt_withRank$region_rank, withR2 = FALSE, lty=2, col="darkgray")

addCorr(x=g2t1_dt_withRank$region_rank, 
        y = g2t2_dt_withRank$region_rank,
        bty="n",
        legPos="bottomright")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


