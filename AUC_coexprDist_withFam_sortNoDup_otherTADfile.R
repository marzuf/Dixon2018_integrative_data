startTime <- Sys.time()
cat(paste0("> Rscript AUC_coexprDist_withFam_sortNoDup_otherTADfile.R\n"))

stop("!!! DON'T USE ME - USE AUC_coexprDist_withFam_sortNoDup_otherTADfile_otherFamFile.R instead \n")

options(scipen=100)

printAndLog <- function(text, logFile = ""){
  cat(text)
  cat(text, append =T , file = logFile)
}

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggstatsplot, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

axisLabSize <- 12
legendSize <- 10
plotTitSize <- 14

mytheme <- theme(
  # top, right, bottom and left
  plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
  plot.title = element_text(hjust = 0.5, face = "bold", size=plotTitSize, vjust=1),
  plot.subtitle = element_text(hjust = 0.5, face = "bold", size=plotTitSize-2, vjust=1),
  panel.background = element_rect(fill = "white", colour = NA), 
  panel.border = element_rect(fill = NA, colour = "grey20"), 
  panel.grid.major = element_line(colour = "grey92"), 
  panel.grid.minor = element_line(colour = "grey92", size = 0.25), 
  strip.background = element_rect(fill = "grey85", colour = "grey20"), 
  #legend.key = element_rect(fill = "white", colour = NA), 
  axis.line.x = element_line(size = .3, color = "black"),
  axis.line.y = element_line(size = .3, color = "black"),
  axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=axisLabSize),
  axis.text.x = element_text(color="black", hjust=0.5,vjust = 1, size=axisLabSize),
  axis.title.y = element_text(color="black", size=axisLabSize+1),
  axis.title.x = element_text(color="black", size=axisLabSize+1),
  legend.text = element_text(size=legendSize),
  legend.key.height = unit(1.5,"cm"),
  legend.key = element_blank()
)


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

### HARD CODED
caller <- "TopDom"
corMethod <- "pearson"
buildTable <- TRUE
# for plotting:
# look at coexpression ~ distance up to distLimit bp
distLimit <- 500 * 10^3
fitMeth <- "loess"

# nbr of points for loess fit to take the AUC
nbrLoessPoints <- 1000

scatterFontSizeLabel <- 14
scatterFontSizeTitle <- 12

# UPDATE 30.06.2018:
# -> check that always $gene1 < $gene2 before left_join !!!


### RETRIEVE FROM COMMAND LINE

# Rscript AUC_coexprDist_withFam_sortNoDup_otherTADfile.R  <TADlist> <dataset> 

# Rscript AUC_coexprDist_withFam_sortNoDup_otherTADfile.R <TADlist> <dataset> <family>

# Rscript AUC_coexprDist_withFam_sortNoDup_otherTADfile.R GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal TCGAgbm_classical_mesenchymal


#### !!! change for otherTADfile: <<< GENE DATA DO NOT CHANGE
# retrieve the sameTAD data frame from:
# file.path("CREATE_SAME_TAD_SORTNODUP", curr_TADlist, "all_TAD_pairs.Rdata")

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 2 | length(args) == 3)

if(length(args) == 3) {
  txt <- paste0("> Parameters retrieved from command line:\n")
  stopifnot(length(args) == 2)
  curr_TADlist <- args[1]
  curr_dataset <- args[2]
  familyData <- args[3]
} else if(length(args) == 2){
  curr_TADlist <- args[1]
  curr_dataset <- args[2]
  txt <- paste0("> Default parameters:\n")
  familyData <- "hgnc"
}

outFold <- file.path("AUC_COEXPRDIST_WITHFAM_SORTNODUP",  curr_TADlist, paste0(curr_dataset,  "_", familyData))
dir.create(outFold, recursive = TRUE)

logFile <- file.path(outFold, paste0("coexpr_dist_withFam_otherTADfile_logFile.txt"))  
system(paste0("rm -f ", logFile))

printAndLog(txt, logFile)
txt <- paste0("... curr_TADlist = ",  curr_TADlist, "\n")
printAndLog(txt, logFile)
txt <- paste0("... curr_dataset = ",  curr_dataset, "\n")
printAndLog(txt, logFile)
txt <- paste0("... familyData = ",  familyData, "\n")
printAndLog(txt, logFile)

txt <- paste0("> ! Hard-coded parameters:\n")
printAndLog(txt, logFile)
txt <- paste0("... caller = ",  caller, "\n")
printAndLog(txt, logFile)
txt <- paste0("... corMethod = ",  corMethod, "\n")
printAndLog(txt, logFile)
txt <- paste0("... buildTable = ",  as.character(buildTable), "\n")
printAndLog(txt, logFile)
txt <- paste0("... distLimit = ",  distLimit, "\n")
printAndLog(txt, logFile)
txt <- paste0("... fitMeth = ",  fitMeth, "\n")
printAndLog(txt, logFile)

# mycols <- c("sameTAD" ="darkorange1" , "diffTAD"="darkslateblue",  "sameFam+sameTAD"="violetred1", "sameFam+diffTAD" = "lightskyblue")
mycols <- c("same TAD" ="darkorange1" , "diff. TAD"="darkslateblue",  "same Fam. + same TAD"="violetred1", "same Fam. + diff. TAD" = "lightskyblue")

sameTADcol <- mycols["same TAD"]
diffTADcol <- mycols["diff. TAD"]
sameFamSameTADcol <- mycols["same Fam. + same TAD"]
sameFamDiffTADcol <- mycols["same Fam. + diff. TAD"]

plotType <- "svg"
# myHeight <- ifelse(plotType == "png", 400, 7)
# myWidth <- ifelse(plotType == "png", 600, 10)
myHeight <- ifelse(plotType == "png", 200, 5)
myWidth <- ifelse(plotType == "png", 350, 6)

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

utilsDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg")
source(file.path(utilsDir, "coreg_utils_ggscatterhist.R"))

distFile <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/CREATE_DIST_SORTNODUP/all_dist_pairs.Rdata"))
stopifnot(file.exists(distFile))

coexprFile <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/CREATE_COEXPR_SORTNODUP"),  paste0(curr_dataset, "_", corMethod), "coexprDT.Rdata")
stopifnot(file.exists(coexprFile))

sameTADfile <- file.path("CREATE_SAME_TAD_SORTNODUP", curr_TADlist, "all_TAD_pairs.Rdata")
stopifnot(file.exists(sameTADfile))


dataset_pipDir <- file.path("PIPELINE", "OUTPUT_FOLDER", curr_TADlist, curr_dataset) # used to retrieve gene list
stopifnot(dir.exists(dataset_pipDir))

################################################ DATA PREPARATION

if(buildTable) {
  
  cat(paste0("... load DIST data\t", Sys.time(), "\t"))
  load(file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/CREATE_DIST_SORTNODUP/all_dist_pairs.Rdata")))
  cat(paste0(Sys.time(), "\n"))
  head(all_dist_pairs)
  nrow(all_dist_pairs)
  all_dist_pairs$gene1 <- as.character(all_dist_pairs$gene1)
  all_dist_pairs$gene2 <- as.character(all_dist_pairs$gene2)
  # UPDATE 30.06.2018
  stopifnot(all_dist_pairs$gene1 < all_dist_pairs$gene2)
  
  cat(paste0("... load TAD data\t", Sys.time(), "\t"))
  ### =>>> CHANGED HERE FOR OTHER TAD FILE !!!
  # load(file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/CREATE_SAME_TAD_SORTNODUP/all_TAD_pairs.Rdata")))
  load(file.path("CREATE_SAME_TAD_SORTNODUP", curr_TADlist, "all_TAD_pairs.Rdata"))
  cat(paste0(Sys.time(), "\n"))
  head(all_TAD_pairs)
  nrow(all_TAD_pairs)
  all_TAD_pairs$gene1 <- as.character(all_TAD_pairs$gene1)
  all_TAD_pairs$gene2 <- as.character(all_TAD_pairs$gene2)
  # UPDATE 30.06.2018
  stopifnot(all_TAD_pairs$gene1 < all_TAD_pairs$gene2)
  
  cat(paste0("... load COEXPR data\t", Sys.time(), "\t"))
  load(file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/CREATE_COEXPR_SORTNODUP"),  paste0(curr_dataset, "_", corMethod), "coexprDT.Rdata"))
  cat(paste0(Sys.time(), "\n"))
  head(coexprDT)
  nrow(coexprDT)
  coexprDT$gene1 <- as.character(coexprDT$gene1)
  coexprDT$gene2 <- as.character(coexprDT$gene2)
  all_TAD_pairs$gene2
  # UPDATE 30.06.2018
  stopifnot(coexprDT$gene1 < coexprDT$gene2)
  
  #============================== RETRIEVE PIPELINE DATA FOR THIS DATASET - USED ONLY FOR GENE LIST
  script0_name <- "0_prepGeneData"
  pipeline_geneList <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "pipeline_geneList.Rdata"))))
  pipeline_geneList <- as.character(pipeline_geneList)
  
  
  dataset_dist_pair <- all_dist_pairs[all_dist_pairs$gene1 %in% pipeline_geneList & 
                                        all_dist_pairs$gene2 %in% pipeline_geneList,]
  
  dataset_dist_pairs_limit <- dataset_dist_pair[dataset_dist_pair$dist <= distLimit,]
  head(dataset_dist_pairs_limit)
  nrow(dataset_dist_pairs_limit)
  
  dataset_TAD_pairs <- all_TAD_pairs[all_TAD_pairs$gene1 %in% pipeline_geneList & 
                                       all_TAD_pairs$gene2 %in% pipeline_geneList,]
  head(dataset_TAD_pairs)
  nrow(dataset_TAD_pairs)
  
  
  # START MERGING DATA 
  
  cat(paste0("... merge DIST - TAD data\t", Sys.time(), "\t"))
  dataset_dist_TAD_DT <- left_join(dataset_dist_pairs_limit, dataset_TAD_pairs, by=c("gene1", "gene2"))
  cat(paste0(Sys.time(), "\n"))
  
  dataset_dist_TAD_DT$sameTAD <- ifelse(is.na(dataset_dist_TAD_DT$region), 0, 1)
  
}

# all_familyData <- paste0(familyData, c("_family", "_family_short"))
all_familyData <- paste0(familyData, c( "_family_short"))

outFold_save <- outFold

for(i_fam in all_familyData) {
  
  
  outFold <- file.path(outFold_save, i_fam)
  dir.create(outFold, recursive = TRUE)
  
  if(buildTable){
    
    
    cat(paste0("... load FAMILY data\t", Sys.time(), "\t"))
    load(file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/CREATE_SAME_FAMILY_SORTNODUP/", i_fam, "_all_family_pairs.Rdata")))
    cat(paste0(Sys.time(), "\n"))
    head(all_family_pairs)
    nrow(all_family_pairs)
    all_family_pairs$gene1 <- as.character(all_family_pairs$gene1)
    all_family_pairs$gene2 <- as.character(all_family_pairs$gene2)
    stopifnot(all_family_pairs$gene1 < all_family_pairs$gene2)
    
    dataset_family_pairs <- all_family_pairs[all_family_pairs$gene1 %in% pipeline_geneList & 
                                               all_family_pairs$gene2 %in% pipeline_geneList,]
    head(dataset_family_pairs)
    
    cat(paste0("... merge FAMILY data\t", Sys.time(), "\t"))
    
    dataset_dist_TAD_fam_DT <- left_join(dataset_dist_TAD_DT, dataset_family_pairs, by=c("gene1", "gene2"))
    cat(paste0(Sys.time(), "\n"))
    
    dataset_dist_TAD_fam_DT$sameFamily <- ifelse(is.na(dataset_dist_TAD_fam_DT$family), 0, 1)
    
    cat(paste0("... merge COEXPR data\t", Sys.time(), "\t"))
    
    dataset_dist_TAD_fam_coexpr_DT <- left_join(dataset_dist_TAD_fam_DT, coexprDT, by=c("gene1", "gene2"))
    cat(paste0(Sys.time(), "\n"))
    
    allData_dt <- dataset_dist_TAD_fam_coexpr_DT
    allData_dt$region <- NULL
    allData_dt$family <- NULL
    allData_dt <- na.omit(allData_dt)
    
    # outFile <-file.path(outFold, paste0(i_fam, "_", "allData_dt.Rdata"))
    outFile <-file.path(outFold, paste0( "allData_dt.Rdata"))
    save(allData_dt, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))    
    
  } else{
    # outFile <-file.path(outFold, paste0(i_fam, "_", "allData_dt.Rdata"))
    outFile <-file.path(outFold, paste0( "allData_dt.Rdata"))
    load(outFile)
    # load("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/COEXPR_DIST_v3/TCGAcrc_msi_mss_hgnc/hgnc_family_allData_dt.Rdata")
    
  }
  
  nrow(allData_dt)
  allData_dt$dist_kb <- allData_dt$dist/1000
  
  allData_dt$curve1 <-  ifelse(allData_dt$sameTAD == "0", "diff. TAD", "same TAD")
  
  allData_dt$curve2 <-  ifelse(allData_dt$sameFamily == "0", NA,
                               ifelse(allData_dt$sameTAD == "0", "same Fam. + diff. TAD", "same Fam. + same TAD"))
  
  sameTAD_DT <- allData_dt[allData_dt$sameTAD == 1,c("gene1", "gene2", "coexpr", "dist", "dist_kb")]
  sameTAD_DT <- na.omit(sameTAD_DT)
  sameTAD_DT <- sameTAD_DT[order(sameTAD_DT$dist_kb),]
  # sameTAD_DT$cumdist <- cumsum(sameTAD_DT$dist_kb)
  sameTAD_DT$nPair <- 1:nrow(sameTAD_DT)
  sameTAD_DT$label <- "same TAD"
  
  diffTAD_DT <- allData_dt[allData_dt$sameTAD == 0,c("gene1", "gene2",  "coexpr", "dist", "dist_kb")]
  diffTAD_DT <- na.omit(diffTAD_DT)
  diffTAD_DT <- diffTAD_DT[order(diffTAD_DT$dist_kb),]
  # diffTAD_DT$cumdist <- cumsum(diffTAD_DT$dist_kb)
  diffTAD_DT$nPair <- 1:nrow(diffTAD_DT)
  diffTAD_DT$label <- "diff. TAD"
  
  sameFam_sameTAD_DT <- allData_dt[allData_dt$sameFamily == 1 & allData_dt$sameTAD == 1 ,c("gene1", "gene2", "coexpr", "dist",  "dist_kb")]
  sameFam_sameTAD_DT <- na.omit(sameFam_sameTAD_DT)
  sameFam_sameTAD_DT <- sameFam_sameTAD_DT[order(sameFam_sameTAD_DT$dist_kb),]
  # sameFam_sameTAD_DT$cumdist <- cumsum(sameFam_sameTAD_DT$dist_kb)
  sameFam_sameTAD_DT$nPair <- 1:nrow(sameFam_sameTAD_DT)
  sameFam_sameTAD_DT$label <- "same Fam. + same TAD"
  
  sameFam_diffTAD_DT <- allData_dt[allData_dt$sameFamily == 1 & allData_dt$sameTAD == 0 ,c("gene1", "gene2",  "coexpr", "dist", "dist_kb")]
  sameFam_diffTAD_DT <- na.omit(sameFam_diffTAD_DT)
  sameFam_diffTAD_DT <- sameFam_diffTAD_DT[order(sameFam_diffTAD_DT$dist_kb),]
  # sameFam_diffTAD_DT$cumdist <- cumsum(sameFam_diffTAD_DT$dist_kb)
  sameFam_diffTAD_DT$nPair <- 1:nrow(sameFam_diffTAD_DT)
  sameFam_diffTAD_DT$label <- "same Fam. + diff. TAD"
  
  stopifnot(is.numeric(sameTAD_DT$dist[1]))
  stopifnot(is.numeric(sameTAD_DT$coexpr[1]))
  stopifnot(is.numeric(diffTAD_DT$dist[1]))
  stopifnot(is.numeric(diffTAD_DT$coexpr[1]))
  stopifnot(is.numeric(sameFam_sameTAD_DT$dist[1]))
  stopifnot(is.numeric(sameFam_sameTAD_DT$coexpr[1]))
  stopifnot(is.numeric(sameFam_diffTAD_DT$dist[1]))
  stopifnot(is.numeric(sameFam_diffTAD_DT$coexpr[1]))
  
  #***
  
  if(fitMeth == "loess") {
    
    my_ylab <- paste0("Gene pair coexpression (", corMethod, ", qqnormDT)")
    my_xlab <- paste0("Distance between the 2 genes (kb)")
    my_sub <- paste0(curr_dataset)
    
    # PREDICT WITH ORIGINAL DISTANCE VALUES
    my_xlab <- paste0("Distance between the 2 genes (bp)")
    
    smooth_vals_sameTAD <- predict(loess(coexpr ~ dist, data = sameTAD_DT), sort(sameTAD_DT$dist))
    smooth_vals_diffTAD <- predict(loess(coexpr ~ dist, data = diffTAD_DT), sort(diffTAD_DT$dist))
    
    auc_diffTAD_obsDist <- auc(x = sort(diffTAD_DT$dist), y = smooth_vals_diffTAD)
    auc_sameTAD_obsDist <- auc(x = sort(sameTAD_DT$dist), y = smooth_vals_sameTAD)
    
    outFile <- file.path(outFold, paste0("sameTAD_diffTAD_loessFit_originalDist", ".", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    plot(NULL,
         xlim = range(allData_dt$dist), 
         ylim = range(c(smooth_vals_sameTAD, smooth_vals_diffTAD)),
         # xlab="", 
         # ylab="",
         xlab=my_xlab, 
         ylab=my_ylab,
         main=paste0(curr_dataset, ": coexpr ~ dist loess fit"))
    mtext(text = "observed distance values", side = 3)
    lines( x = sort(sameTAD_DT$dist), y = smooth_vals_sameTAD, col = sameTADcol)
    lines( x = sort(diffTAD_DT$dist), y = smooth_vals_diffTAD, col = diffTADcol)
    legend("topright", 
           legend=c(paste0("sameTAD\n(AUC=", round(auc_sameTAD_obsDist, 2), ")"), paste0("diffTAD\n(AUC=", round(auc_diffTAD_obsDist, 2))),  
           col = c(sameTADcol, diffTADcol),
           lty=1,
           bty = "n")
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    # PREDICT WITH DISTANCE VECTOR
    distVect <- seq(from=0, to = distLimit, length.out = nbrLoessPoints)
    smooth_vals_sameTAD_distVect <- predict(loess(coexpr ~ dist, data = sameTAD_DT), distVect)
    smooth_vals_diffTAD_distVect <- predict(loess(coexpr ~ dist, data = diffTAD_DT), distVect)
    
    auc_diffTAD_distVect <- auc(x = distVect, y = smooth_vals_diffTAD_distVect)
    auc_sameTAD_distVect <- auc(x = distVect, y = smooth_vals_sameTAD_distVect)
    
    diffTAD_mod <- loess(coexpr ~ dist, data = diffTAD_DT)
    outFile <- file.path(outFold, "diffTAD_mod.Rdata")
    save(diffTAD_mod, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    diffTAD_obsDist <- diffTAD_DT$dist
    outFile <- file.path(outFold, "diffTAD_obsDist.Rdata")
    save(diffTAD_obsDist, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, "auc_diffTAD_distVect.Rdata")
    save(auc_diffTAD_distVect, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, "auc_sameTAD_distVect.Rdata")
    save(auc_sameTAD_distVect, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, "smooth_vals_sameTAD_distVect.Rdata")
    save(smooth_vals_sameTAD_distVect, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, "smooth_vals_diffTAD_distVect.Rdata")
    save(smooth_vals_diffTAD_distVect, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, "distVect.Rdata")
    save(distVect, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, paste0("sameTAD_diffTAD_loessFit_vectDist.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    plot(NULL,
         xlim = range(distVect), 
         ylim = range(c(na.omit(smooth_vals_sameTAD_distVect), na.omit(smooth_vals_diffTAD_distVect))),
         # xlab="", 
         # ylab="",
         xlab=my_xlab,
         ylab=my_ylab,
         main=paste0(curr_dataset, ": coexpr ~ dist loess fit"))
    mtext(text = paste0("distance values seq from 0 to ", distLimit, " (# points = ", nbrLoessPoints, ")"), side = 3)
    lines( x = distVect, y = smooth_vals_sameTAD_distVect, col = sameTADcol)
    lines( x = distVect, y = smooth_vals_diffTAD_distVect, col = diffTADcol)
    legend("topright", 
           legend=c(paste0("sameTAD\n(AUC=", round(auc_sameTAD_distVect, 2), ")"), paste0("diffTAD\n(AUC=", round(auc_diffTAD_distVect, 2))), 
           col = c(sameTADcol, diffTADcol),
           lty=1,
           bty = "n")
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    

    
    
    
    ################################ DO THE SAME FOR SAME FAM SAME TAD VS. SAME FAM DIFF TAD
    
    
    smooth_vals_sameFamSameTAD <- predict(loess(coexpr ~ dist, data = sameFam_sameTAD_DT), sort(sameFam_sameTAD_DT$dist))
    smooth_vals_sameFamDiffTAD <- predict(loess(coexpr ~ dist, data = sameFam_diffTAD_DT), sort(sameFam_diffTAD_DT$dist))
    
    auc_sameFamDiffTAD_obsDist <- auc(x = sort(sameFam_diffTAD_DT$dist), y = smooth_vals_sameFamDiffTAD)
    auc_sameFamSameTAD_obsDist <- auc(x = sort(sameFam_sameTAD_DT$dist), y = smooth_vals_sameFamSameTAD)
    
    outFile <- file.path(outFold, paste0("sameFamsameTAD_sameFamDiffTAD_loessFit_originalDist", ".", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    plot(NULL,
         xlim = range(allData_dt$dist), 
         ylim = range(c(smooth_vals_sameFamSameTAD, smooth_vals_sameFamDiffTAD)),
         # xlab="", 
         # ylab="",
         xlab=my_xlab, 
         ylab=my_ylab,
         main=paste0(curr_dataset, ": coexpr ~ dist loess fit"))
    mtext(text = "observed distance values", side = 3)
    lines( x = sort(sameFam_sameTAD_DT$dist), y = smooth_vals_sameFamSameTAD, col = sameFamSameTADcol)
    lines( x = sort(sameFam_diffTAD_DT$dist), y = smooth_vals_sameFamDiffTAD, col = sameFamDiffTADcol)
    legend("topright", 
           legend=c(paste0("sameFamSameTAD\n(AUC=", round(auc_sameFamSameTAD_obsDist, 2), ")"), 
                    paste0("sameFamDiffTAD\n(AUC=", round(auc_sameFamDiffTAD_obsDist, 2))),  
           col = c(sameFamSameTADcol, sameFamDiffTADcol),
           lty=1,
           bty = "n")
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    # PREDICT WITH DISTANCE VECTOR
    smooth_vals_sameFamSameTAD_distVect <- predict(loess(coexpr ~ dist, data = sameFam_sameTAD_DT), distVect)
    smooth_vals_sameFamDiffTAD_distVect <- predict(loess(coexpr ~ dist, data = sameFam_diffTAD_DT), distVect)
    
    auc_sameFamDiffTAD_distVect <- auc(x = distVect, y = smooth_vals_sameFamDiffTAD_distVect)
    auc_sameFamSameTAD_distVect <- auc(x = distVect, y = smooth_vals_sameFamSameTAD_distVect)
    
    sameFamDiffTAD_mod <- loess(coexpr ~ dist, data = sameFam_diffTAD_DT)
    outFile <- file.path(outFold, "sameFamDiffTAD_mod.Rdata")
    save(sameFamDiffTAD_mod, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    sameFamDiffTAD_obsDist <- sameFam_diffTAD_DT$dist
    outFile <- file.path(outFold, "sameFamDiffTAD_obsDist.Rdata")
    save(sameFamDiffTAD_obsDist, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, "auc_sameFamDiffTAD_distVect.Rdata")
    save(auc_sameFamDiffTAD_distVect, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, "auc_sameFamSameTAD_distVect.Rdata")
    save(auc_sameFamSameTAD_distVect, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, "smooth_vals_sameFamSameTAD_distVect.Rdata")
    save(smooth_vals_sameFamSameTAD_distVect, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, "smooth_vals_sameFamDiffTAD_distVect.Rdata")
    save(smooth_vals_sameFamDiffTAD_distVect, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, paste0("sameFamSameTAD_sameFamDiffTAD_loessFit_vectDist.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    plot(NULL,
         xlim = range(distVect), 
         ylim = range(c(na.omit(smooth_vals_sameFamSameTAD_distVect), na.omit(smooth_vals_sameFamDiffTAD_distVect))),
         # xlab="", 
         # ylab="",
         xlab=my_xlab,
         ylab=my_ylab,
         main=paste0(curr_dataset, ": coexpr ~ dist loess fit"))
    mtext(text = paste0("distance values seq from 0 to ", distLimit, " (# points = ", nbrLoessPoints, ")"), side = 3)
    lines( x = distVect, y = smooth_vals_sameFamSameTAD_distVect, col = sameFamSameTADcol)
    lines( x = distVect, y = smooth_vals_sameFamDiffTAD_distVect, col = sameFamDiffTADcol)
    legend("topright", 
           legend=c(paste0("sameFamSameTAD\n(AUC=", round(auc_sameFamSameTAD_distVect, 2), ")"), 
                    paste0("sameFamDiffTAD\n(AUC=", round(auc_sameFamDiffTAD_distVect, 2))), 
           col = c(sameFamSameTADcol, sameFamDiffTADcol),
           lty=1,
           bty = "n")
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    ################################################################ 
    outFile <- file.path(outFold, paste0("sameTAD_diffTAD_sameFamSameTAD_sameFamDiffTAD_loessFit_vectDist.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    plot(NULL,
         xlim = range(distVect), 
         ylim = range(c(na.omit(smooth_vals_sameTAD_distVect), 
                        na.omit(smooth_vals_sameFamSameTAD_distVect),
                        na.omit(smooth_vals_diffTAD_distVect),
                        na.omit(smooth_vals_sameFamDiffTAD_distVect))),
         # xlab="", 
         # ylab="",
         xlab=my_xlab,
         ylab=my_ylab,
         main=paste0(curr_dataset, ": coexpr ~ dist loess fit"))
    mtext(text = paste0("distance values seq from 0 to ", distLimit, " (# points = ", nbrLoessPoints, ")"), side = 3)
    lines( x = distVect, y = smooth_vals_sameFamSameTAD_distVect, col = sameFamSameTADcol)
    lines( x = distVect, y = smooth_vals_sameFamDiffTAD_distVect, col = sameFamDiffTADcol)
    lines( x = distVect, y = smooth_vals_sameTAD_distVect, col = sameTADcol)
    lines( x = distVect, y = smooth_vals_diffTAD_distVect, col = diffTADcol)
    legend("topright", 
           legend=c(paste0("sameTAD\n(AUC=", round(auc_sameTAD_distVect, 2), ")"), 
                    paste0("diffTAD\n(AUC=", round(auc_diffTAD_distVect, 2)), 
                     paste0("sameFamSameTAD\n(AUC=", round(auc_sameFamSameTAD_distVect, 2), ")"), 
                     paste0("sameFamDiffTAD\n(AUC=", round(auc_sameFamDiffTAD_distVect, 2))), 
           col = c(sameTADcol, diffTADcol,
                   sameFamSameTADcol, sameFamDiffTADcol),
           lty=1,
           bty = "n")
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    ################################################################ 
    auc_values <- list(
      auc_diffTAD_distVect = auc_diffTAD_distVect,
      auc_sameTAD_distVect = auc_sameTAD_distVect,
      auc_ratio_same_over_diff_distVect = auc_sameTAD_distVect/auc_diffTAD_distVect,
      auc_diffTAD_obsDist = auc_diffTAD_obsDist,
      auc_sameTAD_obsDist = auc_sameTAD_obsDist,
      auc_ratio_same_over_diff_obsDist = auc_sameTAD_distVect/auc_diffTAD_obsDist,
      
      auc_sameFamDiffTAD_distVect = auc_sameFamDiffTAD_distVect,
      auc_sameFamSameTAD_distVect = auc_sameFamSameTAD_distVect,
      auc_ratio_sameFam_same_over_diff_distVect = auc_sameFamSameTAD_distVect/auc_sameFamDiffTAD_distVect,
      auc_sameFamDiffTAD_obsDist = auc_sameFamDiffTAD_obsDist,
      auc_sameFamSameTAD_obsDist = auc_sameFamSameTAD_obsDist,
      auc_ratio_sameFam_same_over_diff_obsDist = auc_sameFamSameTAD_distVect/auc_sameFamDiffTAD_obsDist
      
    )
    
    outFile <- file.path(outFold, paste0("auc_values.Rdata"))
    save(auc_values, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
  } else{
    stop("only loess implemented yet\n")
    
  }
} # end iterating over family data
  
  ######################################################################################
  ######################################################################################
  ######################################################################################
  cat(paste0("... written: ", logFile, "\n"))
  ######################################################################################
  cat("*** DONE\n")
  cat(paste0(startTime, "\n", Sys.time(), "\n"))