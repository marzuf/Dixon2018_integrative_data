# look at the number of families by TADs for the consensus and the tissue TADs

# UPDATE 08.01.2019 => with the tissue specific family files !

# Rscript cmp_tissue_consensus_fam_TADs_otherFamFile.R

SSHFS <- F

setDir <- ifelse(SSHFS, "/media/electron", "")
setDir <- ifelse(SSHFS, "~/media/electron", "")

if(SSHFS) setwd("~/media/electron/mnt/etemp/marie/Dixon2018_integrative_data")

source("utils_fct.R")

plotType <- "svg"
myHeight <- ifelse(plotType=="svg", 7, 300) 
myWidth <- myHeight * 1.2 

outFolder <- "CMP_TISSUE_CONSENSUS_FAM"
dir.create(outFolder, recursive=TRUE)

other_tads_folder <- "gene_data_final"
stopifnot(dir.exists(other_tads_folder))

family_tads_folder <- "PREP_GENE_FAMILIES_TAD_DATA"
stopifnot(dir.exists(family_tads_folder))

familyData <- "hgnc"

# gene_data_final/GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal/genes2tad/all_assigned_regions.txt
# gene_data_final/GSE105318_ENCFF439QFU_DLD1/genes2tad/all_assigned_regions.txt
# gene_data_final/GSE105465_ENCFF777DUA_Caki2_vs_GSE105235_ENCFF235TGH_G401/genes2tad/all_assigned_regions.txt
# gene_data_final/GSE105566_ENCFF358MNA_Panc1/genes2tad/all_assigned_regions.txt
# gene_data_final/GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460/genes2tad/all_assigned_regions.txt
# gene_data_final/GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5/genes2tad/all_assigned_regions.txt
# gene_data_final/GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP/genes2tad/all_assigned_regions.txt
# gene_data_final/GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP_vs_GSE105697_ENCFF364CWZ_T47D/genes2tad/all_assigned_regions.txt
tissue_tadPos_files <- list.files(other_tads_folder, recursive = T, pattern="all_assigned_regions.txt")
stopifnot(length(tissue_tadPos_files) > 0)

tissue_g2t_files <- list.files(other_tads_folder, recursive = T, pattern="all_genes_positions.txt", full.names = TRUE)
stopifnot(length(tissue_g2t_files) > 0)

tissue_familyFiles <- list.files(family_tads_folder, recursive=T, pattern=paste0(familyData, "_entrezID_family_TAD_DT.Rdata"), full.names=TRUE)
stopifnot(length(tissue_familyFiles) > 0)

consensus_tadPosFile <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")    
consensus_tadPosDT <- read.delim(consensus_tadPosFile, header=F, stringsAsFactors = FALSE,
                                 col.names=c("chromo", "region", "start", "end"))
head(consensus_tadPosDT)
# chromo       region   start     end
# 1  chr10 chr10_BOUND1       1   80000

consensus_g2tFile <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt")    
consensus_g2t_DT <- read.delim(consensus_g2tFile, header=F, stringsAsFactors = FALSE,
                                 col.names=c("entrezID", "chromo", "start", "end", "region"))
consensus_g2t_DT$entrezID <- as.character(consensus_g2t_DT$entrezID)
head(consensus_g2t_DT)
# entrezID chromo  start    end     region
# 1    347688  chr10  92828  95178 chr10_TAD1
# 2    439945  chr10 125916 132183 chr10_TAD1

caller <- "TopDom"
inFoldFamily <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/", "PREP_GENE_FAMILIES_TAD_DATA"))
consensus_familyDT <- eval(parse(text = load(file.path(inFoldFamily, paste0(familyData, "_entrezID_family_TAD_DT.Rdata")))))
consensus_familyDT$entrezID <- as.character(consensus_familyDT$entrezID)
consensus_familyDT$region <- NULL
head(consensus_familyDT)
# entrezID chromo   start     end                                                                              hgnc_family                      hgnc_family_short
# 1    347688  chr10   92828   95178                                                                                 Tubulins                               Tubulins
# 3     10771  chr10  180405  300577 Zinc fingers MYND-type|PHD finger proteins|PWWP domain containing|Bromodomain containing                 Zinc fingers MYND-type

all_familyVars <- c("hgnc_family", "hgnc_family_short")
familyVar = "hgnc_family"

consensus_fam_DT <- merge(consensus_g2t_DT[,c("entrezID", "region")], consensus_familyDT[,c("entrezID", all_familyVars)], by="entrezID")
# entrezID      region                             hgnc_family                     hgnc_family_short
# 1         1 chr19_TAD74   Immunoglobulin like domain containing Immunoglobulin like domain containing
# 2        10  chr8_TAD21          Arylamine N-acetyltransferases        Arylamine N-acetyltransferases
# nbrFam_consensus_fam_DT<- consensus_fam_DT
# nbrFam_consensus_fam_DT$entrezID <- NULL
# nbrFam_consensus_fam_DT <- aggregate(. ~ region, data=nbrFam_consensus_fam_DT, FUN=function(x) length(unique(x)))

hgncFam_by_consensusTAD_DT <- aggregate(hgnc_family ~ region, data=consensus_fam_DT, FUN=function(x) length(unique(x)))
hgncFamShort_by_consensusTAD_DT <- aggregate(hgnc_family_short ~ region, data=consensus_fam_DT, FUN=function(x) length(unique(x)))
# nbrFam_consensus_fam_DT <- merge(hgncFam_by_consensusTAD_DT, hgncFamShort_by_consensusTAD_DT, by = "region")
# stopifnot(!is.na(nbrFam_consensus_fam_DT))

nEntrez_by_consensusTAD_DT <- aggregate(entrezID ~ region, data=consensus_fam_DT, FUN=function(x) length(unique(x)))


all_nbrs <- list()
all_nbrs[["pipelineConsensus"]] <- list(
  hgnc_fam = hgncFam_by_consensusTAD_DT$hgnc_family,
  hgnc_fam_short = hgncFamShort_by_consensusTAD_DT$hgnc_family_short,
  nEntrez = nEntrez_by_consensusTAD_DT$entrezID
) 

tissueFile = tissue_g2t_files[1]
tissueFile = tissue_g2t_files[7]

for(tissueFile in tissue_g2t_files) {
  
  curr_ds <- basename(dirname((dirname(tissueFile))))
  
  
  
  tissue_familyFile <- tissue_familyFiles[grepl(paste0("/", curr_ds, "/"), tissue_familyFiles)]
  
  cat("tissue_familyFile:\n", tissue_familyFile, "\n")
  
  
  stopifnot(length(tissue_familyFile) == 1)
  
  # cat("curr_ds:\n", curr_ds, "\n")
  # cat("tissue_familyFile:\n", tissue_familyFile, "\n")
  # cat("tissue_familyFiles:\n", tissue_familyFiles, "\n")
  # 
  # 
  # stopifnot(length(tissue_familyFile) == 1)
  # stopifnot(file.exists(tissue_familyFile))
  # cat(tissue_familyFile)
  tissue_familyDT <- eval(parse(text = load(tissue_familyFile)))
  tissue_familyDT$entrezID <- as.character(tissue_familyDT$entrezID)
  
  tissue_g2t_DT <- read.delim(tissueFile, header=F, stringsAsFactors = FALSE,
                                 col.names=c("entrezID", "chromo", "start", "end", "region"))
  tissue_g2t_DT$entrezID <- as.character(tissue_g2t_DT$entrezID)
  head(tissue_g2t_DT)
  
  tissue_fam_DT <- merge(tissue_g2t_DT[,c("entrezID", "region")], tissue_familyDT[,c("entrezID", all_familyVars)], by="entrezID")
  # entrezID      region                             hgnc_family                     hgnc_family_short
  # 1         1 chr19_TAD74   Immunoglobulin like domain containing Immunoglobulin like domain containing
  # 2        10  chr8_TAD21          Arylamine N-acetyltransferases        Arylamine N-acetyltransferases
  # nbrFam_tissue_fam_DT<- tissue_fam_DT
  # nbrFam_tissue_fam_DT$entrezID <- NULL
  # nbrFam_tissue_fam_DT <- aggregate(. ~ region, data=nbrFam_tissue_fam_DT, FUN=function(x) length(unique(x)))
  
  hgncFam_by_tissueTAD_DT <- aggregate(hgnc_family ~ region, data=tissue_fam_DT, FUN=function(x) length(unique(x)))
  hgncFamShort_by_tissueTAD_DT <- aggregate(hgnc_family_short ~ region, data=tissue_fam_DT, FUN=function(x) length(unique(x)))
  # nbrFam_tissue_fam_DT <- merge(hgncFam_by_tissueTAD_DT, hgncFamShort_by_tissueTAD_DT, by = "region")
  # stopifnot(!is.na(nbrFam_tissue_fam_DT))
  
  nEntrez_by_tissueTAD_DT <- aggregate(entrezID ~ region, data=tissue_fam_DT, FUN=function(x) length(unique(x)))
  
  all_nbrs[[paste0(curr_ds)]] <- list(
    hgnc_fam = hgncFam_by_tissueTAD_DT$hgnc_family,
    hgnc_fam_short = hgncFamShort_by_tissueTAD_DT$hgnc_family_short,
    nEntrez = nEntrez_by_tissueTAD_DT$entrezID
  ) 
  
}


outFile <- file.path(outFolder, "all_nbrs.Rdata")
save(all_nbrs, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

str(all_nbrs)

tissue_nbrs <- all_nbrs[names(all_nbrs) != "pipelineConsensus"]
str(tissue_nbrs)

pipConsensus_nbrFams <- all_nbrs[["pipelineConsensus"]][["hgnc_fam"]]
pipConsensus_nbrShortFams <- all_nbrs[["pipelineConsensus"]][["hgnc_fam_short"]]
pipConsensus_nbrGenes <- all_nbrs[["pipelineConsensus"]][["nEntrez"]]

head(pipConsensus_nbrFams)
tissue_nbrFams <- unlist(unname(lapply(tissue_nbrs, function(x)x[["hgnc_fam"]])))
tissue_nbrShortFams <- unlist(unname(lapply(tissue_nbrs, function(x)x[["hgnc_fam_short"]])))
tissue_nbrGenes <- unlist(unname(lapply(tissue_nbrs, function(x)x[["nEntrez"]])))

head(unlist(tissue_nbrFams))

outFile <- file.path(outFolder, paste0("nbrFams_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(list(
  nbrFams_pipConsensus=log10(pipConsensus_nbrFams),
  nbrFams_tissues = log10(tissue_nbrFams)
),
plotTit = paste0("# families by TAD"),
my_xlab = paste0("# fam. [log10]")
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("nbrShortFams_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(list(
  nbrShortFams_pipConsensus=log10(pipConsensus_nbrShortFams),
  nbrShortFams_tissues = log10(tissue_nbrShortFams)
),
plotTit = paste0("# short_families by TAD"),
my_xlab = paste0("# short_fam. [log10]")
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("nbrGenes_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(list(
  nbrGenes_pipConsensus=log10(pipConsensus_nbrGenes),
  nbrGenes_tissues = log10(tissue_nbrGenes)
),
plotTit = paste0("# genes by TAD"),
my_xlab = paste0("# genes. [log10]")
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))







