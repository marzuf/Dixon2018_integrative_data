# look at the number of families by TADs for the consensus and the tissue TADs

SSHFS <- TRUE

setDir <- ifelse(SSHFS, "/media/electron", "")

if(SSHFS) setwd("/media/electron/mnt/etemp/marie/Dixon2018_integrative_data")

other_tads_folder <- "gene_data_final"
stopifnot(dir.exists(other_tads_folder))

# gene_data_final/GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal/genes2tad/all_assigned_regions.txt
# gene_data_final/GSE105318_ENCFF439QFU_DLD1/genes2tad/all_assigned_regions.txt
# gene_data_final/GSE105465_ENCFF777DUA_Caki2_vs_GSE105235_ENCFF235TGH_G401/genes2tad/all_assigned_regions.txt
# gene_data_final/GSE105566_ENCFF358MNA_Panc1/genes2tad/all_assigned_regions.txt
# gene_data_final/GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460/genes2tad/all_assigned_regions.txt
# gene_data_final/GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5/genes2tad/all_assigned_regions.txt
# gene_data_final/GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP/genes2tad/all_assigned_regions.txt
# gene_data_final/GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP_vs_GSE105697_ENCFF364CWZ_T47D/genes2tad/all_assigned_regions.txt
tissue_tadPos_files <- list.files(other_tads_folder, recursive = T, patter="all_assigned_regions.txt")
stopifnot(length(tissue_tadPos_files) > 0)


tissue_g2t_files <- list.files(other_tads_folder, recursive = T, pattern="all_genes_positions.txt", full.names = TRUE)
stopifnot(length(tissue_g2t_files) > 0)



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


familyData <- "hgnc"
caller <- "TopDom"
inFoldFamily <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/", "PREP_GENE_FAMILIES_TAD_DATA"))
familyDT <- eval(parse(text = load(file.path(inFoldFamily, paste0(familyData, "_entrezID_family_TAD_DT.Rdata")))))
familyDT$entrezID <- as.character(familyDT$entrezID)
familyDT$region <- NULL
head(familyDT)
# entrezID chromo   start     end                                                                              hgnc_family                      hgnc_family_short
# 1    347688  chr10   92828   95178                                                                                 Tubulins                               Tubulins
# 3     10771  chr10  180405  300577 Zinc fingers MYND-type|PHD finger proteins|PWWP domain containing|Bromodomain containing                 Zinc fingers MYND-type


all_familyVars <- c("hgnc_family", "hgnc_family_short")
familyVar = "hgnc_family"


consensus_fam_DT <- merge(consensus_g2t_DT[,c("entrezID", "region")], familyDT[,c("entrezID", all_familyVars)], by="entrezID")
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
for(tissueFile in tissue_g2t_files) {
  
  tissue_g2t_DT <- read.delim(tissueFile, header=F, stringsAsFactors = FALSE,
                                 col.names=c("entrezID", "chromo", "start", "end", "region"))
  tissue_g2t_DT$entrezID <- as.character(tissue_g2t_DT$entrezID)
  head(tissue_g2t_DT)
  
  tissue_fam_DT <- merge(tissue_g2t_DT[,c("entrezID", "region")], familyDT[,c("entrezID", all_familyVars)], by="entrezID")
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
  
  curr_ds <- basename(dirname((dirname(tissueFile))))
  
  all_nbrs[[paste0(curr_ds)]] <- list(
    hgnc_fam = hgncFam_by_tissueTAD_DT$hgnc_family,
    hgnc_fam_short = hgncFamShort_by_tissueTAD_DT$hgnc_family_short,
    nEntrez = nEntrez_by_tissueTAD_DT$entrezID
  ) 
  
}


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

plot_multiDens(list(
  nbrFams_pipConsensus=log10(pipConsensus_nbrFams),
  nbrFams_tissues = log10(tissue_nbrFams)
))


plot_multiDens(list(
  nbrShortFams_pipConsensus=log10(pipConsensus_nbrShortFams),
  nbrShortFams_tissues = log10(tissue_nbrShortFams)
))


plot_multiDens(list(
  nbrGenes_pipConsensus=log10(pipConsensus_nbrGenes),
  nbrGenes_tissues = log10(tissue_nbrGenes)
))







