
# in this file, settings that are specific for a run on a dataset

# gives path to output folder
pipOutFold <- "OUTPUT_FOLDER/TCGAluad_luad_mutKRAS"

# full path (starting with /mnt/...)
# following format expected for the input
# colnames = samplesID
# rownames = geneID
# !!! geneID are expected not difficulted

# *************************************************************************************************************************
# ************************************ SETTINGS FOR 0_prepGeneData
# *************************************************************************************************************************

rnaseqDT_file <- "/mnt/ed4/marie/other_datasets/TCGA_luad_luad_mutKRAS/exprDT.Rdata"
my_sep <- "\t"
# input is Rdata or txt file ?
# TRUE if the input is Rdata
inRdata <- TRUE

# can be ensemblID, entrezID, geneSymbol
geneID_format <- "geneSymbol"
stopifnot(geneID_format %in% c("ensemblID", "entrezID", "geneSymbol"))

# are geneID rownames ? -> "rn" or numeric giving the column
geneID_loc <- "rn"
stopifnot(geneID_loc == "rn" | is.numeric(geneID_loc))

removeDupGeneID <- TRUE

# *************************************************************************************************************************
# ************************************ SETTINGS FOR 1_runGeneDE
# *************************************************************************************************************************

# labels for conditions
cond1 <- "wt"
cond2 <- "mut"

# path to sampleID for each condition - should be Rdata ( ! sample1 for cond1, sample2 for cond2 ! )
sample1_file <- "/mnt/ed4/marie/other_datasets/TCGA_luad_luad_mutKRAS/wt_ID.Rdata"
sample2_file <- "/mnt/ed4/marie/other_datasets/TCGA_luad_luad_mutKRAS/mut_ID.Rdata"


minCpmRatio <- 20/888 

inputDataType <- "RSEM"

# > file edited: Fri, 16 Nov 2018 15:21:31 +0100 

# path to output folder:
pipOutFold <- "/mnt/etemp/marie/Dixon2018_integrative_data/PIPELINE/OUTPUT_FOLDER/GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460/TCGAluad_luad_mutKRAS"

# OVERWRITE THE DEFAULT SETTINGS FOR INPUT FILES - use TADs from the current Hi-C dataset 
TADpos_file <- paste0(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data/gene_data_final/GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460/genes2tad/all_assigned_regions.txt")
#chr1    chr1_TAD1       750001  1300000
#chr1    chr1_TAD2       2750001 3650000
#chr1    chr1_TAD3       3650001 4150000

gene2tadDT_file <- paste0(setDir, "/mnt/etemp/marie/Dixon2018_integrative_data/gene_data_final/GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460/genes2tad/all_genes_positions.txt")
#LINC00115       chr1    761586  762902  chr1_TAD1
#FAM41C  chr1    803451  812283  chr1_TAD1
#SAMD11  chr1    860260  879955  chr1_TAD1
#NOC2L   chr1    879584  894689  chr1_TAD1

# overwrite main_settings.R: nCpu <- 25
nCpu <- 20

# *************************************************************************************************************************
# ************************************ SETTINGS FOR PERMUTATIONS (5#_, 8c_)
# *************************************************************************************************************************

# number of permutations
nRandomPermut <- 10000
gene2tadAssignMethod <- "maxOverlap"
nRandomPermutShuffle <- 10000
step8_for_permutGenes <- TRUE
step8_for_randomTADsFix <- FALSE
step8_for_randomTADsGaussian <- FALSE
step8_for_randomTADsShuffle <- FALSE
step14_for_randomTADsShuffle <- FALSE

