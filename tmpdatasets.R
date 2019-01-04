topdomPattern <- "_final_domains.txt$"
consensusPattern <- "_conservedTADs.txt$"
setDir=""
#================================================
#================================================ BREAST
#================================================

breastConsensusname <-  "GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP_vs_GSE105697_ENCFF364CWZ_T47D"
breastConsensusFold <- file.path("FIND_CONSENSUS_TADS", breastConsensusname)
breastConsensusFiles <- list.files(breastConsensusFold, full.names=T, pattern = consensusPattern)

cat("breastConsensusname:\t", breastConsensusname, "\n")
cat("breastConsensusFiles[1]:\t", breastConsensusFiles[1], "\n")


mcf7Consensusname <-  "GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP"
mcf7ConsensusFold <- file.path("FIND_CONSENSUS_TADS",mcf7Consensusname)
mcf7ConsensusFiles <- list.files(mcf7ConsensusFold, full.names=T, pattern = consensusPattern)

cat("mcf7Consensusname:\t", mcf7Consensusname, "\n")
cat("mcf7ConsensusFiles[1]:\t", mcf7ConsensusFiles[1], "\n")


# breast/MCF7/GSM1631185_GSE66733/MCF7/TopDom/GSM1631185_MCF7_40kb_chr1_final_domains.txt
breastCL1name <- "GSM1631185_GSE66733"
breastFold1 <- file.path("breast/MCF7/", breastCL1name, "MCF7/TopDom")
breastCL1Files <- list.files(breastFold1, full.names=T, pattern = topdomPattern)

cat("breastCL1name:\t", breastCL1name, "\n")
cat("breastCL1Files[1]:\t", breastCL1Files[1], "\n")


# breast/MCF7/GSM1942100_GSM1942101_GSE75070/shGFP/TopDom/GSE75070_MCF7_shGFP_40kb_chr1_final_domains.txt
breastCL2name <- "GSM1942100_GSM1942101_GSE75070"
breastFold2 <- file.path("breast/MCF7", breastCL2name, "shGFP/TopDom")
breastCL2Files <- list.files(breastFold2, full.names=T, pattern = topdomPattern)

cat("breastCL2name:\t", breastCL2name, "\n")
cat("breastCL2Files[1]:\t", breastCL2Files[1], "\n")


# breast/T47D/ENCSR549MGQ_GSE105697/TopDom/GSE105697_ENCFF364CWZ_T47D_40kb_chr1_final_domains.txt
breastCL3name <- "ENCSR549MGQ_GSE105697"
breastFold3 <- file.path("breast/T47D", breastCL3name, "TopDom")
breastCL3Files <- list.files(breastFold3, full.names=T, pattern = topdomPattern)

cat("breastCL3name:\t", breastCL3name, "\n")
cat("breastCL3Files[1]:\t", breastCL3Files[1], "\n")


#================================================
#================================================ LUNG
#================================================

lungConsensusname <-  "GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460"
lungConsensusFold <- file.path("FIND_CONSENSUS_TADS", lungConsensusname)
lungConsensusFiles <- list.files(lungConsensusFold, full.names=T, pattern = consensusPattern)

cat("lungConsensusname:\t", lungConsensusname, "\n")
cat("lungConsensusFiles[1]:\t", lungConsensusFiles[1], "\n")


# lung/A549/ENCSR444WCZ/TopDom/GSE105600_ENCFF852YOE_A549_40kb_chr1_final_domains.txt
lungCL1name <- "ENCSR444WCZ"
lungFold1 <- file.path("lung/A549", lungCL1name, "TopDom")
lungCL1Files <- list.files(lungFold1, full.names=T, pattern = topdomPattern)

cat("lungCL1name:\t", lungCL1name, "\n")
cat("lungCL1Files[1]:\t", lungCL1Files[1], "\n")


# lung/NCI-H460/ENCSR489OCU/TopDom/GSE105725_ENCFF697NNX_NCIH460_40kb_chr1_final_domains.txt
lungCL2name <- "ENCSR489OCU"
lungFold2 <- file.path("lung/NCI-H460", lungCL2name, "TopDom")
lungCL2Files <- list.files(lungFold2, full.names=T, pattern = topdomPattern)

cat("lungCL2name:\t", lungCL2name, "\n")
cat("lungCL2Files[1]:\t", lungCL2Files[1], "\n")


#================================================
#================================================ PANCREAS
#================================================


#pancreas/Panc1/GSE105566/ENCFF358MNA/TopDom GSE105566_ENCFF358MNA_Panc1_40kb_chr.+final_domains.txt 
pancreasCL1name <- "ENCFF358MNA"
pancreasFold1 <- file.path("pancreas/Panc1/GSE105566", pancreasCL1name, "TopDom")
pancreasCL1Files <- list.files(pancreasFold1, full.names=T, pattern = topdomPattern)

cat("pancreasCL1name:\t", pancreasCL1name, "\n")
cat("pancreasCL1Files[1]:\t", pancreasCL1Files[1], "\n")


#prostate/LNCaP/GSE105557/ENCFF270HJX/TopDom GSE105557_ENCFF270HJX_LNCaP_40kb_chr.+final_domains.txt 
prostateCL1name <- "ENCFF270HJX"
prostateFold1 <- file.path("prostate/LNCaP/GSE105557", prostateCL1name, "TopDom")
prostateCL1Files <- list.files(prostateFold1, full.names=T, pattern = topdomPattern)

cat("prostateCL1name:\t", prostateCL1name, "\n")
cat("prostateCL1Files[1]:\t", prostateCL1Files[1], "\n")


#================================================
#================================================ KIDNEY
#================================================

#kidney/Caki2/GSE105465/ENCFF777DUA/TopDom GSE105465_ENCFF777DUA_Caki2_40kb_chr.+final_domains.txt 
kidneyCL1name <- "ENCFF777DUA"
kidneyFold1 <- file.path("kidney/Caki2/GSE105465", kidneyCL1name, "TopDom")
kidneyCL1Files <- list.files(kidneyFold1, full.names=T, pattern = topdomPattern)

cat("kidneyCL1name:\t", kidneyCL1name, "\n")
cat("kidneyCL1Files[1]:\t", kidneyCL1Files[1], "\n")


#kidney/G401/GSE105235/ENCFF235TGH/TopDom GSE105235_ENCFF235TGH_G401_40kb_chr.+final_domains.txt 
kidneyCL2name <- "ENCFF235TGH"
kidneyFold2 <- file.path("kidney/G401/GSE105235", kidneyCL2name, "TopDom")
kidneyCL2Files <- list.files(kidneyFold2, full.names=T, pattern = topdomPattern)

cat("kidneyCL2name:\t", kidneyCL2name, "\n")
cat("kidneyCL2Files[1]:\t", kidneyCL2Files[1], "\n")


kidneyConsensusname <-  "GSE105465_ENCFF777DUA_Caki2_vs_GSE105235_ENCFF235TGH_G401"
kidneyConsensusFold <- file.path("FIND_CONSENSUS_TADS", kidneyConsensusname)
kidneyConsensusFiles <- list.files(kidneyConsensusFold, full.names=T, pattern = consensusPattern)

cat("kidneyConsensusname:\t", kidneyConsensusname, "\n")
cat("kidneyConsensusFiles[1]:\t", kidneyConsensusFiles[1], "\n")



#================================================
#================================================ SKIN
#================================================


#skin/RPMI-7951/GSE106022/ENCFF614EKT/TopDom GSE106022_ENCFF614EKT_RPMI7951_40kb_chr.+final_domains.txt 
skinCL1name <- "ENCFF614EKT"
skinFold1 <- file.path("skin/RPMI-7951/GSE106022", skinCL1name, "TopDom")
skinCL1Files <- list.files(skinFold1, full.names=T, pattern = topdomPattern)

cat("skinCL1name:\t", skinCL1name, "\n")
cat("skinCL1Files[1]:\t", skinCL1Files[1], "\n")


#skin/SK-MEL-5/GSE105491/ENCFF458OWO/TopDom GSE105491_ENCFF458OWO_SKMEL5_40kb_chr.+final_domains.txt 
skinCL2name <- "ENCFF458OWO"
skinFold2 <- file.path("skin/SK-MEL-5/GSE105491", skinCL2name, "TopDom")
skinCL2Files <- list.files(skinFold2, full.names=T, pattern = topdomPattern)

cat("skinCL2name:\t", skinCL2name, "\n")
cat("skinCL2Files[1]:\t", skinCL2Files[1], "\n")


skinConsensusname <-  "GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5"
skinConsensusFold <- file.path("FIND_CONSENSUS_TADS", skinConsensusname)
skinConsensusFiles <- list.files(skinConsensusFold, full.names=T, pattern = consensusPattern)

cat("skinConsensusname:\t", skinConsensusname, "\n")
cat("skinConsensusFiles[1]:\t", skinConsensusFiles[1], "\n")



#================================================
#================================================ COLORECTAL
#================================================


#/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF439QFU/GSE105318_ENCFF439QFU_chromatin_interactions_hg19_chr1_TopDom.matrix → float
#/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF714TMN/GSE105318_ENCFF714TMN_chromatin_interactions_hg19_chr1_TopDom.matrix → int ???
#"coloCL1" = "GSE105318_ENCFF439QFU",
#"coloCL2" = "GSE105318_ENCFF714TMN", # int



coloCL1name <- "ENCFF439QFU"
coloFold1 <- file.path("colon/DLD1/GSE105318", coloCL1name, "TopDom")
coloCL1Files <- list.files(coloFold1, full.names=T, pattern = topdomPattern)
stopifnot(length(coloCL1Files) > 0)
colo1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(coloCL1Files)))
stopifnot(length(colo1_chromos) > 0)

cat("coloCL1name:\t", coloCL1name, "\n")
cat("coloCL1Files[1]:\t", coloCL1Files[1], "\n")


## => int matrix
#coloCL2name <- "ENCFF714TMN"
#coloFold2 <- file.path("colon/DLD1/GSE105318", coloCL2name, "TopDom")
#coloCL2Files <- list.files(coloFold2, full.names=T, pattern = topdomPattern)
#stopifnot(length(coloCL2Files) > 0)
#colo2_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(coloCL2Files)))
#stopifnot(length(colo2_chromos) > 0)


                    # NOT DONE
                    #coloConsensusname <-  "GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5"
                    #coloConsensusFold <- file.path("FIND_CONSENSUS_TADS", coloConsensusname)
                    #coloConsensusFiles <- list.files(coloConsensusFold, full.names=T, pattern = consensusPattern)
                    #stopifnot(length(coloConsensusFiles) > 0)
                    #colo_consensus_chromos <- unique(gsub("(chr.+)_conservedTADs.txt", "\\1", basename(coloConsensusFiles)))
                    #stopifnot(length(colo_consensus_chromos) > 0)

#================================================
#================================================ ASTROCYTES
#================================================


#astrocyte/cerebellum/GSE105194/ENCFF027IEO/GSE105194_ENCFF027IEO_chromatin_interactions_hg19_chr1_TopDom.matrix → float
#astrocyte/cerebellum/GSE105194/ENCFF122YID/GSE105194_ENCFF122YID_chromatin_interactions_hg19_chr1_TopDom.matrix → int ???
#astrocyte/spinal_cord/GSE105957/ENCFF478UBU/GSE105957_ENCFF478UBU_chromatin_interactions_hg19_chr1_TopDom.matrix → int ???
#astrocyte/spinal_cord/GSE105957/ENCFF715HDW/GSE105957_ENCFF715HDW_chromatin_interactions_hg19_chr1_TopDom.matrix → float

# "astroCL1" = "GSE105194_ENCFF027IEO",
# "astroCL2"= "GSE105194_ENCFF122YID", # int
# "astroCL3"= "GSE105957_ENCFF715HDW",
#"astroCL4"= "GSE105957_ENCFF478UBU", # int

astroCL1name <- "ENCFF027IEO"
astroFold1 <- file.path("astrocyte/cerebellum/GSE105194", astroCL1name, "TopDom")
astroCL1Files <- list.files(astroFold1, full.names=T, pattern = topdomPattern)

cat("astroCL1name:\t", astroCL1name, "\n")
cat("astroCL1Files[1]:\t", astroCL1Files[1], "\n")


# int
#astroCL2name <- "ENCFF122YID"
#astroFold2 <- file.path("astrocyte/cerebellum/GSE105194/ENCFF027IEO", astroCL2name, "TopDom")
#astroCL2Files <- list.files(astroFold2, full.names=T, pattern = topdomPattern)
#stopifnot(length(astroCL2Files) > 0)
#astro1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(astroCL2Files)))
#stopifnot(length(astro2_chromos) > 0)


astroCL3name <- "ENCFF715HDW"
astroFold3 <- file.path("astrocyte/spinal_cord/GSE105957", astroCL3name, "TopDom")
astroCL3Files <- list.files(astroFold3, full.names=T, pattern = topdomPattern)

cat("astroCL3name:\t", astroCL3name, "\n")
cat("astroCL3Files[1]:\t", astroCL3Files[1], "\n")


# int
#astroCL4name <- "ENCFF478UBU"
#astroFold4 <- file.path("astrocyte/spinal_cord/GSE105957", astroCL4name, "TopDom")
#astroCL4Files <- list.files(astroFold4, full.names=T, pattern = topdomPattern)
#stopifnot(length(astroCL4Files) > 0)
#astro1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(astroCL4Files)))
#stopifnot(length(astro4_chromos) > 0)


# based only on astrocL1 and astroCL3 (float matrices) FIND_CONSENSUS_TADS/GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal/chrX_conservedTADs.txt
# NOT DONE
astroConsensusname <-  "GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal"
astroConsensusFold <- file.path("FIND_CONSENSUS_TADS", astroConsensusname)
astroConsensusFiles <- list.files(astroConsensusFold, full.names=T, pattern = consensusPattern)

cat("astroConsensusname:\t", astroConsensusname, "\n")
cat("astroConsensusFiles[1]:\t", astroConsensusFiles[1], "\n")



#================================================
#================================================ LIVER
#================================================



#================================================
#================================================ LYMPHOBlAST k562 - LEUKEMIAS
#================================================

#/mnt/etemp/marie/Dixon2018_integrative_data/leukemia/K562/GSE63525/GSE63525_K562_40kb_ICE_chr1_TopDom.matrix
#"lympho1" = "GSE63525_K562"

#lymphoCL1name <- "GSE63525_K562"
#lymphoFold1 <- file.path("leukemia/K562/GSE63525", lymphoCL1name, "TopDom")
#lymphoCL1Files <- list.files(lymphoFold1, full.names=T, pattern = topdomPattern)
#stopifnot(length(lymphoCL1Files) > 0)
#lympho1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(lymphoCL1Files)))
#stopifnot(length(lympho1_chromos) > 0)


#================================================
#================================================ PIPELINE
#================================================


pipConsensusname <- "pipeline_TopDom"
pipConsensusFold <- file.path(setDir, "/mnt/ed4/marie/TAD_call_pipeline_TopDom", "consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final")
pipConsensusFiles <- list.files(pipConsensusFold, full.names=T, pattern = consensusPattern)

cat("pipConsensusname:\t", pipConsensusname, "\n")
cat("pipConsensusFiles[1]:\t", pipConsensusFiles[1], "\n")

#================================================
#================================================ NAME SETTINGS
#================================================


ds_mapping <- c(
"breastCL1" = "HiCStein-MCF7-WT__hg19__.+TopDom.matrix" , 
"breastCL2" = "GSE75070_HiCStein-MCF7-shGFP_hg19_.+TopDom.matrix",
"breastCL3" = "GSE105697_ENCFF364CWZ", 

"lungCL1" = "GSE105600_ENCFF852YOE",
"lungCL2" = "GSE105725_ENCFF697NNX",

"pancreasCL1" = "GSE105566_ENCFF358MNA", 

"prostateCL1" = "GSE105557_ENCFF270HJX",

"kidneyCL1" = "GSE105465_ENCFF777DUA", 
"kidneyCL2" =  "GSE105235_ENCFF235TGH",

"skinCL1" = "GSE106022_ENCFF614EKT", 
"skinCL2" = "GSE105491_ENCFF458OWO",


"coloCL1" = "GSE105318_ENCFF439QFU",
"coloCL2" = "GSE105318_ENCFF714TMN", # int

"astroCL1" = "GSE105194_ENCFF027IEO",
 "astroCL2"= "GSE105194_ENCFF122YID", # int
"astroCL3"= "GSE105957_ENCFF715HDW",
"astroCL4"= "GSE105957_ENCFF478UBU", # int

"lympho1" = "GSE63525_K562_40kb_ICE"


)

cat("\n")
for(i in seq_along(ds_mapping)) {
  cat("> ", i, "- ", names(ds_mapping)[i], "\t=\t", ds_mapping[i], "\n")
}
