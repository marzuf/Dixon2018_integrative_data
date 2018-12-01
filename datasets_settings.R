breastConsensusname <-  "GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP_vs_GSE105697_ENCFF364CWZ_T47D"
breastConsensusFold <- file.path("FIND_CONSENSUS_TADS", breastConsensusname)
breastConsensusFiles <- list.files(breastConsensusFold, full.names=T, pattern = consensusPattern)
stopifnot(length(breastConsensusFiles) > 0)
breast_consensus_chromos <- unique(gsub("(chr.+)_conservedTADs.txt", "\\1", basename(breastConsensusFiles)))
stopifnot(length(breast_consensus_chromos) > 0)

mcf7Consensusname <-  "GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP"
mcf7ConsensusFold <- file.path("FIND_CONSENSUS_TADS",mcf7Consensusname)
mcf7ConsensusFiles <- list.files(mcf7ConsensusFold, full.names=T, pattern = consensusPattern)
stopifnot(length(mcf7ConsensusFiles) > 0)
mcf7_consensus_chromos <- unique(gsub("(chr.+)_conservedTADs.txt", "\\1", basename(mcf7ConsensusFiles)))
stopifnot(length(mcf7_consensus_chromos) > 0)

lungConsensusname <-  "GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460"
lungConsensusFold <- file.path("FIND_CONSENSUS_TADS", lungConsensusname)
lungConsensusFiles <- list.files(lungConsensusFold, full.names=T, pattern = consensusPattern)
stopifnot(length(lungConsensusFiles) > 0)
lung_consensus_chromos <- unique(gsub("(chr.+)_conservedTADs.txt", "\\1", basename(lungConsensusFiles)))
stopifnot(length(lung_consensus_chromos) > 0)

skinConsensusname <-  "GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5"
skinConsensusFold <- file.path("FIND_CONSENSUS_TADS", skinConsensusname)
skinConsensusFiles <- list.files(skinConsensusFold, full.names=T, pattern = consensusPattern)
stopifnot(length(skinConsensusFiles) > 0)
skin_consensus_chromos <- unique(gsub("(chr.+)_conservedTADs.txt", "\\1", basename(skinConsensusFiles)))
stopifnot(length(skin_consensus_chromos) > 0)

kidneyConsensusname <-  "GSE105465_ENCFF777DUA_Caki2_vs_GSE105235_ENCFF235TGH_G401"
kidneyConsensusFold <- file.path("FIND_CONSENSUS_TADS", kidneyConsensusname)
kidneyConsensusFiles <- list.files(kidneyConsensusFold, full.names=T, pattern = consensusPattern)
stopifnot(length(kidneyConsensusFiles) > 0)
kidney_consensus_chromos <- unique(gsub("(chr.+)_conservedTADs.txt", "\\1", basename(kidneyConsensusFiles)))
stopifnot(length(kidney_consensus_chromos) > 0)



topdomPattern <- "_final_domains.txt$"

# breast/MCF7/GSM1631185_GSE66733/MCF7/TopDom/GSM1631185_MCF7_40kb_chr1_final_domains.txt
breastCL1name <- "GSM1631185_GSE66733"
breastFold1 <- file.path("breast/MCF7/", breastCL1name, "MCF7/TopDom")
breastCL1Files <- list.files(breastFold1, full.names=T, pattern = topdomPattern)
stopifnot(length(breastCL1Files) > 0)
breast1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(breastCL1Files)))
stopifnot(length(breast1_chromos) > 0)

# breast/MCF7/GSM1942100_GSM1942101_GSE75070/shGFP/TopDom/GSE75070_MCF7_shGFP_40kb_chr1_final_domains.txt
breastCL2name <- "GSM1942100_GSM1942101_GSE75070"
breastFold2 <- file.path("breast/MCF7", breastCL2name, "shGFP/TopDom")
breastCL2Files <- list.files(breastFold2, full.names=T, pattern = topdomPattern)
stopifnot(length(breastCL2Files) > 0)
breast2_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(breastCL2Files)))
stopifnot(length(breast2_chromos) > 0)

# breast/T47D/ENCSR549MGQ_GSE105697/TopDom/GSE105697_ENCFF364CWZ_T47D_40kb_chr1_final_domains.txt
breastCL3name <- "ENCSR549MGQ_GSE105697"
breastFold3 <- file.path("breast/T47D", breastCL3name, "TopDom")
breastCL3Files <- list.files(breastFold3, full.names=T, pattern = topdomPattern)
stopifnot(length(breastCL3Files) > 0)
breast3_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(breastCL3Files)))
stopifnot(length(breast3_chromos) > 0)

# lung/A549/ENCSR444WCZ/TopDom/GSE105600_ENCFF852YOE_A549_40kb_chr1_final_domains.txt
lungCL1name <- "ENCSR444WCZ"
lungFold1 <- file.path("lung/A549", lungCL1name, "TopDom")
lungCL1Files <- list.files(lungFold1, full.names=T, pattern = topdomPattern)
stopifnot(length(lungCL1Files) > 0)
lung1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(lungCL1Files)))
stopifnot(length(lung1_chromos) > 0)

# lung/NCI-H460/ENCSR489OCU/TopDom/GSE105725_ENCFF697NNX_NCIH460_40kb_chr1_final_domains.txt
lungCL2name <- "ENCSR489OCU"
lungFold2 <- file.path("lung/NCI-H460", lungCL2name, "TopDom")
lungCL2Files <- list.files(lungFold2, full.names=T, pattern = topdomPattern)
stopifnot(length(lungCL2Files) > 0)
lung2_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(lungCL2Files)))
stopifnot(length(lung2_chromos) > 0)



#pancreas/Panc1/GSE105566/ENCFF358MNA/TopDom GSE105566_ENCFF358MNA_Panc1_40kb_chr.+final_domains.txt 
pancreasCL1name <- "ENCFF358MNA"
pancreasFold1 <- file.path("pancreas/Panc1/GSE105566", pancreasCL1name, "TopDom")
pancreasCL1Files <- list.files(pancreasFold1, full.names=T, pattern = topdomPattern)
stopifnot(length(pancreasCL1Files) > 0)
pancreas1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(pancreasCL1Files)))
stopifnot(length(pancreas1_chromos) > 0)


#prostate/LNCaP/GSE105557/ENCFF270HJX/TopDom GSE105557_ENCFF270HJX_LNCaP_40kb_chr.+final_domains.txt 
prostateCL1name <- "ENCFF270HJX"
prostateFold1 <- file.path("prostate/LNCaP/GSE105557", prostateCL1name, "TopDom")
prostateCL1Files <- list.files(prostateFold1, full.names=T, pattern = topdomPattern)
stopifnot(length(prostateCL1Files) > 0)
prostate1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(prostateCL1Files)))
stopifnot(length(prostate1_chromos) > 0)

#kidney/Caki2/GSE105465/ENCFF777DUA/TopDom GSE105465_ENCFF777DUA_Caki2_40kb_chr.+final_domains.txt 
kidneyCL1name <- "ENCFF777DUA"
kidneyFold1 <- file.path("kidney/Caki2/GSE105465", kidneyCL1name, "TopDom")
kidneyCL1Files <- list.files(kidneyFold1, full.names=T, pattern = topdomPattern)
stopifnot(length(kidneyCL1Files) > 0)
kidney1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(kidneyCL1Files)))
stopifnot(length(kidney1_chromos) > 0)

#kidney/G401/GSE105235/ENCFF235TGH/TopDom GSE105235_ENCFF235TGH_G401_40kb_chr.+final_domains.txt 
kidneyCL2name <- "ENCFF235TGH"
kidneyFold2 <- file.path("kidney/G401/GSE105235", kidneyCL2name, "TopDom")
kidneyCL2Files <- list.files(kidneyFold2, full.names=T, pattern = topdomPattern)
stopifnot(length(kidneyCL2Files) > 0)
kidney2_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(kidneyCL2Files)))
stopifnot(length(kidney2_chromos) > 0)

#skin/RPMI-7951/GSE106022/ENCFF614EKT/TopDom GSE106022_ENCFF614EKT_RPMI7951_40kb_chr.+final_domains.txt 
skinCL1name <- "ENCFF614EKT"
skinFold1 <- file.path("skin/RPMI-7951/GSE106022", skinCL1name, "TopDom")
skinCL1Files <- list.files(skinFold1, full.names=T, pattern = topdomPattern)
stopifnot(length(skinCL1Files) > 0)
skin1_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(skinCL1Files)))
stopifnot(length(skin1_chromos) > 0)

#skin/SK-MEL-5/GSE105491/ENCFF458OWO/TopDom GSE105491_ENCFF458OWO_SKMEL5_40kb_chr.+final_domains.txt 
skinCL2name <- "ENCFF458OWO"
skinFold2 <- file.path("skin/SK-MEL-5/GSE105491", skinCL2name, "TopDom")
skinCL2Files <- list.files(skinFold2, full.names=T, pattern = topdomPattern)
stopifnot(length(skinCL2Files) > 0)
skin2_chromos <- unique(gsub(".+(chr.+)_final_domains.txt$", "\\1", basename(skinCL2Files)))
stopifnot(length(skin2_chromos) > 0)

#stop("--ok--\n")


pipConsensusname <- "pipeline_TopDom"
pipConsensusFold <- file.path(setDir, "/mnt/ed4/marie/TAD_call_pipeline_TopDom", "consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final")
pipConsensusFiles <- list.files(pipConsensusFold, full.names=T, pattern = consensusPattern)
pip_consensus_chromos <- unique(gsub("(chr.+)_conservedTADs.txt", "\\1", basename(pipConsensusFiles)))

intersectChromos <- Reduce(intersect, list(
  pip_consensus_chromos,
  breast_consensus_chromos, mcf7_consensus_chromos, 
  lung_consensus_chromos,
  breast1_chromos, breast2_chromos, breast3_chromos,
  lung1_chromos, lung2_chromos,
  pancreas1_chromos,
  prostate1_chromos,
  kidney1_chromos, kidney2_chromos,
  skin1_chromos, skin2_chromos
))

all_ds <- c(
  "pipConsensus",
  "breastConsensus", "mcf7Consensus", "lungConsensus", "kidneyConsensus", "skinConsensus",
  "breastCL1", "breastCL2",
  "lungCL1", "lungCL2",
  "pancreasCL1", 
  "prostateCL1",
  "kidneyCL1", "kidneyCL2",
  "skinCL1", "skinCL2"
)

