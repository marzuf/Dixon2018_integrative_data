#!/bin/bash


# ./cmds_TADcall.sh

start_time=$(date -R)    

#*********************** TO SET ***********************
dataID="$1"
#**********************************************

maxJobs=40
maxLoad=70


topdom_script="/mnt/ed2/shared/TADcompare/Software/TopDom/run_TopDom.R"
window_size=5


if [[ $dataID == "GSM1631185" ]]; then
# /mnt/etemp/marie/Dixon2018_integrative_data/breast/MCF7/GSM1631185_GSE66733/Hi-C_MCF7_MCF10A_processed_HiCfiles/Heatmaps/chrxchr/40kb/MCF7/HiCStein-MCF7-WT__hg19__chr1__C-40000-iced_TopDom.matrix
    inFold="/mnt/etemp/marie/Dixon2018_integrative_data/breast/MCF7/GSM1631185_GSE66733/Hi-C_MCF7_MCF10A_processed_HiCfiles/Heatmaps/chrxchr/40kb/MCF7"
    outFold="/mnt/etemp/marie/Dixon2018_integrative_data/breast/MCF7/GSM1631185_GSE66733/MCF7/TopDom"
    mkdir -p outFold
    inFilePrefix="HiCStein-MCF7-WT__hg19__"
    inFileSuffix="__C-40000-iced_TopDom.matrix"
    outFilePrefix="GSM1631185_MCF7_40kb_"

elif [[ $dataID == "GSE75070" ]]; then
# /mnt/etemp/marie/Dixon2018_integrative_data/breast/MCF7/GSM1942100_GSM1942101_GSE75070/shGFP/GSE75070_HiCStein-MCF7-shGFP_hg19_chr1_C-40000-iced_TopDom.matrix
    inFold="/mnt/etemp/marie/Dixon2018_integrative_data/breast/MCF7/GSM1942100_GSM1942101_GSE75070/shGFP/"
    outFold="/mnt/etemp/marie/Dixon2018_integrative_data/breast/MCF7/GSM1942100_GSM1942101_GSE75070/shGFP/TopDom"
    mkdir -p outFold
    inFilePrefix="GSE75070_HiCStein-MCF7-shGFP_hg19_"
    inFileSuffix="_C-40000-iced_TopDom.matrix"
    outFilePrefix="GSE75070_MCF7_shGFP_40kb_"

elif [[ $dataID == "GSE105697" ]]; then
#/mnt/etemp/marie/Dixon2018_integrative_data/breast/T47D/ENCSR549MGQ_GSE105697/GSE105697_ENCFF364CWZ_chromatin_interactions_hg19_chr1_TopDom.matrix
    inFold="/mnt/etemp/marie/Dixon2018_integrative_data/breast/T47D/ENCSR549MGQ_GSE105697"
    outFold="/mnt/etemp/marie/Dixon2018_integrative_data/breast/T47D/ENCSR549MGQ_GSE105697/TopDom"
    mkdir -p outFold
    inFilePrefix="GSE105697_ENCFF364CWZ_chromatin_interactions_hg19_"
    inFileSuffix="_TopDom.matrix"
    outFilePrefix="GSE105697_ENCFF364CWZ_T47D_40kb_"

elif [[ $dataID == "GSE105600" ]]; then
#/mnt/etemp/marie/Dixon2018_integrative_data/lung/A549/ENCSR444WCZ/GSE105600_ENCFF852YOE_chromatin_interactions_hg19_chr1_TopDom.matrix
    inFold="/mnt/etemp/marie/Dixon2018_integrative_data/lung/A549/ENCSR444WCZ"
    outFold="/mnt/etemp/marie/Dixon2018_integrative_data/lung/A549/ENCSR444WCZ/TopDom"
    mkdir -p outFold
    inFilePrefix="GSE105600_ENCFF852YOE_chromatin_interactions_hg19_"
    inFileSuffix="_TopDom.matrix"
    outFilePrefix="GSE105600_ENCFF852YOE_A549_40kb_"


elif [[ $dataID == "GSE105725" ]]; then
#/mnt/etemp/marie/Dixon2018_integrative_data/lung/NCI-H460/ENCSR489OCU/GSE105725_ENCFF697NNX_chromatin_interactions_hg19_chr2_TopDom.matrix
    inFold="/mnt/etemp/marie/Dixon2018_integrative_data/lung/NCI-H460/ENCSR489OCU"
    outFold="/mnt/etemp/marie/Dixon2018_integrative_data/lung/NCI-H460/ENCSR489OCU/TopDom"
    mkdir -p outFold
    inFilePrefix="GSE105725_ENCFF697NNX_chromatin_interactions_hg19_"
    inFileSuffix="_TopDom.matrix"
    outFilePrefix="GSE105725_ENCFF697NNX_NCIH460_40kb_"


elif [[ $dataID == "GSE105566" ]]; then
#/mnt/etemp/marie/Dixon2018_integrative_data/pancreas/Panc1/GSE105566/ENCFF358MNA/GSE105566_ENCFF358MNA_chromatin_interactions_hg19_chr21_TopDom.matrix
    inFold="/mnt/etemp/marie/Dixon2018_integrative_data/pancreas/Panc1/GSE105566/ENCFF358MNA"
    outFold="/mnt/etemp/marie/Dixon2018_integrative_data/pancreas/Panc1/GSE105566/ENCFF358MNA/TopDom"
    mkdir -p outFold
    inFilePrefix="GSE105566_ENCFF358MNA_chromatin_interactions_hg19_"
    inFileSuffix="_TopDom.matrix"
    outFilePrefix="GSE105566_ENCFF358MNA_Panc1_40kb_"

elif [[ $dataID == "GSE105557" ]]; then
#/mnt/etemp/marie/Dixon2018_integrative_data/prostate/LNCaP/GSE105557/ENCFF270HJX/GSE105557_ENCFF270HJX_chromatin_interactions_hg19_chr21_TopDom.matrix
    inFold="/mnt/etemp/marie/Dixon2018_integrative_data/prostate/LNCaP/GSE105557/ENCFF270HJX"
    outFold="/mnt/etemp/marie/Dixon2018_integrative_data/prostate/LNCaP/GSE105557/ENCFF270HJX/TopDom"
    mkdir -p outFold
    inFilePrefix="GSE105557_ENCFF270HJX_chromatin_interactions_hg19_"
    inFileSuffix="_TopDom.matrix"
    outFilePrefix="GSE105557_ENCFF270HJX_LNCaP_40kb_"

elif [[ $dataID == "GSE105465" ]]; then
#/mnt/etemp/marie/Dixon2018_integrative_data/kidney/Caki2/GSE105465/ENCFF777DUA/GSE105465_ENCFF777DUA_chromatin_interactions_hg19_chr2_TopDom.matrix
    inFold="/mnt/etemp/marie/Dixon2018_integrative_data/kidney/Caki2/GSE105465/ENCFF777DUA"
    outFold="/mnt/etemp/marie/Dixon2018_integrative_data/kidney/Caki2/GSE105465/ENCFF777DUA/TopDom"
    mkdir -p outFold
    inFilePrefix="GSE105465_ENCFF777DUA_chromatin_interactions_hg19_"
    inFileSuffix="_TopDom.matrix"
    outFilePrefix="GSE105465_ENCFF777DUA_Caki2_40kb_"

elif [[ $dataID == "GSE105235" ]]; then
#/mnt/etemp/marie/Dixon2018_integrative_data/kidney/G401/GSE105235/ENCFF235TGH/GSE105235_ENCFF235TGH_chromatin_interactions_hg19_chr9_TopDom.matrix
    inFold="/mnt/etemp/marie/Dixon2018_integrative_data/kidney/G401/GSE105235/ENCFF235TGH"
    outFold="/mnt/etemp/marie/Dixon2018_integrative_data/kidney/G401/GSE105235/ENCFF235TGH/TopDom"
    mkdir -p outFold
    inFilePrefix="GSE105235_ENCFF235TGH_chromatin_interactions_hg19_"
    inFileSuffix="_TopDom.matrix"
    outFilePrefix="GSE105235_ENCFF235TGH_G401_40kb_"

elif [[ $dataID == "GSE106022" ]]; then
#/mnt/etemp/marie/Dixon2018_integrative_data/skin/RPMI-7951/GSE106022/ENCFF614EKT/GSE106022_ENCFF614EKT_chromatin_interactions_hg19_chr21_TopDom.matrix
    inFold="/mnt/etemp/marie/Dixon2018_integrative_data/skin/RPMI-7951/GSE106022/ENCFF614EKT"
    outFold="/mnt/etemp/marie/Dixon2018_integrative_data/skin/RPMI-7951/GSE106022/ENCFF614EKT/TopDom"
    mkdir -p outFold
    inFilePrefix="GSE106022_ENCFF614EKT_chromatin_interactions_hg19_"
    inFileSuffix="_TopDom.matrix"
    outFilePrefix="GSE106022_ENCFF614EKT_RPMI7951_40kb_"

elif [[ $dataID == "GSE105491" ]]; then
#/mnt/etemp/marie/Dixon2018_integrative_data/skin/SK-MEL-5/GSE105491/ENCFF458OWO/GSE105491_ENCFF458OWO_chromatin_interactions_hg19_chr21_TopDom.matrix
    inFold="/mnt/etemp/marie/Dixon2018_integrative_data/skin/SK-MEL-5/GSE105491/ENCFF458OWO"
    outFold="/mnt/etemp/marie/Dixon2018_integrative_data/skin/SK-MEL-5/GSE105491/ENCFF458OWO/TopDom"
    mkdir -p outFold
    inFilePrefix="GSE105491_ENCFF458OWO_chromatin_interactions_hg19_"
    inFileSuffix="_TopDom.matrix"
    outFilePrefix="GSE105491_ENCFF458OWO_SKMEL5_40kb_"


fi

echo "... inFold = $inFold"
echo "... outFold = $outFold"
echo "... inFilePrefix = $inFilePrefix"
echo "... inFileSuffix = $inFileSuffix"
echo "... outFilePrefix = $outFilePrefix"



#all_chromos=( "chr"{1..22} "chrX" )
all_chromos=( "chr1" )
#all_chromos=( "chr21" )


parallel -i -j $maxJobs -l $maxLoad sh -c "echo Rscript $topdom_script -i $inFold/${inFilePrefix}{}${inFileSuffix} -o $outFold/${outFilePrefix}{} -w $window_size" -- ${all_chromos[@]}
parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript $topdom_script -i $inFold/${inFilePrefix}{}${inFileSuffix} -o $outFold/${outFilePrefix}{} -w $window_size" -- ${all_chromos[@]}


###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################
echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0


