#!/usr/bin/bash

# ./draw_matrix.sh



#Rscript draw_matrix.R \
#-m /mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_aggregMaps/AGGREG_MAPS_RESCALE_TOTCOUNT_NORMmat/chr21/chr21_aggregMap.txt \
#-o DRAW_METATAD/RESCALE_TOTCOUNT_NORMmat/chr21_draw_test_meanInterCount_hierarchLev25_35280001_40560000.png \
#-k 0 \
#-b 40000 \
#-s 35280001 \
#-e 40560000 \
#-c chr21 \
#-d 25

# !!! warning hard-coded settings in draw_matrix_metaTAD.R !!!
#- matrixHeader 
#- tadHeader
#- metaTADheader 
#- featureHeader 

maxJobs=40
maxLoad=70

script_name="draw_matrix.R"


binSizeKb=40
colToSkip=3


#all_chromos=( "chr1" "chr10" "chr21" )
#all_chromos=( "chr5" )
all_chromos=( "chr"{1..22} ) 



#matrixFolder="/mnt/etemp/marie/GSE87112/SB_GSM2322555/SAM_DanieleParam/danieleParam_20kb_matrix_correct"
#matrixPrefix="GSM2322555_danieleParam_${binSizeKb}kb_"
#matrixSuffix="_correct_dekker_matrix.txt"

#/mnt/etemp/marie/Dixon2018_integrative_data/breast/T47D/ENCSR549MGQ_GSE105697/GSE105697_ENCFF364CWZ_chromatin_interactions_hg19_chr5_TopDom.matrix

#matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/breast/T47D/ENCSR549MGQ_GSE105697"
#matrixPrefix="GSE105697_ENCFF364CWZ_chromatin_interactions_hg19_"
#matrixSuffix="_TopDom.matrix"
#outFolder="DRAW_ENCODE_MATRIX/GSE105697_ENCFF364CWZ_T47D"

#matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/breast/MCF7/GSM1631185_GSE66733/Hi-C_MCF7_MCF10A_processed_HiCfiles/Heatmaps/chrxchr/40kb/MCF7"
#matrixPrefix="HiCStein-MCF7-WT__hg19__"
#matrixSuffix="__C-40000-iced_TopDom.matrix"
#outFolder="DRAW_ENCODE_MATRIX/GSM1631185_MCF7"

#matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/breast/MCF7/GSM1942100_GSM1942101_GSE75070/shGFP"
#matrixPrefix="GSE75070_HiCStein-MCF7-shGFP_hg19_"
#matrixSuffix="_C-40000-iced_TopDom.matrix"
#outFolder="DRAW_ENCODE_MATRIX/GSE75070_MCF7_shGFP"

#matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/lung/A549/ENCSR444WCZ"
#matrixPrefix="GSE105600_ENCFF852YOE_chromatin_interactions_hg19_"
#matrixSuffix="_TopDom.matrix"
#outFolder="DRAW_ENCODE_MATRIX/GSE105600_ENCFF852YOE_A549"

matrixFolder="/mnt/etemp/marie/Dixon2018_integrative_data/lung/NCI-H460/ENCSR489OCU"
matrixPrefix="GSE105725_ENCFF697NNX_chromatin_interactions_hg19_"
matrixSuffix="_TopDom.matrix"
outFolder="DRAW_ENCODE_MATRIX/GSE105725_ENCFF697NNX_NCIH460"


parallel -i -j $maxJobs -l $maxLoad sh -c "echo Rscript $script_name  -m $matrixFolder/$matrixPrefix{}$matrixSuffix -o $outFolder/$matrixPrefix{}${matrixSuffix}.png -k $colToSkip -b ${binSizeKb}000 -c {}" -- ${all_chromos[@]}
parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript $script_name  -m $matrixFolder/$matrixPrefix{}$matrixSuffix -o $outFolder/$matrixPrefix{}${matrixSuffix}.png -k $colToSkip -b ${binSizeKb}000 -c {}" -- ${all_chromos[@]}


start_positions=( "35280001" )
end_positions=( "40560000" )

for i in "${!start_positions[@]}"; do
    start_pos="${start_positions[i]}"
    end_pos="${end_positions[i]}"

    parallel -i -j $maxJobs -l $maxLoad sh -c "echo Rscript $script_name  -m $matrixFolder/$matrixPrefix{}$matrixSuffix -o $outFolder/$matrixPrefix{}${matrixSuffix}_${start_pos}_${end_pos}.png -k $colToSkip -b ${binSizeKb}000 -c {} -s $start_pos -e $end_pos" -- ${all_chromos[@]}
    parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript $script_name  -m $matrixFolder/$matrixPrefix{}$matrixSuffix -o $outFolder/$matrixPrefix{}${matrixSuffix}_${start_pos}_${end_pos}.png -k $colToSkip -b ${binSizeKb}000 -c {} -s $start_pos -e $end_pos" -- ${all_chromos[@]}

done


