#!/bin/bash


start_time=$(date -R)    


convert_script="/mnt/etemp/marie/Dixon2018_integrative_data/convert_TopDom_matrix.R"

#matrixDir="/mnt/etemp/marie/Dixon2018_integrative_data/breast/MCF7/GSM1631185/Hi-C_MCF7_MCF10A_processed_HiCfiles/Heatmaps/chrxchr/40kb/MCF7"
#filePrefix="$matrixDir/HiCStein-MCF7-WT__hg19__"
#fileSuffix="__C-40000-iced.matrix"
#fileOutSuffix="__C-40000-iced_TopDom.matrix"


matrixDir="/mnt/etemp/marie/Dixon2018_integrative_data/breast/MCF7/GSM1942100/shGFP"
filePrefix="$matrixDir/GSE75070_HiCStein-MCF7-shGFP_hg19_"
fileSuffix="_C-40000-iced.matrix"
fileOutSuffix="_C-40000-iced_TopDom.matrix"

all_chromos=( "chr"{1..22} "chrX" )
#all_chromos=( "chr1" )

maxJobs=40
maxLoad=70

parallel -i -j $maxJobs -l $maxLoad sh -c "echo Rscript $convert_script ${filePrefix}{}${fileSuffix} ${filePrefix}{}${fileOutSuffix}" -- ${all_chromos[@]}
parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript $convert_script ${filePrefix}{}${fileSuffix} ${filePrefix}{}${fileOutSuffix}" -- ${all_chromos[@]}



###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################
echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0


