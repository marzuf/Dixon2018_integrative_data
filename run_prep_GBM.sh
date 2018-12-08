#!/bin/bash


# ./run_prep_GBM.sh

start_time=$(date -R)    


gbm_script="prep_GBM_data.R"


all_chromos=( "chr"{1..22} "chrX" )
#all_chromos=( "chr21" )
#all_chromos=( "chr9" )

maxJobs=40
maxLoad=70

parallel -i -j $maxJobs -l $maxLoad sh -c "echo Rscript $gbm_script glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176966_GB176_{}_40kb_norm_int_count.txt {} 40000 40000 glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176966_GB176_{}_40kb_norm_int_count_TopDom.matrix" -- ${all_chromos[@]}
parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript $gbm_script glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176966_GB176_{}_40kb_norm_int_count.txt {} 40000 40000 glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176966_GB176_{}_40kb_norm_int_count_TopDom.matrix" -- ${all_chromos[@]}

parallel -i -j $maxJobs -l $maxLoad sh -c "echo Rscript $gbm_script glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176967_GB180_{}_40kb_norm_int_count.txt {} 40000 40000 glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176967_GB180_{}_40kb_norm_int_count_TopDom.matrix" -- ${all_chromos[@]}
parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript $gbm_script glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176967_GB180_{}_40kb_norm_int_count.txt {} 40000 40000 glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176967_GB180_{}_40kb_norm_int_count_TopDom.matrix" -- ${all_chromos[@]}

parallel -i -j $maxJobs -l $maxLoad sh -c "echo Rscript $gbm_script glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176968_GB182_{}_40kb_norm_int_count.txt {} 40000 40000 glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176968_GB182_{}_40kb_norm_int_count_TopDom.matrix" -- ${all_chromos[@]}
parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript $gbm_script glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176968_GB182_{}_40kb_norm_int_count.txt {} 40000 40000 glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176968_GB182_{}_40kb_norm_int_count_TopDom.matrix" -- ${all_chromos[@]}

parallel -i -j $maxJobs -l $maxLoad sh -c "echo Rscript $gbm_script glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176969_GB183_{}_40kb_norm_int_count.txt {} 40000 40000 glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176969_GB183_{}_40kb_norm_int_count_TopDom.matrix" -- ${all_chromos[@]}
parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript $gbm_script glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176969_GB183_{}_40kb_norm_int_count.txt {} 40000 40000 glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176969_GB183_{}_40kb_norm_int_count_TopDom.matrix" -- ${all_chromos[@]}

parallel -i -j $maxJobs -l $maxLoad sh -c "echo Rscript $gbm_script glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176970_GB238_{}_40kb_norm_int_count.txt {} 40000 40000 glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176970_GB238_{}_40kb_norm_int_count_TopDom.matrix" -- ${all_chromos[@]}
parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript $gbm_script glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176970_GB238_{}_40kb_norm_int_count.txt {} 40000 40000 glioblastoma/GSE81879/COUNTS/GSE81879_GSM2176970_GB238_{}_40kb_norm_int_count_TopDom.matrix" -- ${all_chromos[@]}


###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################
echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0








