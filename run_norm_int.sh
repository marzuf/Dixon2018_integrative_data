#!/bin/bash


# ./run_norm_int.sh

start_time=$(date -R)    


norm_script="norm_int_data.R"


all_chromos=( "chr"{1..22} "chrX" )
#all_chromos=( "chr21" )
#all_chromos=( "chr9" )

maxJobs=40
maxLoad=70

parallel -i -j $maxJobs -l $maxLoad sh -c "echo Rscript $norm_script colon/DLD1/GSE105318/ENCFF714TMN/GSE105318_ENCFF714TMN_chromatin_interactions_hg19_{}_TopDom.matrix {} 40000 40000 colon/DLD1/GSE105318/ENCFF714TMN/GSE105318_ENCFF714TMN_chromatin_interactions_hg19_{}_ICE_TopDom.matrix" -- ${all_chromos[@]}
parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript $norm_script colon/DLD1/GSE105318/ENCFF714TMN/GSE105318_ENCFF714TMN_chromatin_interactions_hg19_{}_TopDom.matrix {} 40000 40000 colon/DLD1/GSE105318/ENCFF714TMN/GSE105318_ENCFF714TMN_chromatin_interactions_hg19_{}_ICE_TopDom.matrix" -- ${all_chromos[@]}


#parallel -i -j $maxJobs -l $maxLoad sh -c "echo Rscript $norm_script astrocyte/cerebellum/GSE105194/ENCFF122YID/GSE105194_ENCFF122YID_chromatin_interactions_hg19_{}_TopDom.matrix {} 40000 40000 astrocyte/cerebellum/GSE105194/ENCFF122YID/GSE105194_ENCFF122YID_chromatin_interactions_hg19_{}_ICE_TopDom.matrix" -- ${all_chromos[@]}
#parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript $norm_script astrocyte/cerebellum/GSE105194/ENCFF122YID/GSE105194_ENCFF122YID_chromatin_interactions_hg19_{}_TopDom.matrix {} 40000 40000 astrocyte/cerebellum/GSE105194/ENCFF122YID/GSE105194_ENCFF122YID_chromatin_interactions_hg19_{}_ICE_TopDom.matrix" -- ${all_chromos[@]}


#parallel -i -j $maxJobs -l $maxLoad sh -c "echo Rscript $norm_script astrocyte/spinal_cord/GSE105957/ENCFF478UBU/GSE105957_ENCFF478UBU_chromatin_interactions_hg19_{}_TopDom.matrix {} 40000 40000 astrocyte/spinal_cord/GSE105957/ENCFF478UBU/GSE105957_ENCFF478UBU_chromatin_interactions_hg19_{}_ICE_TopDom.matrix" -- ${all_chromos[@]}
#parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript $norm_script astrocyte/spinal_cord/GSE105957/ENCFF478UBU/GSE105957_ENCFF478UBU_chromatin_interactions_hg19_{}_TopDom.matrix {} 40000 40000 astrocyte/spinal_cord/GSE105957/ENCFF478UBU/GSE105957_ENCFF478UBU_chromatin_interactions_hg19_{}_ICE_TopDom.matrix" -- ${all_chromos[@]}


# Rscript norm_int_data.R colon/DLD1/GSE105318/ENCFF714TMN/GSE105318_ENCFF714TMN_chromatin_interactions_hg19_chr1_TopDom.matrix chr1 40000 40000 colon/DLD1/GSE105318/ENCFF714TMN/GSE105318_ENCFF714TMN_chromatin_interactions_hg19_chr1_ICE_TopDom.matrix

# Rscript norm_int_data.R astrocyte/cerebellum/GSE105194/ENCFF122YID/GSE105194_ENCFF122YID_chromatin_interactions_hg19_chr1_TopDom.matrix chr1 40000 40000 astrocyte/cerebellum/GSE105194/ENCFF122YID/GSE105194_ENCFF122YID_chromatin_interactions_hg19_chr1_ICE_TopDom.matrix

# Rscript norm_int_data.R astrocyte/spinal_cord/GSE105957/ENCFF478UBU/GSE105957_ENCFF478UBU_chromatin_interactions_hg19_chr1_TopDom.matrix chr1 40000 40000 astrocyte/spinal_cord/GSE105957/ENCFF478UBU/GSE105957_ENCFF478UBU_chromatin_interactions_hg19_chr1_ICE_TopDom.matrix



###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################
echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0








