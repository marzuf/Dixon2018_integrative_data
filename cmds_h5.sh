#!/bin/bash


start_time=$(date -R)    


convert_script="/mnt/etemp/marie/Dixon2018_integrative_data/convert_h5_TopDom_matrix.R"


all_genome_files=(
"/mnt/etemp/marie/Dixon2018_integrative_data/prostate/LNCaP/GSE105557/ENCFF270HJX/GSE105557_ENCFF270HJX_chromatin_interactions_hg19.h5" 
)

#all_genome_files=(
#"/mnt/etemp/marie/Dixon2018_integrative_data/kidney/G401/GSE105235/ENCFF235TGH/GSE105235_ENCFF235TGH_chromatin_interactions_hg19.h5"
#)

##all_genome_files=(
##"/mnt/etemp/marie/Dixon2018_integrative_data/skin/RPMI-7951/GSE106022/ENCFF614EKT/GSE106022_ENCFF614EKT_chromatin_interactions_hg19.h5"
##)
#all_chromos=( "chr"{7..8} )
#all_chromos=( "chr"{19..20} )
#all_chromos=( "chr1" )

#all_genome_files=(
#"/mnt/etemp/marie/Dixon2018_integrative_data/skin/SK-MEL-5/GSE105491/ENCFF458OWO/GSE105491_ENCFF458OWO_chromatin_interactions_hg19.h5"
#)

#all_genome_files=(
#"/mnt/etemp/marie/Dixon2018_integrative_data/pancreas/Panc1/GSE105566/ENCFF358MNA/GSE105566_ENCFF358MNA_chromatin_interactions_hg19.h5" 
#)

#all_genome_files=(
#"/mnt/etemp/marie/Dixon2018_integrative_data/kidney/Caki2/GSE105465/ENCFF777DUA/GSE105465_ENCFF777DUA_chromatin_interactions_hg19.h5"
#)



#all_chromos=( "chr"{1..22} "chrX" )
#all_chromos=( "chr1" )
#all_chromos=( "chr"{2..5} )
#all_chromos=( "chr22" )
all_chromos=( "chr1" )

maxJobs=1
maxLoad=70

for genome_h5file in "${all_genome_files[@]}"; do

    parallel -i -j $maxJobs -l $maxLoad sh -c "echo Rscript $convert_script $genome_h5file {}" -- ${all_chromos[@]}
    parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript $convert_script $genome_h5file {}" -- ${all_chromos[@]}

done

###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################
echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

############################################################################################################################################################################ (already run)

#****** LUNG
#all_genome_files=(
#"/mnt/etemp/marie/Dixon2018_integrative_data/lung/A549/ENCSR444WCZ/GSE105600_ENCFF852YOE_chromatin_interactions_hg19.h5"
#)
#all_chromos=( "chr7" "chr8" "chr12" "chr17" "chr18" "chr21" )
#all_chromos=( "chr2" )

#all_genome_files=(
#"/mnt/etemp/marie/Dixon2018_integrative_data/lung/NCI-H460/ENCSR489OCU/GSE105725_ENCFF697NNX_chromatin_interactions_hg19.h5"
#)
#all_chromos=( "chr1" "chr4" "chr5" "chr6" "chr7" "chr9" "chr10" "chr12" "chr14" "chr15" "chr16" "chr17" "chr18" "chr20" "chr21" "chr22" )
#all_chromos=( "chr1" "chr7" "chr15" "chr20" )

#****** BREAST
#all_genome_files=(
#"/mnt/etemp/marie/Dixon2018_integrative_data/breast/T47D/ENCSR549MGQ_GSE105697/GSE105697_ENCFF364CWZ_chromatin_interactions_hg19.h5"
#)
#all_chromos=( "chr10" "chr12" "chr14" "chr16" "chr17" "chr21" "chr22" )
#all_chromos=( "chr14" "chr17" "chr19" )


#****** PANCREAS

#all_genome_files=(
#"/mnt/etemp/marie/Dixon2018_integrative_data/pancreas/Panc1/GSE105566/ENCFF276YKV/GSE105566_ENCFF276YKV_chromatin_interactions_hg19.h5" 
#"/mnt/etemp/marie/Dixon2018_integrative_data/pancreas/Panc1/GSE105566/ENCFF358MNA/GSE105566_ENCFF358MNA_chromatin_interactions_hg19.h5" 
#"/mnt/etemp/marie/Dixon2018_integrative_data/pancreas/Panc1/GSE105566/ENCFF463HGQ/GSE105566_ENCFF463HGQ_chromatin_interactions_hg19.h5" 
#)

#****** PROSTATE

#all_genome_files=(
#"/mnt/etemp/marie/Dixon2018_integrative_data/prostate/LNCaP/GSE105557/ENCFF270HJX/GSE105557_ENCFF270HJX_chromatin_interactions_hg19.h5" 
#"/mnt/etemp/marie/Dixon2018_integrative_data/prostate/LNCaP/GSE105557/ENCFF284HII/GSE105557_ENCFF284HII_chromatin_interactions_hg19.h5" 
#"/mnt/etemp/marie/Dixon2018_integrative_data/prostate/LNCaP/GSE105557/ENCFF702VDY/GSE105557_ENCFF702VDY_chromatin_interactions_hg19.h5" 
#)

#****** KIDNEY

#all_genome_files=(
#"/mnt/etemp/marie/Dixon2018_integrative_data/kidney/Caki2/GSE105465/ENCFF762DPV/GSE105465_ENCFF762DPV_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/kidney/Caki2/GSE105465/ENCFF777DUA/GSE105465_ENCFF777DUA_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/kidney/Caki2/GSE105465/ENCFF796ONA/GSE105465_ENCFF796ONA_chromatin_interactions_hg19.h5"
#)


#all_genome_files=(
#"/mnt/etemp/marie/Dixon2018_integrative_data/kidney/G401/GSE105235/ENCFF235TGH/GSE105235_ENCFF235TGH_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/kidney/G401/GSE105235/ENCFF298ZFN/GSE105235_ENCFF298ZFN_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/kidney/G401/GSE105235/ENCFF905WIG/GSE105235_ENCFF905WIG_chromatin_interactions_hg19.h5"
#)


#***** SKIN
#all_genome_files=(
#"/mnt/etemp/marie/Dixon2018_integrative_data/skin/RPMI-7951/GSE106022/ENCFF357UBI/GSE106022_ENCFF357UBI_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/skin/RPMI-7951/GSE106022/ENCFF449CKH/GSE106022_ENCFF449CKH_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/skin/RPMI-7951/GSE106022/ENCFF614EKT/GSE106022_ENCFF614EKT_chromatin_interactions_hg19.h5"
#)


#all_genome_files=(
#"/mnt/etemp/marie/Dixon2018_integrative_data/skin/SK-MEL-5/GSE105491/ENCFF064NPT/GSE105491_ENCFF064NPT_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/skin/SK-MEL-5/GSE105491/ENCFF458OWO/GSE105491_ENCFF458OWO_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/skin/SK-MEL-5/GSE105491/ENCFF605CAZ/GSE105491_ENCFF605CAZ_chromatin_interactions_hg19.h5"
#)




