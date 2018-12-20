#!/bin/bash

# ./cmds_h5.sh

start_time=$(date -R)    


convert_script="/mnt/etemp/marie/Dixon2018_integrative_data/convert_h5_TopDom_matrix.R"

all_genome_files=(
"astrocyte/cerebellum/GSE105194/ENCFF027IEO/GSE105194_ENCFF027IEO_chromatin_interactions_hg19.h5"
"astrocyte/cerebellum/GSE105194/ENCFF031NDI/GSE105194_ENCFF031NDI_chromatin_interactions_hg19.h5"
"astrocyte/cerebellum/GSE105194/ENCFF094JAG/GSE105194_ENCFF094JAG_chromatin_interactions_hg19.h5"
"astrocyte/cerebellum/GSE105194/ENCFF122YID/GSE105194_ENCFF122YID_chromatin_interactions_hg19.h5"
"astrocyte/cerebellum/GSE105194/ENCFF241JZG/GSE105194_ENCFF241JZG_chromatin_interactions_hg19.h5"
"astrocyte/cerebellum/GSE105194/ENCFF363ZZX/GSE105194_ENCFF363ZZX_chromatin_interactions_hg19.h5"
"astrocyte/cerebellum/GSE105194/ENCFF497EDU/GSE105194_ENCFF497EDU_chromatin_interactions_hg19.h5"
"astrocyte/cerebellum/GSE105194/ENCFF526GSE/GSE105194_ENCFF526GSE_chromatin_interactions_hg19.h5"
"astrocyte/cerebellum/GSE105194/ENCFF652CHM/GSE105194_ENCFF652CHM_chromatin_interactions_hg19.h5"
"astrocyte/cerebellum/GSE105194/ENCFF693YEZ/GSE105194_ENCFF693YEZ_chromatin_interactions_hg19.h5"
"astrocyte/cerebellum/GSE105194/ENCFF739HLL/GSE105194_ENCFF739HLL_chromatin_interactions_hg19.h5"
"astrocyte/cerebellum/GSE105194/ENCFF754UZR/GSE105194_ENCFF754UZR_chromatin_interactions_hg19.h5"
"astrocyte/cerebellum/GSE105194/ENCFF782WTU/GSE105194_ENCFF782WTU_chromatin_interactions_hg19.h5"
"astrocyte/cerebellum/GSE105194/ENCFF870NPA/GSE105194_ENCFF870NPA_chromatin_interactions_hg19.h5"
)


all_chromos=( "chr"{1..22} "chrX" )
#all_chromos=( "chr1" )
#all_chromos=( "chr"{2..5} )
#all_chromos=( "chr22" )
#all_chromos=( "chr1" )

maxJobs=1
maxLoad=70

#for genome_h5file in "${all_genome_files[@]}"; do

#    parallel -i -j $maxJobs -l $maxLoad sh -c "echo Rscript $convert_script $genome_h5file {}" -- ${all_chromos[@]}
#    parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript $convert_script $genome_h5file {}" -- ${all_chromos[@]}

#done


for chromo in "${all_chromos[@]}"; do
for genome_h5file in "${all_genome_files[@]}"; do

    echo Rscript $convert_script $genome_h5file $chromo
    Rscript $convert_script $genome_h5file $chromo

done
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

    
#***** COLORECTAL
#all_genome_files=(
#"/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF115ORD/GSE105318_ENCFF115ORD_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF207ZJW/GSE105318_ENCFF207ZJW_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF210VTY/GSE105318_ENCFF210VTY_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF413KFQ/GSE105318_ENCFF413KFQ_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF439QFU/GSE105318_ENCFF439QFU_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF467CIO/GSE105318_ENCFF467CIO_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF572WDP/GSE105318_ENCFF572WDP_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF590XFM/GSE105318_ENCFF590XFM_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF714TMN/GSE105318_ENCFF714TMN_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF734ZXG/GSE105318_ENCFF734ZXG_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF811XOQ/GSE105318_ENCFF811XOQ_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF846AVK/GSE105318_ENCFF846AVK_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF872NJE/GSE105318_ENCFF872NJE_chromatin_interactions_hg19.h5"
#"/mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF993AZB/GSE105318_ENCFF993AZB_chromatin_interactions_hg19.h5"
#)


#***** ASTROCYTES - CEREBELLUM



#***** ASTROCYTES - SPINAL CORD
#all_genome_files=(
#"astrocyte/spinal_cord/GSE105957/ENCFF009UQL/GSE105957_ENCFF009UQL_chromatin_interactions_hg19.h5"
#"astrocyte/spinal_cord/GSE105957/ENCFF018DYF/GSE105957_ENCFF018DYF_chromatin_interactions_hg19.h5"
#"astrocyte/spinal_cord/GSE105957/ENCFF032URD/GSE105957_ENCFF032URD_chromatin_interactions_hg19.h5"
#"astrocyte/spinal_cord/GSE105957/ENCFF256BWW/GSE105957_ENCFF256BWW_chromatin_interactions_hg19.h5"
#"astrocyte/spinal_cord/GSE105957/ENCFF478UBU/GSE105957_ENCFF478UBU_chromatin_interactions_hg19.h5"
#"astrocyte/spinal_cord/GSE105957/ENCFF497BFM/GSE105957_ENCFF497BFM_chromatin_interactions_hg19.h5"
#"astrocyte/spinal_cord/GSE105957/ENCFF533CBL/GSE105957_ENCFF533CBL_chromatin_interactions_hg19.h5"
#"astrocyte/spinal_cord/GSE105957/ENCFF600FGH/GSE105957_ENCFF600FGH_chromatin_interactions_hg19.h5"
#"astrocyte/spinal_cord/GSE105957/ENCFF604UTQ/GSE105957_ENCFF604UTQ_chromatin_interactions_hg19.h5"
#"astrocyte/spinal_cord/GSE105957/ENCFF712GRL/GSE105957_ENCFF712GRL_chromatin_interactions_hg19.h5"
#"astrocyte/spinal_cord/GSE105957/ENCFF715HDW/GSE105957_ENCFF715HDW_chromatin_interactions_hg19.h5"
#"astrocyte/spinal_cord/GSE105957/ENCFF849SCW/GSE105957_ENCFF849SCW_chromatin_interactions_hg19.h5"
#"astrocyte/spinal_cord/GSE105957/ENCFF941LBD/GSE105957_ENCFF941LBD_chromatin_interactions_hg19.h5"
#"astrocyte/spinal_cord/GSE105957/ENCFF969PYA/GSE105957_ENCFF969PYA_chromatin_interactions_hg19.h5"
#)
