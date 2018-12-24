#!/usr/bin/bash

# ./assign_genes.sh

set -e

start_time=$(date -R)    

Rexec="Rscript"

## PARALLELIZATION
maxJobs=10
maxLoad=70

# general settings
chromo=( "chr"{1..22} "chrX" ) 
#chromo=( "chr21" ) 
ncpu="5"


mainInFolder="/mnt/etemp/marie/Dixon2018_integrative_data/FIND_CONSENSUS_TADS"

all_datasets=(
"GSE105566_ENCFF358MNA_Panc1"
)


# /mnt/etemp/marie/Dixon2018_integrative_data/colon/DLD1/GSE105318/ENCFF439QFU/TopDom/GSE105318_ENCFF439QFU_40kb_chr6.bed

### done in "gene_data_final":
#GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460
#GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP
#GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP_vs_GSE105697_ENCFF364CWZ_T47D
#GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal
#GSE105465_ENCFF777DUA_Caki2_vs_GSE105235_ENCFF235TGH_G401
#GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5
#GSE105318_ENCFF439QFU_DLD1

step5_script="/mnt/ed4/marie/gene_data_final/v2_assign/gene2TAD_consensus_version2.R"



infold_genes="/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"



### STEP5: FILES AND FOLDERS
bin_size="40000"



# HARD CODED IN R SCRIPT (used at the end of STEP5 for concatenating files)
g2t_prefix="tmp_g2t"
reg_prefix="tmp_assigned"

    for dataset in "${all_datasets[@]}"; do

        inFold="$mainInFolder/$dataset"

        echo "*** START dataset: $dataset"

        outfolder_tmp="/mnt/etemp/marie/Dixon2018_integrative_data/gene_data_final/$dataset/genes2tad/tmp"
        outfolder_final="/mnt/etemp/marie/Dixon2018_integrative_data/gene_data_final/$dataset/genes2tad"
        mkdir -p $outfolder_final
        mkdir -p $outfolder_tmp


	    parallel -i -j $maxJobs -l $maxLoad sh -c "echo Rscript $step5_script -f $infold_genes -t $inFold/{}_conservedTADs.txt -c {} -o $outfolder_tmp -b $bin_size" -- ${chromo[@]}
	    parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript $step5_script -f $infold_genes -t $inFold/{}_conservedTADs.txt -c {} -o $outfolder_tmp -b $bin_size" -- ${chromo[@]}

        echo "cat $outfolder_tmp/$g2t_prefix* > $outfolder_final/all_genes_positions.txt"
        cat $outfolder_tmp/$g2t_prefix* > $outfolder_final/all_genes_positions.txt

        echo "cat $outfolder_tmp/$reg_prefix* > $outfolder_final/all_assigned_regions.txt"
        cat $outfolder_tmp/$reg_prefix* > $outfolder_final/all_assigned_regions.txt

    done


###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time


