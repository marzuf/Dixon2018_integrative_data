#!/usr/bin/bash




#### run 22.12.2018: (14 ds)

# GBM
# ./run_pipeline.sh GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal TCGAgbm_classical_mesenchymal
# ./run_pipeline.sh GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal TCGAgbm_classical_neural
# ./run_pipeline.sh GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal TCGAgbm_classical_proneural

# COLORECTAL
# ./run_pipeline.sh GSE105318_ENCFF439QFU_DLD1 TCGAcoad_msi_mss

# BREAST
# ./run_pipeline.sh GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP TCGAbrca_lum_bas

# KIDNEY
# ./run_pipeline.sh GSE105465_ENCFF777DUA_Caki2_vs_GSE105235_ENCFF235TGH_G401 TCGAkich_norm_kich

# LUNG
# ./run_pipeline.sh GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460 TCGAluad_mutKRAS_mutEGFR
# ./run_pipeline.sh GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460 TCGAluad_nonsmoker_smoker
# ./run_pipeline.sh GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460 TCGAluad_wt_mutKRAS
# ./run_pipeline.sh GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460 TCGAlusc_norm_lusc

# SKIN
# ./run_pipeline.sh GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5 TCGAskcm_lowInf_highInf
# ./run_pipeline.sh GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5 TCGAskcm_wt_mutBRAF
# ./run_pipeline.sh GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5 TCGAskcm_wt_mutCTNNB1


# PANCREAS

# ./run_pipeline.sh GSE105566_ENCFF358MNA_Panc1 TCGApaad_wt_mutKRAS

# PROSTATE


start_time=$(date -R)    
set -e

if [[ $# != 2 ]]; then
    echo "invalid # of arguments"
    exit 1
fi


runDir="/mnt/etemp/marie/Dixon2018_integrative_data"

geneDataDir="$runDir/gene_data_final"

TAD_DE_pipDir="/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom"
TAD_DE_script="./zzz_run_given_step_given_data_v2.sh"

function mvBack {
  echo "... go back to my folder"
  cd $runDir  
}
trap mvBack EXIT


echo_and_launch() {
   echo "> $1"
   $1
}

#cmd="pwd"
#echo_and_launch "pwd"
#exit 0

#all_expr_datasets=(
#"TCGAbrca_lum_bas"
#"GSE58135_ERpos_tripleNeg"
#"GSE58135_ERpos_adjERpos"
#"GSE58135_tripleNeg_adjTripleNeg"
#)

#all_hic_datasets=(
#"GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460"
#"GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP"
#"GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP_vs_GSE105697_ENCFF364CWZ_T47D"
#)

hic_dataset="$1"
expr_dataset="$2"

echo "*** START ***"
echo "... > Hi-C dataset: $hic_dataset"
echo "... > Gene expression dataset: $expr_dataset"

old_inputFolder="/mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput"


#********************** HARD-CODED SETTINGS ********************************************
Rexec=`which Rscript`

new_inputFolder="$runDir/PIPELINE/INPUT_FILES/$hic_dataset"
mkdir -p $new_inputFolder
outputFolder="$runDir/PIPELINE/OUTPUT_FOLDER/$hic_dataset/$expr_dataset"
mkdir -p $outputFolder

nPermut="10000"
ncpu="20"

step1=1     # prepare setting file
step2=1     # run the pipeline

#TAD_DE_pipSteps=( "0cleanInputTCGA" "1cleanInput" )
TAD_DE_pipSteps=( "0cleanInputTCGA" "1cleanInput" "2" "3" "5" "4" "6" "7" "8c" "9" "10" "11" "13cleanInput" "14f2" "170revision2EZH2" )
#TAD_DE_pipSteps=( "0cleanInput" "1cleanInput" "2" "3" "5" )
#TAD_DE_pipSteps=( "0cleanInput" )
#TAD_DE_pipSteps=( "13cleanInput" )

#****************************************************************************************


old_setting_file="$old_inputFolder/run_settings_${expr_dataset}.R"

if [[ ! -f $old_setting_file ]] ; then
    echo "$old_setting_file does not exist"
    exit 1
fi

new_setting_file="$new_inputFolder/run_settings_${expr_dataset}.R"

echo_and_launch "cp $old_setting_file $new_inputFolder"


# $geneDataDir/$hic_dataset/genes2tad/all_assigned_regions.txt
# /mnt/etemp/marie/Dixon2018_integrative_data/gene_data_final/GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP_vs_GSE105697_ENCFF364CWZ_T47D/genes2tad/all_assigned_regions.txt
# TADposDT.txt
#chr1    chr1_TAD1       750001  1300000
new_TADpos_file="$geneDataDir/$hic_dataset/genes2tad/all_assigned_regions.txt"

# $geneDataDir/$hic_dataset/genes2tad/all_genes_positions.txt
# /mnt/etemp/marie/Dixon2018_integrative_data/gene_data_final/GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP_vs_GSE105697_ENCFF364CWZ_T47D/genes2tad/all_genes_positions.txt
#gene2tadDT.txt
#12893       chr1    761586  762902  chr1_TAD1
new_gene2tad_file="$geneDataDir/$hic_dataset/genes2tad/all_genes_positions.txt"

if [[ ! -f  $new_TADpos_file ]]; then
    echo "... new_TADpos_file $new_TADpos_file does not exist !"
    exit 1
fi

if [[ ! -f  $new_gene2tad_file ]]; then
    echo "... new_gene2tad_file $new_gene2tad_file does not exist !"
    exit 1
fi

if [[ "$step1" -eq 1 ]] ; then

	cat >> ${new_setting_file} <<- EOM

			# > file edited: `date -R` 

			# path to output folder:
			pipOutFold <- "${outputFolder}"

			# OVERWRITE THE DEFAULT SETTINGS FOR INPUT FILES - use TADs from the current Hi-C dataset 
			TADpos_file <- paste0(setDir, "`realpath $new_TADpos_file`")
							#chr1    chr1_TAD1       750001  1300000
							#chr1    chr1_TAD2       2750001 3650000
							#chr1    chr1_TAD3       3650001 4150000

			gene2tadDT_file <- paste0(setDir, "`realpath $new_gene2tad_file`")
							#LINC00115       chr1    761586  762902  chr1_TAD1
							#FAM41C  chr1    803451  812283  chr1_TAD1
							#SAMD11  chr1    860260  879955  chr1_TAD1
							#NOC2L   chr1    879584  894689  chr1_TAD1

			# overwrite main_settings.R: nCpu <- 25
			nCpu <- ${ncpu}

			# *************************************************************************************************************************
			# ************************************ SETTINGS FOR PERMUTATIONS (5#_, 8c_)
			# *************************************************************************************************************************

			# number of permutations
			nRandomPermut <- $nPermut
			gene2tadAssignMethod <- "maxOverlap"
			nRandomPermutShuffle <- $nPermut
			step8_for_permutGenes <- TRUE
			step8_for_randomTADsFix <- FALSE
			step8_for_randomTADsGaussian <- FALSE
			step8_for_randomTADsShuffle <- FALSE
			step14_for_randomTADsShuffle <- FALSE

	EOM
	echo "WRITTEN: ${new_setting_file}"


echo "> END STEP1:" $(date -R)

fi # end if STEP1

###################################################################################################################################################
### STEP2: MOVE TO THE TAD DE PIPELINE DIRECTORY TO LAUNCH THE TAD DE PIPELINE
###################################################################################################################################################
if [[ "$step2" -eq 1 ]] ; then
	echo "> START STEP2:" $(date -R)
	echo_and_launch "cd $TAD_DE_pipDir"

	echo $TAD_DE_script ${new_setting_file} `echo ${TAD_DE_pipSteps[*]}`
	$TAD_DE_script ${new_setting_file} `echo ${TAD_DE_pipSteps[*]}`

	cd $runDir
	echo "> END STEP2:" $(date -R) 
fi


###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################
echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

