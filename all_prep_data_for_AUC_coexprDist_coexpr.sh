#!/usr/bin/bash

# ./all_prep_data_for_AUC_coexprDist_coexpr.sh

start_time=$(date -R)   

scriptCoexpr="create_coexpr_sortNoDup_otherTADfile.R"

all_TAD_files_ds=(
"GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal TCGAgbm_classical_mesenchymal"
"GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal TCGAgbm_classical_neural"
"GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal TCGAgbm_classical_proneural"

#"GSE105318_ENCFF439QFU_DLD1 TCGAcoad_msi_mss"

"GSE105465_ENCFF777DUA_Caki2_vs_GSE105235_ENCFF235TGH_G401 TCGAkich_norm_kich"

"GSE105566_ENCFF358MNA_Panc1 TCGApaad_wt_mutKRAS"

"GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460 TCGAluad_mutKRAS_mutEGFR"
"GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460 TCGAluad_nonsmoker_smoker"
"GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460 TCGAluad_wt_mutKRAS"
"GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460 TCGAlusc_norm_lusc"

"GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5 TCGAskcm_lowInf_highInf"
"GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5 TCGAskcm_wt_mutBRAF"
"GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5 TCGAskcm_wt_mutCTNNB1"

"GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP TCGAbrca_lum_bas"

"GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP_vs_GSE105697_ENCFF364CWZ_T47D TCGAbrca_lum_bas"
)


for ds in "${all_TAD_files_ds[@]}"; do
    
    echo Rscript $scriptCoexpr $ds
    Rscript $scriptCoexpr $ds
done


###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

