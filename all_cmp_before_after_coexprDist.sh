#!/usr/bin/bash

# ./all_cmp_before_after_coexprDist.sh

# Rscript cmp_before_after_coexprDist.R GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal TCGAgbm_classical_mesenchymal #
# Rscript cmp_before_after_coexprDist.R GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal TCGAgbm_classical_neural #
# Rscript cmp_before_after_coexprDist.R GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal TCGAgbm_classical_proneural #

# COLORECTAL
# Rscript cmp_before_after_coexprDist.R GSE105318_ENCFF439QFU_DLD1 TCGAcoad_msi_mss # no change

# BREAST
# Rscript cmp_before_after_coexprDist.R GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP TCGAbrca_lum_bas #

# KIDNEY
# Rscript cmp_before_after_coexprDist.R GSE105465_ENCFF777DUA_Caki2_vs_GSE105235_ENCFF235TGH_G401 TCGAkich_norm_kich #

# LUNG
# Rscript cmp_before_after_coexprDist.R GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460 TCGAluad_mutKRAS_mutEGFR #
# Rscript cmp_before_after_coexprDist.R GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460 TCGAluad_nonsmoker_smoker #
# Rscript cmp_before_after_coexprDist.R GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460 TCGAluad_wt_mutKRAS #
# Rscript cmp_before_after_coexprDist.R GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460 TCGAlusc_norm_lusc #

# SKIN
# Rscript cmp_before_after_coexprDist.R GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5 TCGAskcm_lowInf_highInf #
# Rscript cmp_before_after_coexprDist.R GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5 TCGAskcm_wt_mutBRAF #
# Rscript cmp_before_after_coexprDist.R GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5 TCGAskcm_wt_mutCTNNB1 #

# PANCREAS
# Rscript cmp_before_after_coexprDist.R GSE105566_ENCFF358MNA_Panc1 TCGApaad_wt_mutKRAS


all_ds=(
# "GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal TCGAgbm_classical_mesenchymal"
# "GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal TCGAgbm_classical_neural"
# "GSE105194_ENCFF027IEO_astroCerebellum_vs_GSE105957_ENCFF715HDW_astroSpinal TCGAgbm_classical_proneural"

# COLORECTAL
"GSE105318_ENCFF439QFU_DLD1 TCGAcoad_msi_mss" #no change

# BREAST
# "GSM1631185_MCF7_vs_GSE75070_MCF7_shGFP TCGAbrca_lum_bas"
# 
# # KIDNEY
# "GSE105465_ENCFF777DUA_Caki2_vs_GSE105235_ENCFF235TGH_G401 TCGAkich_norm_kich"
# 
# # LUNG
# "GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460 TCGAluad_mutKRAS_mutEGFR"
# "GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460 TCGAluad_nonsmoker_smoker"
# "GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460 TCGAluad_wt_mutKRAS"
# "GSE105600_ENCFF852YOE_A549_vs_GSE105725_ENCFF697NNX_NCIH460 TCGAlusc_norm_lusc"
# 
# # SKIN
# "GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5 TCGAskcm_lowInf_highInf"
# "GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5 TCGAskcm_wt_mutBRAF"
# "GSE106022_ENCFF614EKT_RPMI7951_vs_GSE105491_ENCFF458OWO_SKMEL5 TCGAskcm_wt_mutCTNNB1"
# 
# # PANCREAS
# "GSE105566_ENCFF358MNA_Panc1 TCGApaad_wt_mutKRAS"
)

for ds in "${all_ds[@]}"; do
  cmd="Rscript cmp_before_after_coexprDist.R $ds"
  echo $cmd
  $cmd

done