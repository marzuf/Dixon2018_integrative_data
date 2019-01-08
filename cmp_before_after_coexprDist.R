
# GBM
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

dataset="GSE105318_ENCFF439QFU_DLD1"
exprds="TCGAcoad_msi_mss"

cat("> START ", "cmp_before_after_coexprDist.R", "\n")
# Rscript cmp_before_after_coexprDist.R GSE105318_ENCFF439QFU_DLD1 TCGAcoad_msi_mss_hgnc

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
dataset=args[1]
exprds=args[2]

exprds=paste0(exprds, "_hgnc")

cat("... START: ", dataset, " - ", exprds, "\n")

cat("load genefam0 \n")
genefam0 = eval(parse(text=load("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/PREP_GENE_FAMILIES_TAD_DATA/hgnc_entrezID_family_TAD_DT.Rdata")))
cat("load samefam0 \n")
samefam0 = eval(parse(text=load("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CREATE_SAME_FAMILY_SORTNODUP/hgnc_family_all_family_pairs.Rdata")))
cat("load dist0 \n")
dist0 = eval(parse(text=load("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CREATE_DIST_SORTNODUP/all_dist_pairs.Rdata")))
cat("load coexpr0 \n")
coexpr0 = eval(parse(text=load(paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/CREATE_COEXPR_SORTNODUP/", gsub("hgnc", "", exprds), "pearson/coexprDT.Rdata"))))

cat("load genefam1 \n")
genefam1 = eval(parse(text=load(paste0("PREP_GENE_FAMILIES_TAD_DATA/", dataset, "/hgnc_entrezID_family_TAD_DT.Rdata"))))
cat("load samefam1 \n")
samefam1 = eval(parse(text=load(paste0("CREATE_SAME_FAMILY_SORTNODUP/", dataset, "/hgnc_family_short_all_family_pairs.Rdata"))))
cat("load dist1 \n")
dist1 = eval(parse(text=load(paste0( "CREATE_DIST_SORTNODUP/", dataset, "/all_dist_pairs.Rdata"))))
cat("load coexpr1 \n")
coexpr1 = eval(parse(text=load(paste0( "CREATE_COEXPR_SORTNODUP/", dataset, "/", gsub("hgnc", "", exprds), "pearson/coexprDT.Rdata"))))

cat("... dim(genefam0) = ", dim(genefam0), "\n")
cat("... dim(genefam1) = ", dim(genefam1), "\n")

cat("... dim(samefam0) = ", dim(samefam0), "\n")
cat("... dim(samefam1) = ", dim(samefam1), "\n")

cat("... dim(dist0) = ", dim(dist0), "\n")
cat("... dim(dist1) = ", dim(dist1), "\n")

cat("... dim(coexpr0) = ", dim(coexpr0), "\n")
cat("... dim(coexpr1) = ", dim(coexpr1), "\n")



### FOR THE RATIO RESULTS
cat("load ratio0 \n")
# AUC_COEXPRDIST_WITHFAM_SORTNODUP_BEFORE08.01.19_sameTAD_sameFamFile/GSE105318_ENCFF439QFU_DLD1/TCGAcoad_msi_mss_hgnc/hgnc_family_short/auc_values.Rdata
ratioFile0 <- paste0("AUC_COEXPRDIST_WITHFAM_SORTNODUP_BEFORE08.01.19_sameTAD_sameFamFile/", dataset, "/", exprds, "/hgnc_family_short/auc_values.Rdata")
ratio0 = eval(parse(text=load(ratioFile0)))
# AUC_COEXPRDIST_WITHFAM_SORTNODUP_BEFORE08.01.19_sameTAD_sameFamFile/GSE105318_ENCFF439QFU_DLD1/TCGAcoad_msi_mss_hgnc_hgnc/hgnc_family_short/auc_values.Rdata
# AUC_COEXPRDIST_WITHFAM_SORTNODUP_BEFORE08.01.19_sameTAD_sameFamFile/GSE105318_ENCFF439QFU_DLD1/TCGAcoad_msi_mss_hgnc/hgnc_family_short/auc_values.Rdata")


cat("load ratio1 \n")
ratio1 = eval(parse(text=load(paste0("AUC_COEXPRDIST_WITHFAM_SORTNODUP/", dataset, "/", exprds, "/hgnc_family_short/auc_values.Rdata"))))

all_vars <- c(
  "auc_diffTAD_distVect",                    
  "auc_sameTAD_distVect",                     
  "auc_ratio_same_over_diff_distVect", 
  "auc_diffTAD_obsDist",                
  "auc_sameTAD_obsDist",                 
  "auc_ratio_same_over_diff_obsDist",         
  "auc_sameFamDiffTAD_distVect",              
  "auc_sameFamSameTAD_distVect",              
  "auc_ratio_sameFam_same_over_diff_distVect",
  "auc_sameFamDiffTAD_obsDist",               
  "auc_sameFamSameTAD_obsDist",               
  "auc_ratio_sameFam_same_over_diff_obsDist"
) 

var="auc_ratio_sameFam_same_over_diff_obsDist"

for(var in all_vars) {
  if(ratio0[[paste0(var)]] != ratio1[[paste0(var)]]){
    cat(paste0("...... ", var, "\nratio0=", ratio0[[paste0(var)]] , "\nratio1=",ratio1[[paste0(var)]], "\n" ))
  } else{
    cat(paste0("...... ", var, "\nratio0==ratio1\n" ))
  }
}





