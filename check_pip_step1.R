
DE_dt1_f <- "/media/electron/mnt/etemp/marie/Dixon2018_integrative_data/PIPELINE/OUTPUT_FOLDER/1_runGeneDE/DE_topTable.Rdata"
DE_dt2_f <- "/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER/TCGAbrca_lum_bas/1_runGeneDE/DE_topTable.Rdata"

DE_dt1 <- eval(parse(text = load(DE_dt1_f)))
DE_dt2 <- eval(parse(text = load(DE_dt2_f)))

DE_dt1[1:5,1:5]
DE_dt2[1:5,1:5]
stopifnot(all.equal(DE_dt1,DE_dt2))
