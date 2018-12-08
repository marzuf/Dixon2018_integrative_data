#!/bin/bash


# ./run_prep_Rao.sh

start_time=$(date -R)    


rao_script="prep_Rao_data.R"


#all_chromos=( "chr"{1..22} "chrX" )
#all_chromos=( "chr21" )
all_chromos=( "chr9" )

maxJobs=40
maxLoad=70

parallel -i -j $maxJobs -l $maxLoad sh -c "echo Rscript $rao_script leukemia/K562/GSE63525/K562/10kb_resolution_intrachromosomal/{}/MAPQGE30/{}_10kb.RAWobserved {} 10000 40000 leukemia/K562/GSE63525/GSE63525_K562_40kb_ICE_{}_TopDom.matrix" -- ${all_chromos[@]}
parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript $rao_script leukemia/K562/GSE63525/K562/10kb_resolution_intrachromosomal/{}/MAPQGE30/{}_10kb.RAWobserved {} 10000 40000 leukemia/K562/GSE63525/GSE63525_K562_40kb_ICE_{}_TopDom.matrix" -- ${all_chromos[@]}



###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################
echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0








