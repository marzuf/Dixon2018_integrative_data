#!/bin/bash

exit 0

logfile="liver_buildHiC_logfile.txt"
rm -f $logfile
echo "--------------------------------------------------------------------------------------" >> ${logfile}
start_time=$(date -R)   
echo $start_time 
echo $start_time >> ${logfile}
echo "----- Script starts -----" >> ${logfile}


### Actual script


genomeFile="/mnt/ndata/daniele/hpv_integration/Data/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa"

bwaBin="/mnt/ed4/marie/software/bwa-0.7.17/bwa"


# -k INT        minimum seed length [19]
#       -t INT        number of threads [1]
#       -m INT        perform at most INT rounds of mate rescues for each read [50]


all_files_1=("ENCLB022KPF_ENCFF419ZIV_a1" "ENCLB022KPF_ENCFF564ZRQ_b1" "ENCLB022KPF_ENCFF652ZBK_c1" "ENCLB022KPF_ENCFF613EQX_d1" "ENCLB022KPF_ENCFF779SUW_e1" "ENCLB625TGE_ENCFF932XVZ_a1" "ENCLB625TGE_ENCFF181XUK_b1" "ENCLB625TGE_ENCFF385SWQ_c1" "ENCLB625TGE_ENCFF973IEP_d1" "ENCLB625TGE_ENCFF846FHX_e1" )
all_files_2=("ENCLB022KPF_ENCFF122SLQ_a2" "ENCLB022KPF_ENCFF428VHR_b2" "ENCLB022KPF_ENCFF568YFJ_c2" "ENCLB022KPF_ENCFF969ZME_d2" "ENCLB022KPF_ENCFF936TOC_e2" "ENCLB625TGE_ENCFF024QBA_a2" "ENCLB625TGE_ENCFF574TXX_b2" "ENCLB625TGE_ENCFF424TRG_c2" "ENCLB625TGE_ENCFF602QHH_d2" "ENCLB625TGE_ENCFF616FUR_e2" )

len1=${#all_files_1[@]}
len2=${#all_files_2[@]}

if [[ $len1 -ne $len2 ]]; then
    echo "ERROR: not same number of files"
    exit 1
fi


for i in "${!all_files_1[@]}"; do 


    file1="${all_files_1[$i]}"
    file2="${all_files_2[$i]}"



    echo "*** start $file1"
    echo "*** start $file2" >> $logfile

    # Daniele params
#    echo $bwaBin mem -t 16 -M -k 20 $genomeFile ${run}_${nbr}.fastq -o ${run}_${nbr}.sam
#    $bwaBin mem -t 16 -M -k 20 $genomeFile ${run}_${nbr}.fastq -o ${run}_${nbr}.sam

            #$ bwa mem -A1 -B4  -E50 -L0  index_path \
            #    -U mate_R1.fastq.gz ) 2>>mate_R1.log | samtools view -Shb - > mate_R1.bam

    echo "$bwaBin mem -A1 -B4  -E50 -L0  $genomeFile ${file1}.fastq | samtools view -Shb - > ${run}.bam"
    $bwaBin mem -A1 -B4  -E50 -L0  $genomeFile ${file1}.fastq | samtools view -Shb - > ${run}.bam
    done


done 

echo "*** DONE"
end_time=$(date -R)    
echo $end_time
echo $end_time >> ${logfile}
echo "start: $start_time" >> ${logfile}
echo "end: $end_time" >> ${logfile}
echo "----- Script ends -----" >> ${logfile}
exit 0


#*************** OTHER COMMANDS



# qc_template.html was missing in /home/marie/.local/lib/python2.7/site-packages/hicexplorer/
# download from github and copy paste in the folder

#Schmitt data restriction enzyme: HindIII AAGCTT



## ! wrong order -> SHOULD USE cat NOT MERGE TO AVOID REORDERING
##samtools merge hicExpParam_31_32_33_34_1.bam *_1.bam
##samtools merge hicExpParam_31_32_33_34_2.bam *_2.bam
## hicBuildMatrix --samFiles hicExpParam_31_32_33_34_1.bam hicExpParam_31_32_33_34_2.bam \
##                 --outBam GSM2322555_hicExpParam.bam \
##                 --binSize 10000 \
##                 --restrictionSequence AAGCTT \
##                 --outFileName GSM2322555_hicExpParam_10kb.npz \
##                --QCfolder GSM2322555_hicExpParam_hicQC > GSM2322555_hicExpParam.log
## /home/marie/.local/lib/python2.7/site-packages/hicexplorer/hicBuildMatrix.py", line 652, in main
##    "the --reorder option".format(mate1.qname, mate2.qname)
##AssertionError: FATAL ERROR SRR4272032.2 SRR4272031.2 Be sure that the sam files have the same read order If using Bowtie2 or Hisat2 add the --reorder option

#rm -f hicExpParam_31_32_33_34_1.bam
#rm -f hicExpParam_31_32_33_34_2.bam

#************************************************************
#********** CONCATENATE THE BAM FILES FROM DIFFERENT RUNS - HiCExplorer parameters
#************************************************************

#samtools cat -o cat_hicExpParam_31_32_33_34_1.bam SRR*_1.bam
#samtools cat -o cat_hicExpParam_31_32_33_34_2.bam SRR*_2.bam

#************************************************************
#*************** BUILD MATRIX OF A GIVEN BIN SIZE FROM BAM FILES - HiCExplorer parameters
#************************************************************

# hicBuildMatrix --samFiles cat_hicExpParam_31_32_33_34_1.bam cat_hicExpParam_31_32_33_34_2.bam \
#                 --outBam GSM2322555_cat_hicExpParam.bam \
#                 --binSize 10000 \
#                 --restrictionSequence AAGCTT \
#                 --outFileName GSM2322555_cat_hicExpParam_10kb.npz \
#                --QCfolder GSM2322555_cat_hicExpParam_hicQC > GSM2322555_cat_hicExpParam.log

# hicBuildMatrix --samFiles cat_hicExpParam_31_32_33_34_1.bam cat_hicExpParam_31_32_33_34_2.bam \
#                 --outBam GSM2322555_cat_hicExpParam_20kb.bam \
#                 --binSize 20000 \
#                 --restrictionSequence AAGCTT \
#                 --outFileName GSM2322555_cat_hicExpParam_20kb.npz \
#                --QCfolder GSM2322555_cat_hicExpParam_hicQC_20kb > GSM2322555_cat_hicExpParam_20kb.log


#************************************************************
#*************** CONVERT NPZ.H5 TO DEKKER FORMAT - HiCExplorer parameters
#************************************************************

#hicExport --inFile GSM2322555_cat_hicExpParam_20kb.npz.h5 \
#-o hicExpParam_20kb_matrix/GSM2322555_hicExpParam_20kb_dekker_matrix.txt --outputFormat dekker


hicExport --inFile GSM2322555_cat_hicExpParam_20kb.npz.h5 \
--chrNameList 1 -o hicExpParam_20kb_matrix/GSM2322555_hicExpParam_20kb_chr1_dekker_matrix.txt --outputFormat dekker
# running


# -> OK ! - final file: /mnt/etemp/marie/GSE87112/SB_GSM2322555/hicExpParam_20kb_matrix/GSM2322555_hicExpParam_20kb_dekker_matrix.txt.gz
gunzip  /mnt/etemp/marie/GSE87112/SB_GSM2322555/hicExpParam_20kb_matrix/GSM2322555_hicExpParam_20kb_dekker_matrix.txt.gz

hicExport --inFile GSM2322555_cat_hicExpParam_10kb.npz.h5 \
-o hicExpParam_10kb_matrix/GSM2322555_hicExpParam_10kb_dekker_matrix.txt --outputFormat dekker
# => running

#====================================================================================================================================================================================
#====================================================================================================================================================================================
#====================================================================================================================================================================================
#====================================================================================================================================================================================
#====================================================================================================================================================================================

#************************************************************
#********** CONVERT BAM TO SAM FILES - Daniele parameters
#************************************************************

# in SAM_DanieleParam folder
#samtools view -Sb -o SRR4272031_1.bam SRR4272031_1.sam
#samtools view -Sb -o SRR4272031_2.bam SRR4272031_2.sam
#samtools view -Sb -o SRR4272032_1.bam SRR4272032_1.sam
#samtools view -Sb -o SRR4272032_2.bam SRR4272032_2.sam
#samtools view -Sb -o SRR4272033_1.bam SRR4272033_1.sam
#samtools view -Sb -o SRR4272033_2.bam SRR4272033_2.sam
#samtools view -Sb -o SRR4272034_1.bam SRR4272034_1.sam
#samtools view -Sb -o SRR4272034_2.bam SRR4272034_2.sam

#************************************************************
#********** CONCATENATE THE BAM FILES FROM DIFFERENT RUNS - Daniele parameters
#************************************************************

#samtools cat -o cat_danieleParam_31_32_33_34_1.bam SRR*_1.bam
#samtools cat -o cat_danieleParam_31_32_33_34_2.bam SRR*_2.bam


#************************************************************
#*************** BUILD MATRIX OF A GIVEN BIN SIZE FROM BAM FILES - Daniele parameters
#************************************************************

# hicBuildMatrix --samFiles cat_danieleParam_31_32_33_34_1.bam cat_danieleParam_31_32_33_34_2.bam \
#                 --outBam GSM2322555_cat_danieleParam_20kb.bam \
#                 --binSize 20000 \
#                 --restrictionSequence AAGCTT \
#                 --outFileName GSM2322555_cat_danieleParam_20kb.npz \
#                --QCfolder GSM2322555_cat_danieleParam_hicQC_20kb > GSM2322555_cat_danieleParam_20kb.log

# hicBuildMatrix --samFiles cat_danieleParam_31_32_33_34_1.bam cat_danieleParam_31_32_33_34_2.bam \
#                 --outBam GSM2322555_cat_danieleParam_10kb.bam \
#                 --binSize 10000 \
#                 --restrictionSequence AAGCTT \
#                 --outFileName GSM2322555_cat_danieleParam_10kb.npz \
#                --QCfolder GSM2322555_cat_danieleParam_hicQC_10kb > GSM2322555_cat_danieleParam_10kb.log


#************************************************************
#*************** CONVERT NPZ.H5 TO DEKKER FORMAT - Daniele parameters
#************************************************************

#hicExport --inFile GSM2322555_cat_danieleParam_20kb.npz.h5 \
#-o danieleParam_20kb_matrix/GSM2322555_danieleParam_20kb_dekker_matrix.txt --outputFormat dekker


hicExport --inFile GSM2322555_cat_danieleParam_20kb.npz.h5 \
--chrNameList 1 -o danieleParam_20kb_matrix/GSM2322555_danieleParam_20kb_chr1_dekker_matrix.txt --outputFormat dekker

# -> OK ! - final file: /mnt/etemp/marie/GSE87112/SB_GSM2322555/SAM_DanieleParam/danieleParam_20kb_matrix/GSM2322555_danieleParam_20kb_dekker_matrix.txt.gz
gunzip  /mnt/etemp/marie/GSE87112/SB_GSM2322555/SAM_DanieleParam/danieleParam_20kb_matrix/GSM2322555_danieleParam_20kb_dekker_matrix.txt.gz

hicExport --inFile GSM2322555_cat_danieleParam_10kb.npz.h5 \
-o danieleParam_10kb_matrix/GSM2322555_danieleParam_10kb_dekker_matrix.txt --outputFormat dekker

# ERROR:
#saving...
#Traceback (most recent call last):
#  File "/home/marie/.local/bin/hicExport", line 7, in <module>
#    main()
#  File "/home/marie/.local/lib/python2.7/site-packages/hicexplorer/hicExport.py", line 198, in main
#    hic_ma.save_dekker(args.outFileName)
#  File "/home/marie/.local/lib/python2.7/site-packages/hicexplorer/HiCMatrix.py", line 1010, in save_dekker
#    raise msg
#TypeError: exceptions must be old-style classes or derived from BaseException, not str


