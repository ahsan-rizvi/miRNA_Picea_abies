#!/bin/bash
module load PDC matplotlib
module load R/4.4.1-cpeGNU-23.12
clear
echo 'Bismillah'
#############################
# This pipeline is constructed to analyse the results 
# of miRNA results and analysis
# ###########################


#############################
Pip_results_analysis='no' #Pipeline output Analysis
if [ $Pip_results_analysis = 'yes' ] # 
   then
   input=$PWD/'input/pipOutput_ZE_filtered_fasta_list.txt'
   #input=$PWD/'input/pipOutput_SE_filtered_fasta_list.txt'
   output=$PWD/'output/pipOutput_ZE_filtered_fasta_list.freq.csv'
   #output=$PWD/'output/pipOutput_SE_filtered_fasta_list.freq.csv'
      while IFS= read -r line
            do
            $PWD/script/fasta_seq_size_checker.py $line
      done < "$input" > $output
fi
#############################


############################# DGE ANalysis
DESeq2_SE_miRNA='yes'
if [ $DESeq2_SE_miRNA = 'yes' ] # 
   then
   $PWD/script/DESeq_SE.R \
           $PWD/input/SE_p_countTable.tsv \
           $PWD/input/Design_SE_miRNA.csv \
	   $PWD/output/
fi
echo $PWD/script/DESeq_SE.R

DESeq2_SE_miRNA='no'
if [ $DESeq2_SE_miRNA = 'yes' ] #
   then
   $PWD/script/DESeq_ZE.R \
           $PWD/input/ZE_p_countTable.tsv \
           $PWD/input/Design_ZE_miRNA.csv \
           $PWD/output/
fi
#############################

