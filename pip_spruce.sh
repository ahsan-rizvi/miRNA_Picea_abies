#!/bin/bash
clear
echo 'Bismillah'
#############################
# This pipeline is constructed to analyse the results 
# of miRNA results and analysis
# ###########################


#############################
Pip_results_analysis='yes' #Pipeline output Analysis
input=$PWD/'input/pipOutput_ZE_filtered_fasta_list.txt' # Input file should be constructed with the path of fasta file
#input=$PWD/'input/pipOutput_SE_filtered_fasta_list.txt'
output=$PWD/'output/pipOutput_ZE_filtered_fasta_list.freq.csv' # Define output file in CSV file
#output=$PWD/'output/pipOutput_SE_filtered_fasta_list.freq.csv'
if [ $Pip_results_analysis = 'yes' ] # 
   then
      while IFS= read -r line
            do
	    $PWD/script/fasta_seq_size_checker.py $line
      done < "$input" > $output
      echo 'Input file list for frequency check:' $input
      echo 'Output of frequency check:' $output
fi
#############################


