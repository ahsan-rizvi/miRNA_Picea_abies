#!/bin/bash -l
#SBATCH -A naiss2024-5-440
#SBATCH -p core -n 1 
#SBATCH -t 80:00:00
#SBATCH -J Long_mapping
#SBATCH -o /proj/uppstore2017145/V2/Ulrika_proj_1/err/ah-%A.out 
##############################
Proc=18
###SBATCH -N 1
##############################
Exp=$1
#Exp='SE'
#Exp='ZE'
############################## (A)
Step1='no'  #pipeline run per sample
############################## (A2)
Step2='no' #BlastN to miRBase 
Step3='no' #Count file generation per sample
############################## (B)
Step4='no'  # pre DGE analysis
DGE_miRNA='no'  #miRNA DGE analysis
Step5='no'  # post DGE analysis
############################## Combine ZE and SE analysis
Comb='no'
############################## (C)
Step6='no' #Review article check, and known miRBASE check
Heatmap='no' #heatmap, input is from Step 6
############################## (D)
mRNA_p02='no' #TAXIMPORT for mRNA (SE ZE) count table generation
Heatmap_mRNA='no'
Target='yes'
############################## (E)
if [ $Exp = 'SE' ]
  then
  fastq_list='/proj/uppstore2017145/V2/Ulrika_proj_1/sRNA_data/S39.txt'
  pip_script='miRDP2-v1.1.4_pipelineT1.bash'
  DGE_script='DESeq5_SE.R'
  DGE_design='design.txt'
  fastq_list2='/proj/uppstore2017145/V2/Ulrika_proj_1/sRNA_data/mereged_fasta.txt'
  Res_fold='/proj/uppstore2017145/V2/Ulrika_proj_1/results3'
  Res_fold22='/proj/uppstore2017145/V2/Ulrika_proj_1/results3'
  count_list='/proj/uppstore2017145/V2/Ulrika_proj_1/sRNA_data/DGE_input.m.txt'
  count_table='/proj/uppstore2017145/V2/Ulrika_proj_1/results3/count_table.m.txt'
  P_prediction_list='SE_p_prediction_list.xyz'
  DGE_gene_list='/proj/uppstore2017145/V2/Ulrika_proj_1/scripts/DGE/DGE_SE.txt'
fi
if [ $Exp = 'ZE' ]
  then
  sample='S52.txt'
  fastq_list='/proj/uppstore2017145/V2/Ulrika_proj_1/sRNA_dataZE/'$sample
  pip_script='miRDP2-v1.1.4_pipelineZ1.bash'
  fastq_list2='/proj/uppstore2017145/V2/Ulrika_proj_1/sRNA_dataZE/merged_fasta.txt'
  Res_fold='/proj/uppstore2017145/V2/Ulrika_proj_1/results2'
  Res_fold22='/proj/uppstore2017145/V2/Ulrika_proj_1/results2'
  count_list='/proj/uppstore2017145/V2/Ulrika_proj_1/sRNA_dataZE/ZEcountfilePath.txt'
  count_table='/proj/uppstore2017145/V2/Ulrika_proj_1/results2/count_table.Z.txt'
  P_prediction_list='ZE_p_prediction_list.xyz'
  DGE_script='DESeq3_ZE.R'
  DGE_design='ZE_DesignA2.csv'
  DGE_gene_list='/proj/uppstore2017145/V2/Ulrika_proj_1/scripts/DGE/DGE_ZE.txt'
fi
#############################

###########################
DGE='/proj/uppstore2017145/V2/Ulrika_proj_1/scripts/DGE'
genome='/proj/uppstore2017145/V2/users/ahsan/Picab02.fa'
script='/home/ahsan/bin/1.1.4/scripts'
index='/home/ahsan/bin/1.1.4/scripts/index'
script_loc='/proj/uppstore2017145/V2/Ulrika_proj_1/scripts/pip'
mature_loc='/proj/uppstore2017145/V2/Ulrika_proj_1/results/miRNA_predicted'
#mature_loc='/proj/uppstore2017145/V2/Ulrika_proj_1/results/miRNA_predicted2'
result='/proj/uppstore2017145/V2/Ulrika_proj_1/results/miRNA_predicted'
predict='/proj/uppstore2017145/V2/Ulrika_proj_1/results/combined_predictions_list.xyz'
target_folder='/proj/uppstore2017145/V2/Ulrika_proj_1/scripts/Target_predict'
##mkdir $Res_fold
cd $Res_fold
module load bioinfo-tools python biopython/1.80-py3.10.8 perl/5.32.1 BioPerl/1.7.8-perl5.32.1 bowtie/1.2.3 bowtie2/2.5.2 samtools/1.19 RNAfold/2.4.17 
module load blast/2.14.1+ FastQC/0.11.9 cutadapt/4.5 RNAfold/2.4.17 miRDP2/1.1.4
module load seqtk/1.2-r101 R_packages/4.3.1
clear
###########################


###########################
#Step 1: PIPELINE for miRNAs
###########################
if [ $Step1 = 'yes' ]
  then
     echo '#################PRE PROCESSING BEFORE PIPELINE ########'
     echo '#1. If sequence is in Fastq then convert it in FASTA by:'
     echo '# module load  bioinfo-tools seqtk/1.2-r101 '
     echo '# seqtk seq -a in.fastq.gz > out.fasta '
     echo '# Total count for sequence: grep -c > out.fa'
     echo '#2. Make non-redundent file for pipeline input by:'
     echo '# ../scripts/pip/pre_processing2.py out.fa '
     echo '3. list non-redundent file for pipeline into: ' $fastq_list
     echo '########################################################'
     $script_loc/$pip_script \
            -g $genome -x $genome -f -b $fastq_list \
            -o $Res_fold \
            -p $Proc 
     ls -sh $Res_fold
fi
###########################
###########################



###########################
#Step2: BLASTx with MiRBase to detect known miRNA
########################### 
Res_fold=$result #Directory for results
cd $Res_fold
miRDP2='Picea'
if [ $Step2 = 'yes' ]
   then
   echo 'HOW TO UPDATE mirbase
   #wget https://www.mirbase.org/download/mature.fa
   #bowtie-build --threads $thread $index/mature.fa $index/mature_index
   #fasta_U2T.pl $index/mature.fa $index/mature_wo_U.fa
   #unique_fasta_v1.3.pl $index/mature_wo_U.fa $index/mature_wo_U_uniq.fa mature_uniq
   #makeblastdb -in $index/mature_wo_U_uniq.fa -dbtype 'nucl''
   echo '##############################'
   parse_miRDP2_prediction.pl $predict $miRDP2
   blastn  -db $index/mature_wo_U_uniq.fa \
           -query $Res_fold/$miRDP2'_mature.fa' \
           -out $Res_fold/$miRDP2'_mature_blastn.txt' \
           -word_size 4 -num_alignments 1 #-evalue 0.05 
   general_blast_parser.pl \
           $Res_fold/$miRDP2'_mature_blastn.txt' \
           $Res_fold/$miRDP2'_mature_blastn_parsed.txt'
   parse_parsed_blast_known_plants.pl \
           $Res_fold/$miRDP2'_mature_blastn_parsed.txt' \
           $miRDP2'_mature'
   cut -f2 $Res_fold/$miRDP2'_mature_known.txt' | sort | uniq \
         > $Res_fold/$miRDP2'_mature_known_id.txt'
   cut -f2 $Res_fold/$miRDP2'_mature_variant.txt' | sort | uniq \
         > $Res_fold/temp.txt
   filter_lines_by_key_words_list.pl \
           $Res_fold/temp.txt \
           $Res_fold/$miRDP2'_mature_variant_id.txt' \
           $Res_fold/$miRDP2'_mature_known_id.txt' 0
fi
##########################



##########################
#Step3: Mapping to generate count file
########################## 
if [ $Step3 = 'yes' ]
   then
      echo 'Begin of Step3'
      bowtie-build $Res_fold/$miRDP2'_mature.fa' $Res_fold/miRNA_ref.fa
      input=$fastq_list2
         while IFS= read -r line
            do
               bowtie -p $Proc -v 0 --norc -S $Res_fold/miRNA_ref.fa \
                      -f $line | samtools view -Sb - \
                       > $line'.bam'
               bam2ref_counts.pl -bam $line'.bam' \
                      -f $result/$miRDP2'_mature.fa' \
                       > $line'_count.m.txt'  
             echo $line'_count.m.txt'
      done < "$input"

   $script_loc/combine_htseq_counts.pl $count_list $count_table
fi
#############################

############################# 08Nov2024
#TAXIMPORT for mRNA (SE ZE) count table generation
# Solman data are given by Chmelia and Elena
#############################
if [ $mRNA_p02 = 'yes' ]
   then
###########################################
Convert='/proj/uppstore2017145/V2/Ulrika_proj_1/scripts/pip/ConvertTaximport2Count.py'
input='/proj/uppstore2017145/V2/Ulrika_proj_1/salmon_SE/ListQuantFile.txt'
#echo -e '\e[0;36mInput file for above script is:\e[0m' $input
##    while IFS= read -r line
  ##        do
    ##      $Convert $line
 ##   done < "$input"
ListCountSE='/proj/uppstore2017145/V2/Ulrika_proj_1/salmon_SE/ListCountFile.txt'
##$script_loc/combine_htseq_counts.pl $ListCountSE $ListCountSE'.table'
input='/proj/uppstore2017145/V2/Ulrika_proj_1/salmon_ZE/ListCountFile.txt'
#echo -e '\e[0;36mInput file for above script is:\e[0m' $input
##    while IFS= read -r line
  ##        do
    ##      $Convert $line
##    done < "$input"
ListCountZE='/proj/uppstore2017145/V2/Ulrika_proj_1/salmon_ZE/ListCountFile.txt'
#$script_loc/combine_htseq_counts.pl $ListCountZE $ListCountZE'.table'
### find /proj/uppstore2017145/V2/Ulrika_proj_1/salmon_SE/ListCountFile.txt.table -type f -exec sed -i 's/'#'/'-'/g' {} \;
###########################################



###########################################
   if [ $1 = 'SE' ]
       then
         Script_mRNA_p02='/proj/uppstore2017145/V2/Ulrika_proj_1/scripts/pip/DESeq4_SE_mRNA_p02.R'
         DesignSE='/proj/uppstore2017145/V2/Ulrika_proj_1/scripts/pip/design_p02.txt'
         $Script_mRNA_p02 $ListCountSE'.table' $DesignSE
         #echo -e '\e[0;32mCount files list mRNA_p02 SE: \e[0m' $ListCountSE  
         echo -e '\e[0;31mScript for mRNA_p02 SE: \e[0m' $Script_mRNA_p02
         #echo -e '\e[0;32mInput count table SE: \e[0m' $ListCountSE'.table'
         #echo -e '\e[0;32mDesign file: \e[0m' $DesignSE
         echo -e '\e[0;33mResult of SE Log2FC mRNA for Target analysis: \e[0m' '/proj/uppstore2017145/V2/Ulrika_proj_1/scripts/DGE/SE_Log2FC_mRNA_p02.csv'
         
   fi

   if [ $1 = 'ZE' ]
       then
         Script_mRNA_p02_ZE='/proj/uppstore2017145/V2/Ulrika_proj_1/scripts/pip/DESeq3_ZE_mRNA_p02.R'
         mRNA_p02_ZE_design='/proj/uppstore2017145/V2/Ulrika_proj_1/scripts/pip/ZE_DesignA2_mRNA_p02.csv'
         $Script_mRNA_p02_ZE $ListCountZE'.table' $script_loc/$DGE_design
         #echo -e '\e[0;32mCount files list mRNA_p02 ZE: \e[0m' $ListCountZE
         echo -e '\e[0;35mScript for mRNA_p02 ZE: \e[0m' $Script_mRNA_p02_ZE
         #echo -e '\e[0;32mInput count table ZE: \e[0m' $ListCountZE'.table'
         #echo -e '\e[0;32mDesign file for mRNA_p02 ZE: \e[0m' $mRNA_p02_ZE_design
         echo -e '\e[0;33mResult of ZE Log2FC mRNA for target analysis: \e[0m' '/proj/uppstore2017145/V2/Ulrika_proj_1/scripts/DGE/ZE_Log2FC_mRNA_p02.csv'
   fi
    echo -e '\e[0;35mScript for miRNA DGE: \e[0m' $script_loc/$DGE_script
echo -e '\e[0;33m-------------------------------------------------------------------------\e[0m'
fi
#############################
############################# 08Nov2024


#############################
#Step 4: Pre DGE analysis 
#############################
Input_count_cut=$P_prediction_list'.id.txt'
count_table2=$mature_loc/$Input_count_cut'countTable.txt'
if [ $Step4 = 'yes' ]
   then
   echo -e '\e[0;33mStep4: Plant criteria based count file extraction ------------------\e[0m'
    #$script_loc/combine_htseq_counts.pl $count_list $count_table
    echo -e '\e[0;36mCount table full: \e[0m' $count_table 
    $script_loc/plot_mature_analysis4.py \
        $mature_loc/Picea_mature.fa \
        $mature_loc/Picea_mature_known_id.txt \
        $mature_loc/Picea_mature_variant_id.txt \
        $mature_loc/$P_prediction_list \
        $mature_loc/Picea_mature_blastn_parsed.txt
    echo -e '\e[0;36mScript used for plant criteria ID generation and analysis: \e[0m' $script_loc/plot_mature_analysis4.py 
    #echo -e '\e[0;36mPlant criteria files list: \e[0m' $mature_loc/$P_prediction_list 
    Input_count_cut=$P_prediction_list'.id.txt'
    $script_loc/cut_countRable.py \
              $count_table $mature_loc/$Input_count_cut
    echo -e '\e[0;36mScript to cut count table: \e[0m' $script_loc/cut_countRable.py
    echo -e '\e[0;36mmiRNA IDs with p criteria: \e[0m' $mature_loc/$Input_count_cut
    count_table2=$mature_loc/$Input_count_cut'countTable.txt'
    #echo -e '\e[0;36mScript used to cut count table for plant criteria: \e[0m' $script_loc/cut_countRable.py
    echo -e '\e[0;36mCount Table with plant criteria: \e[0m' $count_table2
    echo -e '\e[0;33------------------------------------------------------------------- ------------------\e[0m'
fi
##############################


##############################
#DGE for Count table with plant criteria
##############################
if [ $DGE_miRNA = 'yes' ]
   then
   echo -e '\e[0;33Section: DGE for miRNA with plant criteria -------------------------------------------\e[0m'
   $script_loc/$DGE_script $count_table2 $script_loc/$DGE_design
   echo -e '\e[0;36mCount data files list: \e[0m' $count_list
   echo -e '\e[0;36mScript used for DGE of miRNA with plant criteria: \e[0m' $script_loc/$DGE_script
   echo -e '\e[0;36mCount table input for script: \e[0m' $count_table2
   echo -e '\e[0;36mDesign file input for script: \e[0m' $script_loc/$DGE_design
   echo -e '\e[0;33m------------------------------------------------------------------- ------------------\e[0m'  
fi
######################################


######################################
if [ $Step5 = 'yes' ]
   then
   echo -e '\e[0;33mSection: Step5, Post DGE analysis and files generation -------------------------------\e[0m'
    $script_loc/extract_miRNA.py \
       $mature_loc/Picea_mature.fa $DGE_gene_list # FASTA for DGE List
    #echo -e '\e[0;36mScript to generate Input FASTA file for Target prediction: \e[0m' $script_loc/extract_miRNA.py 
    ################################## 
    $script_loc/2_mature_analysis.py \
       $mature_loc/Picea_mature.fa \
       $mature_loc/Picea_mature_known_id.txt \
       $mature_loc/Picea_mature_variant_id.txt \
       $mature_loc/$P_prediction_list \
       $mature_loc/Picea_mature_blastn_parsed.txt $DGE_gene_list'.fa'
    echo -e '\e[0;36mScript for post DGE analysis: \e[0m' $script_loc/2_mature_analysis.py
    echo -e '\e[0;33mSubsection, Input DGE FASTA file generation for Target prediction------------------\e[0m'
    echo -e '\e[0;36mScript to generate Input FASTA file for Target prediction: \e[0m' $script_loc/extract_miRNA.py
    $script_loc/extract_miRNA.py \
       $mature_loc/Picea_mature.fa $DGE_gene_list # FASTA for DGE List 
echo -e '\e[0;33m------------------------------------------------------------------- ------------------\e[0m'
#echo '-----------pre DGE steps checkouts----------------------------------------'
#echo '(A) Combined ZE+SE list for  miRNA annotation: ' $predict
#echo '(B) Plant criteria list from pipeline: ' $mature_loc/$P_prediction_list
#echo '-----------cross check and make sure that experiment is correct-----------'
#echo '(1) Fasta file list used for count table: ' $fastq_list2
#echo '(2) Input list of count file to make combined count table: ' $count_list
#echo '(3) Original count table genetated at: ' $count_table
#echo '(4) Trunicated count table with P criteria, used for DEseq2: '$count_table2
#echo '(5) DGE design file for DESeq2: ' $script_loc/$DGE_design
#echo '(6) DGE script used: ' $script_loc/$DGE_script
#echo '(7) P crietria passed list: '  $mature_loc/$P_prediction_list
#echo '--------------------------------------------------------------------------'
fi
######################################


######################################
#MIRNA BASE analysis 
######################################
if [ $Step6 = 'yes' ]
   then
      echo -e '\e[0;33mSection: Step6, ---------------------------------------------------------------\e[0m'
   echo -e '\e[0;35mmirBASEID generation and  Ana Alvis review input generation \e[0m'
      #Making FASTA files ################################################
      #$script_loc/extract_miRNA.py $mature_loc/Picea_mature.fa \
       #       $mature_loc/Picea_mature_known_id.txt #MAKE FASTA file of known one
      #$script_loc/extract_miRNA.py $mature_loc/Picea_mature.fa \
       #       $mature_loc/Picea_mature_variant_id.txt #MAKE FASTA file of Variant oni
      #$script_loc/extract_miRNA.py $mature_loc/Picea_mature.fa \
       #       $DGE_gene_list
      #grep -h mature_seq $Res_fold22/P*/*_predictions > $Res_fold22/predicted_seq.txt
      #######################MIR BASE analysis ##########################
      $script_loc/miR_countsOriginal.py \
                 $mature_loc/Picea_mature_blastn.txt \
                 $mature_loc/Picea_mature_known_id.txt
      echo -e '\e[0;32mScript for MirBASE count generation \e[0m' $script_loc/miR_countsOriginal.py
      #echo -e '\e[0;33mInputs for MirBASE count generation \e[0m' $mature_loc/Picea_mature_blastn.txt $mature_loc/Picea_mature_known_id.txt
      #echo -e '\e[0;35mOutput KnownID Vs MirBASE ID (Important) \e[0m' /proj/uppstore2017145/V2/Ulrika_proj_1/results/miRNA_predicted/Picea_mature_known_id.txt.MiRbase.O.csv
      #echo -e '\e[0;33m------------------------------------------------------------------------------\e[0m'
      # Extract known miRBase Id expressed, plant criteria
      $script_loc/miR_counts4.py \
                 $mature_loc/Picea_mature_blastn.txt \
                 $mature_loc/$P_prediction_list*known.id.txt
      echo -e '\e[0;32mScript to Extract known miRBase Id expressed, plant criteria \e[0m' $script_loc/miR_counts4.py
      #echo -e '\e[0;33mInputs for miR_counts4.py \e[0m' $mature_loc/Picea_mature_blastn.txt $mature_loc/$P_prediction_list*known.id.txt
      #echo -e '\e[0;35mOutput KnownID Vs MirBASE ID (Plant Criteria) \e[0m' $mature_loc/$P_prediction_list*known.id.txt.MiRbase.csv
      #echo -e '\e[0;33m------------------------------------------------------------------------------\e[0m'
      $script_loc/Merge_miRBase2.py $mature_loc/SE_p_prediction_list.xyz.known.id.txt.MiRbase.csv \
                 $mature_loc/ZE_p_prediction_list.xyz.known.id.txt.MiRbase.csv \
                 $mature_loc/SE_ZE_MiRBase_freq.csv
      #echo -e '\e[0;35mScript to check MirBASE Id Vs known miRNA, its used for checking \e[0m' $script_loc/Merge_miRBase2.py 
      #echo -e '\e[0;33m------------------------------------------------------------------------------\e[0m'
      echo -e '\e[0;32mScript to calculate mirBASE frequency and Ana ALvis review litrature miRBASE related miRNA \e[0m' $script_loc/Merge_miRBaseOriginal.py
      #echo -e '\e[0;32mInputs for Merge_miRBaseOriginal.py \e[0m' $mature_loc/Picea_mature_known_id.txt.MiRbase.O.csv 
      $script_loc/Merge_miRBaseOriginal.py \
                 $mature_loc/Picea_mature_known_id.txt.MiRbase.O.csv \
                 $mature_loc/Picea_mature_known_id.txt.MiRbase.O.csv \
                 $mature_loc/SE_ZE_MiRBase_freq.O.csv
      #echo -e '\e[1;35mOutput of Merge_miRBaseOriginal.py \e[0m' $mature_loc/SE_ZE_MiRBase_freq.O.csv 
      #echo -e '\e[1;32mOutput of Merge_miRBaseOriginal.py \e[0m' /proj/uppstore2017145/V2/Ulrika_proj_1/scripts/heatmap/Ana_alvis_list2.csv
      #echo -e '\e[0;33m------------------------------------------------------------------------------\e[0m'
      mirBASE_full_list='/proj/uppstore2017145/V2/Ulrika_proj_1/results/miRNA_predicted/Picea_mature_known_id.txt.MiRbase.O.csv'
      #mirBASE_full_list=$mature_loc/'Picea_mature_known_id.txt.MiRbase.O.csv'
      echo -e '\e[1;35mmirBASE_full_list used for heatmap annotation\e[0m' $mirBASE_full_list
      echo -e '\e[0;33m------------------------------------------------------------------------------\e[0m'
fi
########################################




########################################
# Section for combined analysis of SE and ZE
########################################
if [ $Comb = 'yes' ]
   then
   DGE_loc='/proj/uppstore2017145/V2/Ulrika_proj_1/scripts/DGE'
   DE_SE='DGE_SE.txtmiR.csv'
   DE_ZE='DGE_ZE.txtmiR.csv'
   cat $DGE_loc/$DE_SE $DGE_loc/$DE_ZE > $DGE_loc/'DE_ZE_SE_comb.csv'
   DE_SE_Total='DGE_SE.txtmiR_Total.csv'
   DE_ZE_Total='DGE_ZE.txtmiR_Total.csv'
   cat $DGE_loc/$DE_SE_Total $DGE_loc/$DE_ZE_Total > $DGE_loc/'DGE_SE_ZE.txtmiR_Total.csv'
   $script_loc/combAnalysis.py \
              $DGE_loc/'DE_ZE_SE_comb.csv' \
              $DGE_loc/'DGE_SE_ZE.txtmiR_Total.csv' \
              $DGE_loc/'DGE_SE.txtmiR_Total_count.csv' \
              $DGE_loc/'DGE_ZE.txtmiR_Total_count.csv' 
fi
#########################################



#########################################
#MRNA (OLD) section
#########################################
#target_folder='/proj/uppstore2017145/V2/Ulrika_proj_1/scripts/Target_predict'
#if [ $mRNA = 'yes' ]
 #  then
  # echo 'mRNA analysis ---------------'
   #$script_loc/target1.R $target_folder/dds_genes_TEs.rda
   #echo $script_loc/target1.R
#fi
#########################################
#########################################



#########################################
#Target analysis section
#########################################
if [ $Target = 'yes' ]
   then
   echo -e '\e[0;33mSection: Target analysis-----------------------------------------------------\e[0m'
   echo -e '\e[0;31mTarget analysis section, 1) Run DGE_miRNA section and perform FASTA file generation from DGE gene list of miRNA, 2) Take VST files from DGE section of miRNA and mRNA, 3) Use target file given by NAT, Upload on Target analysis website https://www.zhaolab.org/psRNATarget/, 5) Results are uploaded on UPPMAX and then run this section \e[0m'
   
 if [ $1 = 'SE' ]
   then
   echo -e '\e[0;35mScript for SE miRNA DGE gene list and VST file generation: \e[0m' $script_loc/'DESeq5_SE.R'
   echo -e '\e[0;35mScript for SE mRNA DGE gene list and VST file generation: \e[0m' $script_loc/'DESeq4_SE_mRNA_p02.R'
   $script_loc/extract_miRNA.py \
           $mature_loc/Picea_mature.fa $DGE_gene_list #FASTA DGE
   VST_SE_mRNA=$DGE/'mRNA_DGE_SE_VST.txt'
   VST_SE_miRNA=$DGE/'DGE_SE_VST.txt'
   Corr_pair=$DGE/'VST_SE_CC_negative_pair.csv'
   Target_Predict_result='SE_Pab02_psRNATargetJob-1730223103211742.txt2'
   echo -e '\e[0;32mInput1: Target file \e[0m' $target_folder/$Target_Predict_result
   echo -e '\e[0;32mInput2: VST_SE_miRNA \e[0m' $VST_SE_miRNA
   echo -e '\e[0;32mInput3: VST_SE_mRNA \e[0m' $VST_SE_mRNA
   echo -e '\e[0;35mScript for VST correlation analysis: \e[0m' $script_loc/VST_correlation_check4.py
   #$script_loc/VST_correlation_check4.py \
    #       $target_folder/$Target_Predict_result \
     #      $VST_SE_miRNA \
      #     $VST_SE_mRNA  \
       #    $Corr_pair
   echo -e '\e[0;33mOutput1: miRNA-mRNA pairs passed the test: \e[0m' $Corr_pair
   mirBASE_full_list='/proj/uppstore2017145/V2/Ulrika_proj_1/results/miRNA_predicted/Picea_mature_known_id.txt.MiRbase.O.csv' # Generated from Step6
   mirBASE_full_Novel=$script_loc/'VST_CCpair_miRNA.py'
   echo -e '\e[0;35mScript for heatmap generation: \e[0m' $mirBASE_full_Novel    
 #  $mirBASE_full_Novel \
  #      $mirBASE_full_list \
   #     $VST_SE_miRNA \
    #    $Corr_pair
   echo -e '\e[0;33mOutput2: Heatmap of -ve CC: \e[0m' '/home/ahsan/SE_CC_miRNA.pdf'
   Log_file=$DGE/'VST_SE_CC_negative_pair.csv.log'   
   echo -e '\e[0;33mOutput3: Log file for CC analysis: \e[0m' $Log_file
   Log_analysis=$script_loc/'VST_CC_Log_analysis.py'
   echo -e '\e[0;35mScript for Log analysis: \e[0m' $Log_analysis
   $Log_analysis $Log_file       
 fi

 


   #SE_mRNA_FCtable_p02=$DGE/'SE_Log2FC_mRNA_p02.csv'
   #SE_miRNA_FCtable_p02=$DGE/'SE_Log2FC_miRNA.csv'
   if [ $1 = 'ZE' ]
   then
   echo -e '\e[0;35mScript for ZE miRNA DGE gene list and VST file generation: \e[0m' $script_loc/'DESeq3_ZE.R'
   echo -e '\e[0;35mScript for ZE mRNA DGE gene list and VST file generation: \e[0m' $script_loc/'DESeq3_ZE_mRNA_p02.R'
   $script_loc/extract_miRNA.py \
           $mature_loc/Picea_mature.fa $DGE_gene_list #FASTA DGE
   Target_Predict_result='ZE_psRNATargetJob-1733149339678709.txt'
   echo -e '\e[0;32mInput1: Target file \e[0m' $target_folder/$Target_Predict_result 
   VST_ZE_miRNA=$DGE/'DGE_ZE_VST.txt'
   VST_ZE_mRNA=$DGE/'mRNA_DGE_ZE_VST.txt'
   Corr_pair=$DGE/'VST_ZE_CC_negative_pair.csv'
   echo -e '\e[0;32mInput2: VST_ZE_miRNA \e[0m' $VST_ZE_miRNA 
   echo -e '\e[0;32mInput3: VST_ZE_mRNA \e[0m' $VST_ZE_mRNA
   echo -e '\e[0;35mScript for VST corr analysis: \e[0m' $script_loc/VST_correlation_check_ZE.py
  # $script_loc/VST_correlation_check_ZE.py \
   #        $target_folder/$Target_Predict_result \
    #       $VST_ZE_miRNA \
     #      $VST_ZE_mRNA \
      #     $Corr_pair
   mirBASE_full_list='/proj/uppstore2017145/V2/Ulrika_proj_1/results/miRNA_predicted/Picea_mature_known_id.txt.MiRbase.O.csv' # Generated from Step6
   mirBASE_full_Novel=$script_loc/'VST_CCpair_miRNA_ZE.py'
   echo -e '\e[0;35mScript for heatmap generation: \e[0m' $mirBASE_full_Novel 
   fi

fi
#########################################








#########################################28/11/2024
#HeatMap:: This section is for heatmap generation
#########################################
if [ $Heatmap = 'yes' ]
   then
   echo -e '\e[0;33mSection: Heatmap SE and ZE ---------------------------------------------------\e[0m'
   mirBASE_full_list='/proj/uppstore2017145/V2/Ulrika_proj_1/results/miRNA_predicted/Picea_mature_known_id.txt.MiRbase.O.csv' # Generated from Step6
   mirBASE_full_Novel=$script_loc/'Ulrika_plot_SEz_full_novel.py'
   mirBASE_full_Novel_ZE=$script_loc/'Ulrika_plot_ZEz_full_novel.py'
   echo -e '\e[3;32mScript for mirBASE_full_Novel_SE:\e[0m' $mirBASE_full_Novel
   echo -e '\e[3;32mScript for mirBASE_full_Novel_ZE:\e[0m' $mirBASE_full_Novel_ZE
   $mirBASE_full_Novel $mirBASE_full_list $DGE/DGE_SE_VST.txt 
   $mirBASE_full_Novel_ZE $mirBASE_full_list $DGE/DGE_ZE_VST.txt 
   SE_DF=$DGE/'DGE_SE_VST.txt.SE.csv'
   ZE_DF=$DGE/'DGE_ZE_VST.txt.ZE.csv'
   Combined_heatmap=$script_loc/'Ulrika_heatmap_combined2.py'
   echo -e '\e[3;32mScript for combined heatmap SE+ZE:\e[0m' $Combined_heatmap
   $Combined_heatmap $mirBASE_full_list $SE_DF $ZE_DF
   echo -e  '\e[1;33mOutput1 Heatmaps saved at:\e[0m' '/home/ahsan/*.pdf'
   VST_combined='/proj/uppstore2017145/V2/Ulrika_proj_1/results/miRNA_predicted/Picea_mature_known_id.txt.MiRbase.O.csv.DF.ZE_SE.csv'
   Split_heatmap=$script_loc/'Ulrika_heatmap_split.py'
   echo -e '\e[3;32mScript to split heatmap:\e[0m' $Split_heatmap
   $Split_heatmap $VST_combined
   echo -e '\e[0;33m--------------------------------------------------------------------------\e[0m'
fi

if [ $Heatmap_mRNA = 'yes' ]
   then
   if [ $1 = 'SE' ]
      then
      Script_mRNA_DGE_SE=$script_loc/'DESeq4_SE_mRNA_p02.R'
      echo -e '\e[3;32mScript for Script_mRNA_DGE_SE:\e[0m' $script_loc/'DESeq4_SE_mRNA_p02.R'
      mRNA_VST_SE=$DGE/'mRNA_DGE_SE_VST.txt'     
      echo -e '\e[3;31mInput1: mRNA_VST_SE:\e[0m' $DGE/'mRNA_DGE_SE_VST.txt'
      mRNA_Heatmap_SE=$script_loc/'mRNA_Heatmap_SE.py'
      echo -e '\e[3;32mScript for mRNA_Heatmap_SE:\e[0m' $mRNA_Heatmap_SE
      VST_mRNA_SE=$DGE/'mRNA_DGE_SE_VST.txt.mRNA.SE.csv'
      #$mRNA_Heatmap_SE $DGE/'mRNA_DGE_SE_VST.txt'
      echo -e '\e[3;33mOutput1: VST_mRNA_SE:\e[0m' $VST_mRNA_SE
   fi



   if [ $1 = 'ZE' ]
      then
      Script_mRNA_DGE_ZE=$script_loc/'DESeq3_ZE_mRNA_p02.R' 
      echo -e '\e[3;32mScript for Script_mRNA_DGE_ZE:\e[0m' $Script_mRNA_DGE_ZE
      mRNA_VST_ZE=$DGE/'mRNA_DGE_ZE_VST.txt' 
      echo -e '\e[3;31mInput1: mRNA_VST_ZE:\e[0m' $mRNA_VST_ZE
      mRNA_Heatmap_ZE=$script_loc/'Ulrika_plot_ZEz_mRNA.py'    
      echo -e '\e[3;32mScript for mRNA_Heatmap_ZE:\e[0m' $mRNA_Heatmap_ZE  
      $mRNA_Heatmap_ZE $mRNA_VST_ZE $mRNA_VST_ZE      
   fi


   if [ $1 = 'CMB' ]
      then
      VST_mRNA_SE=$DGE/'mRNA_DGE_SE_VST.txt.mRNA.SE.csv'
      #VST_mRNA_ZE=$DGE/'mRNA_DGE_ZE_VST.txt.ZE.mRNA.csv' 
       VST_mRNA_ZE=$DGE/'mRNA_DGE_ZE_VST.txt.ZE.mRNA.csv'
      echo -e '\e[3;31mInput1: VST_mRNA_SE:\e[0m' $VST_mRNA_SE
      echo -e '\e[3;31mInput1: VST_mRNA_ZE:\e[0m' $VST_mRNA_ZE
      Combined_heatmap=$script_loc/'Ulrika_heatmap_comb_mRNA.py'
      echo -e '\e[3;32mScript for comb heatmap SE+ZE(mRNA):\e[0m' $Combined_heatmap 
      $Combined_heatmap $VST_mRNA_SE $VST_mRNA_SE $VST_mRNA_ZE     
      
   fi

fi



#########################################
