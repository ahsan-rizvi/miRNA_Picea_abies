#!/sw/apps/conda/latest/rackham_stage/bin/python3
#- Dr Ahsan Z Rizvi
#-------------------------
import numpy
import sys
import os
import numpy as np
#import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import mstats
import pandas as pd
#-------------------------


miRNA_exp=list()
with open(sys.argv[2]) as file:
    for line in file:
        #print(line)
        miRNA_exp.append(line)
        
#sys.exit()
mRNA_exp=list()
with open(sys.argv[3]) as file:
    for line in file:
#        print(line)
        mRNA_exp.append(line)


#sys.exit()


#Select miRNA and mRNA significant pairs(E>4)
miRNA_pair=list()
mRNA_pair=list()
E_value=list()
i=0
total=0
total_le4=0
break1=0


with open(sys.argv[1]) as file:
    for line in file:
        #############
        #print(line)
        #break1=break1+1
        #if break1>5000:
         #  break
        #############
        temp1=(line.split('\t'))
        if len(temp1)>1 and i>1:
           total=total+1
           if float(temp1[2]) > 4:
              total_le4=total_le4+1   
              miRNA_pair.append(temp1[0])
              mRNA_pair.append(temp1[1])
              E_value.append(temp1[2])
        i=i+1
print('Total pairs present:',total)
print('Total pairs Evalue> 4:',total_le4)




#sys.exit()
############################################
#correlation check
P_value_lim=0.05
temp22=sys.argv[4]+'.ZE.txt'
f  = open(temp22, "w")
print('list of miRNA-mRNA pair passed the CC test: ',temp22)
#sys.exit()
#temp33=sys.argv[4]+'.log'
#f2  = open(temp33, "w")
Count_CC2=0
Count_CC=0
CC_Pvalue_miRNA=list()
CC_Pvalue_mRNA=list()
break1=0
for i in range(0,len(miRNA_pair)):
    break1=break1+1
    ###############
    #if break1>50:
     #  break
    ###############
    temp1=miRNA_pair[i]
    temp2=mRNA_pair[i]
    #print(temp1,temp2)
    miRNAid_in_Target=''
    mRNAid_in_Target=''
    for j in range(0,len(miRNA_exp)):    
        temp3=miRNA_exp[j].split('\t')
        temp4=temp3[0].split(',')
        miRNAid_in_Target=temp4[0]
        if miRNAid_in_Target==temp1:
           miRNAid_in_Target_matched=miRNA_exp[j]
           #print(miRNAid_in_Target_matched)
    for k in range(0,len(mRNA_exp)):
        temp5=mRNA_exp[k].split('\t')
        mRNAid_in_Target=temp5[0]
        if temp2==mRNAid_in_Target: 
           mRNAid_in_Target_matched=mRNA_exp[k]
           #print(mRNA_exp[k])
    #print('---------------------------------------------------------------------')
    miRNA_temp=miRNAid_in_Target_matched.split('\n')
    miRNA_temp2=miRNA_temp[0].split('\t')

    miRNA_VST=[float(miRNA_temp2[10]),float(miRNA_temp2[12]),float(miRNA_temp2[14]),float(miRNA_temp2[15]),float(miRNA_temp2[18]),float(miRNA_temp2[20]), float(miRNA_temp2[23]),float(miRNA_temp2[25]),float(miRNA_temp2[27]),float(miRNA_temp2[29]),float(miRNA_temp2[30]),float(miRNA_temp2[32]),  float(miRNA_temp2[34]),float(miRNA_temp2[36]),float(miRNA_temp2[38]),float(miRNA_temp2[40]),float(miRNA_temp2[-6]),float(miRNA_temp2[-7]), float(miRNA_temp2[-9]),float(miRNA_temp2[-10]),float(miRNA_temp2[-12]),float(miRNA_temp2[-13])]

    #print('##########################')
    mRNA_temp=mRNAid_in_Target_matched.split('\n')     
    mRNA_temp2=mRNA_temp[0].split('\t')

    mRNA_VST=[float(mRNA_temp2[10]),float(mRNA_temp2[12]),float(mRNA_temp2[14]),float(mRNA_temp2[15]),float(mRNA_temp2[18]),float(mRNA_temp2[20]), float(mRNA_temp2[23]),float(mRNA_temp2[25]),float(mRNA_temp2[27]),float(mRNA_temp2[29]),float(mRNA_temp2[30]),float(mRNA_temp2[32]),  float(mRNA_temp2[34]),float(mRNA_temp2[36]),float(mRNA_temp2[38]),float(mRNA_temp2[40]),float(mRNA_temp2[-6]),float(mRNA_temp2[-7]), float(mRNA_temp2[-9]),float(mRNA_temp2[-10]),float(mRNA_temp2[-12]),float(mRNA_temp2[-13])]

    #-------------------------------------------------
    [Corr,P_value]=mstats.pearsonr(mRNA_VST,miRNA_VST)
    #print(Corr,P_value)
    P_value_lim=0.05
    if Corr<0:
       Count_CC2=Count_CC2+1
    if Corr<0 and P_value < P_value_lim:
       #print(miRNA_temp2[0],mRNA_temp2[0],E_value[i],Corr,P_value)
       #f.writelines(miRNA_temp2[0]+'\t'+mRNA_temp2[0]+'\n')
       f.writelines(miRNA_temp2[0]+'\t'+mRNA_temp2[0]+'\t'+E_value[i]+'\t'+str(Corr)+'\t'+str(P_value)+'\n')
       #f2.writelines(miRNA_temp2[0]+'\t'+mRNA_temp2[0]+'\t'+E_value[i]+'\t'+str(Corr)+'\t'+str(P_value)+'\n')
       CC_Pvalue_miRNA.append(miRNA_temp2[0])
       CC_Pvalue_mRNA.append(mRNA_temp2[0])
       Count_CC=Count_CC+1
    #-------------------------------------------------
       
f.close()
########################################
#print('Mirna-mRNA pairs in negative CC values: ',Count_CC2) 
#print('P value cut off: ', P_value_lim)        
#print('Mirna-mRNA pairs in negative CC values and passed P value test: ',Count_CC)   













