#!/sw/apps/conda/latest/rackham_stage/bin/python3
# Dr Ahsan Z Rizvi
#########################
import numpy
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
#########################
#./plot_mature_analysis.py Picea_mature.fa Picea_mature_known_id.txt Picea_mature_variant_id.txt P_prediction_list2.txt Picea_mature_blastn_parsed.txt

#########################

######################################
#Check Picea_mature.fa
######################################
seq1=list()
id1=list()
with open(sys.argv[1]) as file:
    for line in file:
        #print(line)
        temp=line.split('\t')
        temp2=temp[1]
        lenTemp=len(temp2)
        temp3=temp2[0:lenTemp-1]
        seq1.append(temp3)
        #print(temp3,len(temp3))
P_uniq3,P_counts3 = numpy.unique(seq1,return_counts=True)
print('Non redundent P seqs: ',len(P_uniq3))


temp2=sys.argv[1]+'.csv'
f = open(temp2, "w")
for i in range(0, len(P_uniq3)):
    #print(i,P_uniq3[i])
    f.writelines(P_uniq3[i]+'\n')
f.close()

print('NonRedundent Pridicted seq saved at: ',temp2)
#print('----------------------------------------------------------------------------')

temp2=sys.argv[2]+'.Seq.csv'
f = open(temp2, "w")
i=0
seqKnown=list()
with open(sys.argv[2]) as file:
    for line in file:
        i=i+1
        lenN=len(line)
        line2=line[0:(lenN-1)]
        lenN2=len(line2)
        #print(lenN2,line2)     
        if (i%2)==0: 
            #print(line2,lenN2)
            f.writelines(line2+'\n')
            seqKnown.append(line2)
f.close()
print('Total Known: ',len(seqKnown))
print('Known seq saved at: ',temp2)


temp2=sys.argv[3]+'.Seq.csv'
f = open(temp2, "w")
i=0
seqVariant=list()
with open(sys.argv[3]) as file:
    for line in file:
        i=i+1
        lenN=len(line)
        line2=line[0:(lenN-1)]
        lenN2=len(line2)
        if (i%2)==0:
            #print(line2)
            f.writelines(line2+'\n')
            seqVariant.append(line2)
f.close()
print('Total Variant: ',len(seqVariant))
print('Variant seq saved at: ',temp2)




k=0
n=0
for i in range(0, len(P_uniq3)):
    for j in  range(0, len(seqKnown)):
        if P_uniq3[i]==seqKnown[j]:
           k=k+1
    for m in  range(0, len(seqVariant)):    
        if P_uniq3[i]==seqVariant[m]:
           n=n+1


#print('NonRedundent Predicted sequence: ',len(P_uniq3))
#print('Known predicted, matched with ZE+SE predicted: ',k)
#print('Known Variant, matched with ZE+SE predicted: ',n)

temp=len(P_uniq3)-(k+n)
#print('Novel predicted: ',temp)




