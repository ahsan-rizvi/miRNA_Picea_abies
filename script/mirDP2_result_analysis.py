#!/sw/apps/conda/latest/rackham_stage/bin/python3
# Dr Ahsan Z Rizvi
#########################
import numpy
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
#########################
#id1 is Picea_mature.fa's ID
######################################
#Check Picea_mature.fa
######################################
seq1=list()
id1=list()
with open(sys.argv[1]) as file:
    for line in file:
        if line[0] == '>':
           line2=line.split(';')
           line3=line2[0].split(' ')
           line4=(line3[0].split('>'))
           id1.append(line4[1])
           #print(line4[1])
        if line[0] != '>':
           line5=line.split('\n')
           line6=line5[0]
           seq1.append(line6)
           #print(line6)
######################################
#OUT:: id1 contains all ids, seq1 contains all seq
######################################

        
######################################
# SCREEN OUT all bigger sequences
######################################
id2=list()
seq2=list()
id3=list()
seq3=list()
for i in range(0, len(seq1)):
    len_seq1=len(seq1[i])
    if (len_seq1>17) and (len_seq1<25):
       #print(id1[i],seq1[i],len_seq1)
       id2.append(id1[i])
       seq2.append(seq1[i])
#####################################


#####################################
# READ Picea_mature_known_id.txt for Known IDs
#####################################
#print(len(id2))
#print(len(seq2))
known_id=list()
varient_id=list()
with open(sys.argv[2]) as file:
    for line in file:
        line3=line.split('\n')
        #print(line3[0])
        known_id.append(line3[0])
####################################
print('Known_id file: ',sys.argv[2]) 

####################################
# READ  Picea_mature_variant_id.txt for varient ID
####################################
with open(sys.argv[3]) as file:
    for line in file:
        line3=line.split('\n')
        #print(line3[0])
        varient_id.append(line3[0])
###################################
print('varient_id file: ',sys.argv[3])

##################################
# READ Picea_mature_blastn_parsed.txt for parsed file
##################################
parse_id=list()
parse_match=list()
with open(sys.argv[5]) as file:
    for line in file:
        temp1=(line.split('\t'))
        #print(line)
        #print(temp1[1],temp1[5])
        parse_id.append(temp1[1])
        parse_match.append(temp1[5]) 
#print(len( parse_id))
#print(len( parse_match))
print('parse_id_id file: ',sys.argv[5])
#exit()
##################################

print('Predicted: ',len(id2))
print('Known: ',len(known_id))
print('Varients: ',len(varient_id))

#sys.exit()

##################################
# MATCH and ADD Known and Vrient
##################################
#print('-------------Known-----------------------------------------')
known_id2=list()
known_id_seq=list()
known_id_seq_len=list()
known_id_seq_fist_base=list()
for i in range(0, len(id2)):
    for j in  range(0, len(known_id)):
        if id2[i] == known_id[j]:
           #print(id2[i],known_id[j],seq2[i],len(seq2[i]),'Known',seq2[i][0])
           id3.append(id2[i])
           known_id2.append(id2[i])
           known_id_seq.append(seq2[i])
           known_id_seq_len.append(len(seq2[i]))
           known_id_seq_fist_base.append(seq2[i][0])
#exit()
#os.system('clear')
known_seq_len_nt, known_seq_len_nt_freq = numpy.unique(known_id_seq_len,return_counts=True)
#print('known_seq_len_nt',known_seq_len_nt)
#print('known_seq_len_nt_freq', known_seq_len_nt_freq)

known_seq_len_nt2, known_seq_len_nt_freq2 = numpy.unique(known_id_seq_fist_base,return_counts=True)
#print(known_seq_len_nt2)
#print(known_seq_len_nt_freq2)


#print('-------------Varient-----------------------------------------')

varient_id2=list()
varient_id_seq=list()
varient_id_seq_len=list()
varient_id_seq_first_base=list()
for i in range(0, len(id2)):
    for j in  range(0, len(varient_id)):
        if id2[i] == varient_id[j]:
#           print(id2[i],seq2[i],len(seq2[i]),'Varient',seq2[i][0])
           id3.append(id2[i])
           varient_id2.append(id2[i])
           varient_id_seq.append(seq2[i])
           varient_id_seq_len.append(len(seq2[i]))
           varient_id_seq_first_base.append(seq2[i][0])
varient_seq_len_nt, varient_seq_len_nt_freq = numpy.unique(varient_id_seq_len,return_counts=True)
#print('varient_seq_len_nt',varient_seq_len_nt)
#print('varient_seq_len_nt_freq', varient_seq_len_nt_freq)


varient_seq_len_nt2, varient_seq_len_nt_freq2 = numpy.unique(varient_id_seq_first_base,return_counts=True)
#print(varient_seq_len_nt2)
#print(varient_seq_len_nt_freq2)
#print('Known + Varient: ',len(id3))

#######################################
#NOVEL Analysis
#######################################
#print('--------------Novel-----------------')
Novel_id=list()
Novel_seq=list()
Novel_seq_len=list()
Novel_seq_first_base=list()
for i in range(0, len(seq2)):
    count=0
    for j in  range(0, len(id3)):#id3 has Known+varient ids
        if id2[i] == id3[j]:
           count=count+1
    if count==0:
       #print(count)
       #print(id2[i],seq2[i],len(seq2[i]),'Novel')
       Novel_id.append(id2[i])
       Novel_seq.append(seq2[i])
       Novel_seq_len.append(len(seq2[i]))
       Novel_seq_first_base.append(seq2[i][0])
       #print(seq2[i],len(seq2[i]),seq2[i][0])
print('Novel',len(Novel_id))
#---------------------------------------
Novel_seq_len_nt, Novel_seq_len_nt_freq = numpy.unique(Novel_seq_len,return_counts=True)
#print('Novel_seq_len_nt',Novel_seq_len_nt)
#print('Novel_seq_len_nt_freq', Novel_seq_len_nt_freq)


Novel_seq_len_nt2, Novel_seq_len_nt_freq2 = numpy.unique(Novel_seq_first_base,return_counts=True)
#print(Novel_seq_len_nt2)
#print(Novel_seq_len_nt_freq2)
########################################




########################################
# PLANT CRITERIA
# Checking the file P_prediction_list2.txt
########################################
#print(sys.argv[4])
#print('-------------P_CRITERIA----------------------------------------------')
P_list=list()
P_list_len=list()
#P_list_nuc=list()
count=0
#print(sys.argv[4])
count2=1
print('P-DGE list: ',sys.argv[6])

with open(sys.argv[6]) as file:
    for line in file:
        if count2 % 2:
           pass
        else:
           temp=line.rstrip()
           #print(temp,len(temp))
           P_list.append(temp)
           P_list_len.append(len(temp))
        count2=count2+1


P_uniq,P_counts = numpy.unique(P_list,return_counts=True)
#print('Total_P_list: ',len(P_list))
#print('Total_P_Uniq: ',len(P_uniq))

#sys.exit()

P_uniq_len,P_counts_len = numpy.unique(P_list_len,return_counts=True)

P_uniq_len=list()
P_uniq_fist_base=list()
for i in range(0, len(P_uniq)):
    #print(i,P_uniq[i],len(P_uniq[i]),P_uniq[i][0])
    P_uniq_len.append(len(P_uniq[i]))
    P_uniq_fist_base.append(P_uniq[i][0])
#print(len(P_uniq_len),len(P_uniq_fist_base))
P_uniq2,P_counts2 = numpy.unique(P_uniq_len,return_counts=True)
#print(P_uniq2)
#print(P_counts2)
P_uniq3,P_counts3 = numpy.unique(P_uniq_fist_base,return_counts=True)
#print(P_uniq3)
#print(P_counts3)
#print('----------------------------------------------------------')



#####################MAIN LOGIC##################
#print(sys.argv[4])

#temp22=sys.argv[4]+'.id.txt2'
#print(temp22)
#f = open(temp22, "w")


Tag_P_criteria=0
Tag_P_criteria2=0
Novel_P_criteria=0
DE_P_known_len=list()
DE_P_known_seq=list()
DE_P_varient_seq=list()
DE_P_novel_seq=list()
DE_P_varient_len=list()
DE_P_novel_len=list()


miR_known=list()
miR_varient=list()


for i in range(0, len(id2)):
    Tag2='-'
    kk=0
    vv=0
    temp1=id2[i]#ID for all detected

    for j in  range(0, len(P_list)):
        if seq2[i] == P_list[j]:#If sequence is match with P_criteria
           Tag2='P_criteria'
           #print(id2[i],seq2[i],P_list[j])
 #          f.writelines(id2[i]+'\n')
    #------------------------------------

    #MATCH with Known ID-----------------
    known_count=0
    for j in  range(0, len(known_id)):
        if id2[i] == known_id[j]:
           known_count= known_count+1
    #------------------------------------

    
    #MATCH with Varients-----------------
    varient_count=0
    for j in  range(0, len(varient_id)):
        if (id2[i] == varient_id[j]):
           varient_count=varient_count+1
    #------------------------------------


    #MATCH with Novel one ---------------
    novel_count=0
    for j in  range(0, len(Novel_id)):
        if (id2[i] == Novel_id[j]):
           novel_count=novel_count+1
    #------------------------------------



    #------------------------------------
    if known_count==1:Tag='Known'
    if varient_count==1:Tag='Varient'
    if novel_count==1:Tag='Novel'        
    #------------------------------------


    if Tag=='Known':
       temp_parse=''
       for k in range(0, len(parse_id)):       
            if id2[i] == parse_id[k]:
               temp_parse=parse_match[k]
               #print(i,Tag,Tag2,seq2[i],parse_id[k],parse_match[k])
               kk=kk+1

       if Tag2=='P_criteria':
          Tag_P_criteria=Tag_P_criteria+1
          DE_P_known_len.append(len(seq2[i]))
          #print(seq2[i])
          DE_P_known_seq.append(seq2[i])
          #print(i,Tag,Tag2,seq2[i],parse_id[i],parse_match[i])
          #print(parse_match[i])
          miR_known.append(parse_id[i]+':'+parse_match[i])


    if Tag=='Varient':
       temp_parse=''
       for k in range(0, len(parse_id)):
            if id2[i] == parse_id[k]:
               temp_parse=parse_match[k]
  #             print(i,Tag,Tag2,seq2[i],parse_id[k],parse_match[k])
               vv=vv+1
       if Tag2=='P_criteria':
          Tag_P_criteria2=Tag_P_criteria2+1
          DE_P_varient_len.append(len(seq2[i]))
          DE_P_varient_seq.append(seq2[i])
          #print(i,Tag,Tag2,seq2[i],parse_id[i],parse_match[i])
          miR_varient.append(parse_id[i]+':'+parse_match[i])
    if Tag=='Novel':
       if Tag2=='P_criteria':
          Novel_P_criteria = Novel_P_criteria +1
          DE_P_novel_len.append(len(seq2[i]))
          DE_P_novel_seq.append(seq2[i])          


#f.close()

#exit()
P_Dknown,P_Dknown2 = numpy.unique(DE_P_known_len,return_counts=True)
P_Dvarient,P_Dvarient2 = numpy.unique(DE_P_varient_len,return_counts=True)
P_Dnovel,P_Dnovel2 = numpy.unique(DE_P_novel_len,return_counts=True)
M=100/(sum(P_Dknown2) + sum(P_Dvarient2) + sum(P_Dnovel2) )

print('-------------------P+DE+known-----------------------------------------------')
#P_Dknown,P_Dknown2 = numpy.unique(DE_P_known_len,return_counts=True)
print(P_Dknown)
print(P_Dknown2*M)
print('-------------------P+DE+varient---------------------------------------------')
#P_Dvarient,P_Dvarient2 = numpy.unique(DE_P_varient_len,return_counts=True)
print(P_Dvarient)
print(P_Dvarient2*M)
print('-------------------P+DE+novel-----------------------------------------------')
#P_Dnovel,P_Dnovel2 = numpy.unique(DE_P_novel_len,return_counts=True)
print(P_Dnovel)
print(P_Dnovel2*M)
#print('--------------------------------------------------------------------------------------')

print('DE+P_criteria: ', len(P_uniq))
print('DE+known+P_criteria: ',Tag_P_criteria)
print('DE+Varient+P_criteria: ',Tag_P_criteria2)
print('DE+Novel_P_criteria: ',Novel_P_criteria)
#print('Novel_P_Uniq: ',len(P_uniq))
#print('--------------------------------------------------------------------------------------')



#rint('---------------Known +P+DE frequency------------------------- --------------')    
A18=0;T18=0;C18=0;G18=0;A19=0;T19=0;C19=0;G19=0;A20=0;T20=0;C20=0;G20=0
A21=0;T21=0;C21=0;G21=0;A22=0;T22=0;C22=0;G22=0;A23=0;T23=0;C23=0;G23=0
A24=0;T24=0;C24=0;G24=0
#DE_P_known_seq

#known_id_seq=DE_P_known_seq
for i in range(0, len(DE_P_known_seq)):

    if len(DE_P_known_seq[i])==18:
       if DE_P_known_seq[i][0] =='A':A18=A18+1
       if DE_P_known_seq[i][0] =='T':T18=T18+1
       if DE_P_known_seq[i][0] =='C':C18=C18+1
       if DE_P_known_seq[i][0] =='G':G18=G18+1
    if len(DE_P_known_seq[i])==19:
       if DE_P_known_seq[i][0]=='A':A19=A19+1
       if DE_P_known_seq[i][0]=='T':T19=T19+1
       if DE_P_known_seq[i][0]=='C':C19=C19+1
       if DE_P_known_seq[i][0]=='G':G19=G19+1    
    if len(DE_P_known_seq[i])==20:
       if DE_P_known_seq[i][0]=='A':A20=A20+1
       if DE_P_known_seq[i][0]=='T':T20=T20+1
       if DE_P_known_seq[i][0]=='C':C20=C20+1
       if DE_P_known_seq[i][0]=='G':G20=G20+1
    if len(DE_P_known_seq[i])==21:
       if DE_P_known_seq[i][0] =='A':A21=A21+1
       if DE_P_known_seq[i][0] =='T':T21=T21+1
       if DE_P_known_seq[i][0] =='C':C21=C21+1
       if DE_P_known_seq[i][0] =='G':G21=G21+1
    if len(DE_P_known_seq[i])==22:
       if DE_P_known_seq[i][0] =='A':A22=A22+1
       if DE_P_known_seq[i][0] =='T':T22=T22+1
       if DE_P_known_seq[i][0] =='C':C22=C22+1
       if DE_P_known_seq[i][0] =='G':G22=G22+1
    if len(DE_P_known_seq[i])==23:
       if DE_P_known_seq[i][0] =='A':A23=A23+1
       if DE_P_known_seq[i][0] =='T':T23=T23+1
       if DE_P_known_seq[i][0] =='C':C23=C23+1
       if DE_P_known_seq[i][0] =='G':G23=G23+1
    if len(DE_P_known_seq[i])==24:
       if DE_P_known_seq[i][0] =='A':A24=A24+1
       if DE_P_known_seq[i][0] =='T':T24=T24+1
       if DE_P_known_seq[i][0] =='C':C24=C24+1
       if DE_P_known_seq[i][0] =='G':G24=G24+1

#print('A18',A18,'C18',C18,'G18',G18,'T18',T18)
#print('A19',A19,'C19',C19,'G19',G19,'T19',T19)
#print('A20',A20,'C20',C20,'G20',G20,'T20',T20)
#print('A21',A21,'C21',C21,'G21',G21,'T21',T21)
#print('A22',A22,'C22',C22,'G22',G22,'T22',T22)
#print('A23',A23,'C23',C23,'G23',G23,'T23',T23)
#print('A24',A24,'C24',C24,'G24',G24,'T24',T24)


#print('-------Variants +P+DE frequency --------------------------------------------')
A18=0;T18=0;C18=0;G18=0;A19=0;T19=0;C19=0;G19=0;A20=0;T20=0;C20=0;G20=0
A21=0;T21=0;C21=0;G21=0;A22=0;T22=0;C22=0;G22=0;A23=0;T23=0;C23=0;G23=0
A24=0;T24=0;C24=0;G24=0

for i in range(0, len(DE_P_varient_seq)):
    #print(i,len(DE_P_varient_seq[i]))
    if len(DE_P_varient_seq[i])==18:
       if DE_P_varient_seq[i][0] =='A':A18=A18+1
       if DE_P_varient_seq[i][0] =='T':T18=T18+1
       if DE_P_varient_seq[i][0] =='C':C18=C18+1
       if DE_P_varient_seq[i][0] =='G':G18=G18+1
    if len(DE_P_varient_seq[i])==19:
       if DE_P_varient_seq[i][0] =='A':A19=A19+1
       if DE_P_varient_seq[i][0] =='T':T19=T19+1
       if DE_P_varient_seq[i][0] =='C':C19=C19+1
       if DE_P_varient_seq[i][0] =='G':G19=G19+1
    if len(DE_P_varient_seq[i])==20:
       if DE_P_varient_seq[i][0] =='A':A20=A20+1
       if DE_P_varient_seq[i][0] =='T':T20=T20+1
       if DE_P_varient_seq[i][0] =='C':C20=C20+1
       if DE_P_varient_seq[i][0] =='G':G20=G20+1
    if len(DE_P_varient_seq[i])==21:
       if DE_P_varient_seq[i][0] =='A':A21=A21+1
       if DE_P_varient_seq[i][0] =='T':T21=T21+1
       if DE_P_varient_seq[i][0] =='C':C21=C21+1
       if DE_P_varient_seq[i][0] =='G':G21=G21+1
    if len(DE_P_varient_seq[i])==22:
       if DE_P_varient_seq[i][0] =='A':A22=A22+1
       if DE_P_varient_seq[i][0] =='T':T22=T22+1
       if DE_P_varient_seq[i][0] =='C':C22=C22+1
       if DE_P_varient_seq[i][0] =='G':G22=G22+1
    if len(DE_P_varient_seq[i])==23:
       if DE_P_varient_seq[i][0] =='A':A23=A23+1
       if DE_P_varient_seq[i][0] =='T':T23=T23+1
       if DE_P_varient_seq[i][0] =='C':C23=C23+1
       if DE_P_varient_seq[i][0] =='G':G23=G23+1
    if len(DE_P_varient_seq[i])==24:
       if DE_P_varient_seq[i][0] =='A':A24=A24+1
       if DE_P_varient_seq[i][0] =='T':T24=T24+1
       if DE_P_varient_seq[i][0] =='C':C24=C24+1
       if DE_P_varient_seq[i][0] =='G':G24=G24+1

#print('A18',A18,'C18',C18,'G18',G18,'T18',T18)
#print('A19',A19,'C19',C19,'G19',G19,'T19',T19)
#print('A20',A20,'C20',C20,'G20',G20,'T20',T20)
#print('A21',A21,'C21',C21,'G21',G21,'T21',T21)
#print('A22',A22,'C22',C22,'G22',G22,'T22',T22)
#print('A23',A23,'C23',C23,'G23',G23,'T23',T23)
#print('A24',A24,'C24',C24,'G24',G24,'T24',T24)


#print('----------------------Novel+P+DE--------------------------------------------')
A18=0;T18=0;C18=0;G18=0;A19=0;T19=0;C19=0;G19=0;A20=0;T20=0;C20=0;G20=0
A21=0;T21=0;C21=0;G21=0;A22=0;T22=0;C22=0;G22=0;A23=0;T23=0;C23=0;G23=0
A24=0;T24=0;C24=0;G24=0

for i in range(0, len(DE_P_novel_seq)):
    if len(DE_P_novel_seq[i])==18:
       if DE_P_novel_seq[i][0] =='A':A18=A18+1
       if DE_P_novel_seq[i][0] =='T':T18=T18+1
       if DE_P_novel_seq[i][0] =='C':C18=C18+1
       if DE_P_novel_seq[i][0] =='G':G18=G18+1
    if len(DE_P_novel_seq[i])==19:
       if DE_P_novel_seq[i][0] =='A':A19=A19+1
       if DE_P_novel_seq[i][0] =='T':T19=T19+1
       if DE_P_novel_seq[i][0] =='C':C19=C19+1
       if DE_P_novel_seq[i][0] =='G':G19=G19+1
    if len(DE_P_novel_seq[i])==20:
       if DE_P_novel_seq[i][0] =='A':A20=A20+1
       if DE_P_novel_seq[i][0] =='T':T20=T20+1
       if DE_P_novel_seq[i][0] =='C':C20=C20+1
       if DE_P_novel_seq[i][0] =='G':G20=G20+1
    if len(DE_P_novel_seq[i])==21:
       if DE_P_novel_seq[i][0] =='A':A21=A21+1
       if DE_P_novel_seq[i][0] =='T':T21=T21+1
       if DE_P_novel_seq[i][0] =='C':C21=C21+1
       if DE_P_novel_seq[i][0] =='G':G21=G21+1
    if len(DE_P_novel_seq[i])==22:
       if DE_P_novel_seq[i][0] =='A':A22=A22+1
       if DE_P_novel_seq[i][0] =='T':T22=T22+1
       if DE_P_novel_seq[i][0] =='C':C22=C22+1
       if DE_P_novel_seq[i][0] =='G':G22=G22+1
    if len(DE_P_novel_seq[i])==23:
       if DE_P_novel_seq[i][0] =='A':A23=A23+1
       if DE_P_novel_seq[i][0] =='T':T23=T23+1
       if DE_P_novel_seq[i][0] =='C':C23=C23+1
       if DE_P_novel_seq[i][0] =='G':G23=G23+1
    if len(DE_P_novel_seq[i])==24:
       if DE_P_novel_seq[i][0] =='A':A24=A24+1
       if DE_P_novel_seq[i][0] =='T':T24=T24+1
       if DE_P_novel_seq[i][0] =='C':C24=C24+1
       if DE_P_novel_seq[i][0] =='G':G24=G24+1


#print('A18',A18,'C18',C18,'G18',G18,'T18',T18)
#print('A19',A19,'C19',C19,'G19',G19,'T19',T19)
#print('A20',A20,'C20',C20,'G20',G20,'T20',T20)
#print('A21',A21,'C21',C21,'G21',G21,'T21',T21)
#print('A22',A22,'C22',C22,'G22',G22,'T22',T22)
#print('A23',A23,'C23',C23,'G23',G23,'T23',T23)
#print('A24',A24,'C24',C24,'G24',G24,'T24',T24)
                 
print('---------------P_criteria+DE------------------------------------------------')
#P_uniq_len=list()
#P_uniq_fist_base=list()

A18=0;T18=0;C18=0;G18=0;A19=0;T19=0;C19=0;G19=0;A20=0;T20=0;C20=0;G20=0
A21=0;T21=0;C21=0;G21=0;A22=0;T22=0;C22=0;G22=0;A23=0;T23=0;C23=0;G23=0
A24=0;T24=0;C24=0;G24=0

for i in range(0, len(P_list)):
    #print(P_list[i],len(P_list[i]),P_list[i][0])
    base=P_list[i][0]
    if len(P_list[i])==18:
       if base=='A':A18=A18+1
       if base=='T':T18=T18+1
       if base=='C':C18=C18+1
       if base=='G':G18=G18+1
    if len(P_list[i])==19:
       if base=='A':A19=A19+1
       if base=='T':T19=T19+1
       if base=='C':C19=C19+1
       if base=='G':G19=G19+1
    if len(P_list[i])==20:
       if base=='A':A20=A20+1
       if base=='T':T20=T20+1
       if base=='C':C20=C20+1
       if base=='G':G20=G20+1
    if len(P_list[i])==21:
       if base=='A':A21=A21+1
       if base=='T':T21=T21+1
       if base=='C':C21=C21+1
       if base=='G':G21=G21+1
    if len(P_list[i])==22:
       if base=='A':A22=A22+1
       if base=='T':T22=T22+1
       if base=='C':C22=C22+1
       if base=='G':G22=G22+1
    if len(P_list[i])==23:
       if base=='A':A23=A23+1
       if base=='T':T23=T23+1
       if base=='C':C23=C23+1
       if base=='G':G23=G23+1
    if len(P_list[i])==24:
       if base=='A':A24=A24+1
       if base=='T':T24=T24+1
       if base=='C':C24=C24+1
       if base=='G':G24=G24+1

       
print('A18',A18,'C18',C18,'G18',G18,'T18',T18)
print('A19',A19,'C19',C19,'G19',G19,'T19',T19)
print('A20',A20,'C20',C20,'G20',G20,'T20',T20)
print('A21',A21,'C21',C21,'G21',G21,'T21',T21)
print('A22',A22,'C22',C22,'G22',G22,'T22',T22)
print('A23',A23,'C23',C23,'G23',G23,'T23',T23)
print('A24',A24,'C24',C24,'G24',G24,'T24',T24)



#print('-----------Known Parse id ------------------')

#print(len(miR_known))

#for i in range(0, len(miR_known)):
 #   temp=miR_known[i].split('&')
    #print(temp)
#    print(i,temp[0],len(temp)-1)

#print('-----------Varient Parse id ------------------')

#for i in range(0, len(miR_varient)):
 #   temp=miR_varient[i].split('&')
    #print(temp)
#    print(i,temp[0],len(temp)-1)
print('############################################################################')
