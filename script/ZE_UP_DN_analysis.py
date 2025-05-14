#!/opt/cray/pe/python/3.11.5/bin/python
#- Dr Ahsan Z Rizvi
#-------------------------
import numpy
import sys
import os
import numpy as np
#import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import mstats
#-------------------------


miRNA_name=list()
miRNA_exp=list()
with open(sys.argv[2]) as file:
    for line in file:
        #print(line)
        miRNA_exp.append(line)
        temp1=line.split('\t')
        #print(temp1[0])
        miRNA_name.append(temp1[0])

#sys.exit()
#known_seq_len_nt, known_seq_len_nt_freq = numpy.unique(miRNA_name,return_counts=True)
#-print(known_seq_len_nt)
#-print(known_seq_len_nt_freq)

#-sys.exit()

mRNA_name=list()
mRNA_exp=list()
with open(sys.argv[3]) as file:
    for line in file:
 #       -print(line)
        mRNA_exp.append(line)
        temp1=line.split('\t')
        temp6=(temp1[0].split('.'))
        #print(temp1[0])
        mRNA_name.append(temp1[0])
#sys.exit()

#-known_seq_len_nt, known_seq_len_nt_freq = numpy.unique(mRNA_name,return_counts=True)
#-print(known_seq_len_nt)
#-print(known_seq_len_nt_freq)

#sys.exit()

miRNA_2=list()
miRNA_3=list()
miRNA_4=list()
miRNA_5=list()
miRNA_6=list()
miRNA_7=list()
miRNA_8=list()
miRNA_9=list()
miRNA_10=list()
miRNA_11=list()
miRNA_12=list()
miRNA_13=list()
miRNA_14=list()
miRNA_15=list()
miRNA_16=list()
miRNA_17=list()



mRNA_2=list()
mRNA_3=list()
mRNA_4=list()
mRNA_5=list()
mRNA_6=list()
mRNA_7=list()
mRNA_8=list()
mRNA_9=list()
mRNA_10=list()
mRNA_11=list()
mRNA_12=list()
mRNA_13=list()
mRNA_14=list()
mRNA_15=list()
mRNA_16=list()
mRNA_17=list()

i=0
miRNA_pair=list()
mRNA_pair=list()
k=0
with open(sys.argv[1]) as file:
    for line in file:
        #if k>500:
         #   break
        #print(k)
        #k=k+1
       # print(line)
        temp1=(line.split('\t'))
        if i>1:
           if len(temp1)>1 and float(temp1[2]) >4:
              
              for j in range(0, len(miRNA_name)):
                  if temp1[0] == miRNA_name[j]:
                     #print(i,temp1[0],temp1[1],miRNA_exp[j])
                     temp2=temp1[0]
                     temp3=miRNA_exp[j]
              for k in range(0, len(mRNA_name)):
                  if temp1[1] == mRNA_name[k]:
                     #print(i,temp1[2],temp1[1],mRNA_exp[k])
                     temp4=temp1[1]
                     temp5=mRNA_exp[k]
              temp31=temp3.split('\t')
              #print((temp31))

              miRNA_2.append(temp31[1])
              #print(temp31[1])
              miRNA_3.append(temp31[2])
              miRNA_4.append(temp31[3])
              miRNA_5.append(temp31[4])
              miRNA_6.append(temp31[5])
              miRNA_7.append(temp31[6])
              miRNA_8.append(temp31[7]) 
              miRNA_9.append(temp31[8])
              miRNA_10.append(temp31[9])
              miRNA_11.append(temp31[10])
              miRNA_12.append(temp31[11])
              miRNA_13.append(temp31[12])
              miRNA_14.append(temp31[13])
              miRNA_15.append(temp31[14])
              miRNA_16.append(temp31[15])
              temp100=(temp31[16].split('\n'))
              miRNA_17.append(temp100[0])
              


              temp51=temp5.split('\t')
              #print(len(temp51))

              mRNA_2.append(temp51[1])
              mRNA_3.append(temp51[2])
              mRNA_4.append(temp51[3])
              mRNA_5.append(temp51[4])
              mRNA_6.append(temp51[5])
              mRNA_7.append(temp51[6])
              mRNA_8.append(temp51[7])
              mRNA_9.append(temp51[8])
              mRNA_10.append(temp51[9])
              mRNA_11.append(temp51[10])
              mRNA_12.append(temp51[11])
              mRNA_13.append(temp51[12])
              mRNA_14.append(temp51[13])
              mRNA_15.append(temp51[14])
              mRNA_16.append(temp51[15])
              temp200=(temp51[16].split('\n'))
              #print(temp200[0])
              mRNA_17.append(temp200[0])
        i=i+1
        #if i>100:
         #  break

#sys.exit()

#--------------------------------------------------
#--ANALYSIS
#--------------------------------------------------


#print('Total pairs: ',i,'Passed CC negative: ',mm,'Passed P test: ',nn)


tempS211=0;tempS311=0;tempS411=0;tempS511=0;tempS611=0;tempS711=0;tempS811=0
tempS911=0;tempS1011=0;tempS1111=0;tempS1211=0;tempS1311=0;tempS1411=0;tempS1511=0
tempS1611=0;tempS1711=0
#Pattern='DN-UP'
Pattern='UP-DN'
print('Pattern wish to check: ',Pattern)
print('Total pairs: ',len(miRNA_2) )
for i in range(0, len(miRNA_2)):
    if float(miRNA_2[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_2[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    tempS2=(temp11+'-'+temp21)
    if tempS2==Pattern:
       tempS21=1
    else:
       tempS21=0
    tempS211=tempS211+tempS21
    if float(miRNA_3[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_3[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    tempS3=(temp11+'-'+temp21)
    if tempS3==Pattern:
       tempS31=1
    else:
       tempS31=0
    tempS311=tempS311+tempS31
    if float(miRNA_4[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_4[i])>0:
       temp21='UP'
    else:
      temp21='DN'
    tempS4=(temp11+'-'+temp21) ##

    if tempS4==Pattern:
       tempS41=1
    else:
       tempS41=0
    tempS411=tempS411+tempS41
    if float(miRNA_5[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_5[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    tempS5=(temp11+'-'+temp21)
    if tempS5==Pattern:
       tempS51=1
    else:
       tempS51=0
    tempS511=tempS511+tempS51
    if float(miRNA_6[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_6[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    tempS6=(temp11+'-'+temp21)

    if tempS6==Pattern:
       tempS61=1
    else:
       tempS61=0
    tempS611=tempS611+tempS61
    if float(miRNA_7[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_7[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    tempS7=(temp11+'-'+temp21) ##

    if tempS7==Pattern:
       tempS71=1
    else:
       tempS71=0
    tempS711=tempS711+tempS71

    #------------------------
    if float(miRNA_8[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_8[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    ##
    tempS8=(temp11+'-'+temp21)
    if tempS8==Pattern:
       tempS81=1
    else:
       tempS81=0
    tempS811=tempS811+tempS81
    #------------------------
    
    
    #------------------------
    if float(miRNA_9[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_9[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    ##
    tempS9=(temp11+'-'+temp21)
    if tempS9==Pattern:
       tempS91=1
    else:
       tempS91=0
    tempS911=tempS911+tempS91
    #------------------------
    #------------------------
    if float(miRNA_10[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_10[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    ##
    tempS10=(temp11+'-'+temp21)
    if tempS10==Pattern:
       tempS101=1
    else:
       tempS101=0
    tempS1011=tempS1011+tempS101
    #------------------------
        #------------------------
    if float(miRNA_11[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_11[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    ##
    tempS11=(temp11+'-'+temp21)
    if tempS11==Pattern:
       tempS111=1
    else:
       tempS111=0
    tempS1111=tempS1111+tempS111
    #------------------------

        #------------------------
    if float(miRNA_12[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_12[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    ##
    tempS12=(temp11+'-'+temp21)
    if tempS12==Pattern:
       tempS121=1
    else:
       tempS121=0
    tempS1211=tempS1211+tempS121
    #------------------------

        #------------------------
    if float(miRNA_13[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_13[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    ##
    tempS13=(temp11+'-'+temp21)
    if tempS13==Pattern:
       tempS131=1
    else:
       tempS131=0
    tempS1311=tempS1311+tempS131
    #------------------------
        #------------------------
    if float(miRNA_14[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_14[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    ##
    tempS14=(temp11+'-'+temp21)
    if tempS14==Pattern:
       tempS141=1
    else:
       tempS141=0
    tempS1411=tempS1411+tempS141
    #------------------------
        #------------------------
    if float(miRNA_15[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_15[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    ##
    tempS15=(temp11+'-'+temp21)
    if tempS15==Pattern:
       tempS151=1
    else:
       tempS151=0
    tempS1511=tempS1511+tempS151
    #------------------------

        #------------------------
    if float(miRNA_16[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_16[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    ##
    tempS16=(temp11+'-'+temp21)
    if tempS16==Pattern:
       tempS161=1
    else:
       tempS161=0
    tempS1611=tempS1611+tempS161
    #------------------------

        #------------------------
    if float(miRNA_17[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_17[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    ##
    tempS17=(temp11+'-'+temp21)
    if tempS17==Pattern:
       tempS171=1
    else:
       tempS171=0
    tempS1711=tempS1711+tempS171
    #------------------------


#sys.exit()

Result=sys.argv[4]+'/ZE_UP_DN.result.txt'
f = open(Result, "w")


f.write(Pattern+'------------------------------------'+'\n')
print(Pattern,'-------------------------------------------------------------------------------------')

f.write(str(tempS211)+'\t'+str(tempS311)+'\t'+str(tempS411)+'\t'+str(tempS511)+'\t'+str(tempS611)+'\t'+str(tempS711)+'\t'+str(tempS811)+'\t'+str(tempS911)+'\t'+str(tempS1011)+'\t'+str(tempS1111)+'\t'+str(tempS1211)+'\t'+str(tempS1311)+'\t'+str(tempS1411)+'\t'+str(tempS1511)+'\t'+str(tempS1611)+'\t'+str(tempS1711)   +'\n')


print(tempS211,tempS311,tempS411,tempS511,tempS611,tempS711,tempS811,tempS911,tempS1011,tempS1111,tempS1211,tempS1311,tempS1411,tempS1511,tempS1611,tempS1711)






tempS211=0;tempS311=0;tempS411=0;tempS511=0;tempS611=0;tempS711=0;tempS811=0
tempS911=0;tempS1011=0;tempS1111=0;tempS1211=0;tempS1311=0;tempS1411=0;tempS1511=0
tempS1611=0;tempS1711=0


Pattern='DN-UP'
#Pattern='UP-DN'
print('Pattern wish to check: ',Pattern)
print('Total pairs: ',len(miRNA_2) )
for i in range(0, len(miRNA_2)):
    if float(miRNA_2[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_2[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    tempS2=(temp11+'-'+temp21)
    if tempS2==Pattern:
       tempS21=1
    else:
       tempS21=0
    tempS211=tempS211+tempS21
    if float(miRNA_3[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_3[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    tempS3=(temp11+'-'+temp21)
    if tempS3==Pattern:
       tempS31=1
    else:
       tempS31=0
    tempS311=tempS311+tempS31
    if float(miRNA_4[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_4[i])>0:
       temp21='UP'
    else:
      temp21='DN'
    tempS4=(temp11+'-'+temp21) ##

    if tempS4==Pattern:
       tempS41=1
    else:
       tempS41=0
    tempS411=tempS411+tempS41
    if float(miRNA_5[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_5[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    tempS5=(temp11+'-'+temp21)
    if tempS5==Pattern:
       tempS51=1
    else:
       tempS51=0
    tempS511=tempS511+tempS51
    if float(miRNA_6[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_6[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    tempS6=(temp11+'-'+temp21)

    if tempS6==Pattern:
       tempS61=1
    else:
       tempS61=0
    tempS611=tempS611+tempS61
    if float(miRNA_7[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_7[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    tempS7=(temp11+'-'+temp21) ##

    if tempS7==Pattern:
       tempS71=1
    else:
       tempS71=0
    tempS711=tempS711+tempS71
    if float(miRNA_8[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_8[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    tempS8=(temp11+'-'+temp21)

    if tempS8==Pattern:
       tempS81=1
    else:
       tempS81=0
    tempS811=tempS811+tempS81

        #------------------------
    if float(miRNA_9[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_9[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    ##
    tempS9=(temp11+'-'+temp21)
    if tempS9==Pattern:
       tempS91=1
    else:
       tempS91=0
    tempS911=tempS911+tempS91
    #------------------------
    #------------------------
    if float(miRNA_10[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_10[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    ##
    tempS10=(temp11+'-'+temp21)
    if tempS10==Pattern:
       tempS101=1
    else:
       tempS101=0
    tempS1011=tempS1011+tempS101
    #------------------------
        #------------------------
    if float(miRNA_11[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_11[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    ##
    tempS11=(temp11+'-'+temp21)
    if tempS11==Pattern:
       tempS111=1
    else:
       tempS111=0
    tempS1111=tempS1111+tempS111
    #------------------------

        #------------------------
    if float(miRNA_12[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_12[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    ##
    tempS12=(temp11+'-'+temp21)
    if tempS12==Pattern:
       tempS121=1
    else:
       tempS121=0
    tempS1211=tempS1211+tempS121
    #------------------------
        #------------------------
    if float(miRNA_13[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_13[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    ##
    tempS13=(temp11+'-'+temp21)
    if tempS13==Pattern:
       tempS131=1
    else:
       tempS131=0
    tempS1311=tempS1311+tempS131
    #------------------------
        #------------------------
    if float(miRNA_14[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_14[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    ##
    tempS14=(temp11+'-'+temp21)
    if tempS14==Pattern:
       tempS141=1
    else:
       tempS141=0
    tempS1411=tempS1411+tempS141
    #------------------------
        #------------------------
    if float(miRNA_15[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_15[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    ##
    tempS15=(temp11+'-'+temp21)
    if tempS15==Pattern:
       tempS151=1
    else:
       tempS151=0
    tempS1511=tempS1511+tempS151
    #------------------------

    if float(miRNA_16[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_16[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    ##
    tempS16=(temp11+'-'+temp21)
    if tempS16==Pattern:
       tempS161=1
    else:
       tempS161=0
    tempS1611=tempS1611+tempS161
    #------------------------

        #------------------------
    if float(miRNA_17[i])>0:
       temp11='UP'
    else:
       temp11='DN'
    if float(mRNA_17[i])>0:
       temp21='UP'
    else:
       temp21='DN'
    ##
    tempS17=(temp11+'-'+temp21)
    if tempS17==Pattern:
       tempS171=1
    else:
       tempS171=0
    tempS1711=tempS1711+tempS171
    #------------------------






print(Pattern,'-------------------------------------------------------------------------------------')
#print(tempS211,tempS311,tempS411,tempS511,tempS611,tempS711,tempS811)



f.write(Pattern+'------------------------------------'+'\n')

#f.write(str(tempS211)+'\t'+str(tempS311)+'\t'+str(tempS411)+'\t'+str(tempS511)+'\t'+str(tempS611)+'\t'+str(tempS711)+'\t'+str(tempS811))


#f.write(str(tempS211)+'\t'+str(tempS311)+'\t'+str(tempS411)+'\t'+str(tempS511)+'\t'+str(tempS611)+'\t'+str(tempS711)+'\t'+str(tempS811)+'\t'+str(tempS911)+'\t'+str(tempS1011)+'\t'+str(tempS1111)+'\t'+str(tempS1211)+'\t'+str(tempS1311)+'\t'+str(tempS1411)+'\t'+str(tempS1511)+'\n')


#print(tempS211,tempS311,tempS411,tempS511,tempS611,tempS711,tempS811,tempS911,tempS1011,tempS1111,tempS1211,tempS1311,tempS1411,tempS1511)


f.write(str(tempS211)+'\t'+str(tempS311)+'\t'+str(tempS411)+'\t'+str(tempS511)+'\t'+str(tempS611)+'\t'+str(tempS711)+'\t'+str(tempS811)+'\t'+str(tempS911)+'\t'+str(tempS1011)+'\t'+str(tempS1111)+'\t'+str(tempS1211)+'\t'+str(tempS1311)+'\t'+str(tempS1411)+'\t'+str(tempS1511)+'\t'+str(tempS1611)+'\t'+str(tempS1711)   +'\n')


print(tempS211,tempS311,tempS411,tempS511,tempS611,tempS711,tempS811,tempS911,tempS1011,tempS1111,tempS1211,tempS1311,tempS1411,tempS1511,tempS1611,tempS1711)











f.close()














