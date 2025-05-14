#!/opt/cray/pe/python/3.11.5/bin/python
#- Dr Ahsan Z Rizvi
#-------------------------
import numpy
import sys
#import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
#from scipy import stats
#from scipy.stats import mstats
import pandas as pd
import scipy.cluster.hierarchy as sch
from inspect import currentframe, getframeinfo
#-------------------------
with open(sys.argv[1]) as file:#GroupA reading
    for line in file:
        GroupA=(line.split(','))
print('Total miRNA member in Group A:',len(GroupA))
with open(sys.argv[2]) as file:#GroupB reading
    for line in file:
        GroupB=(line.split(','))
GroupB=GroupB[:-1]

print('Total miRNA member in Group B:',len(GroupB))
with open(sys.argv[3]) as file:#GroupC reading
    for line in file:
        GroupC=(line.split(','))
print('Total miRNA member in Group C:',len(GroupC))
Known=list()
with open(sys.argv[4]) as file:#Known reading
    for line in file:
        Temp1=(line.split('\n'))
       # print(Temp1[0])
        Known.append(Temp1[0])
        #Known=(line.split(','))
#print('Total  miRNA in Known:',len(Known))
resA = list(set(GroupA) & set(Known))
#print('Known miRNA in GroupA',len(resA))
resB = list(set(GroupB) & set(Known))
#print('Known miRNA in GroupB',len(resB))
resC = list(set(GroupC) & set(Known))
#print('Known miRNA in GroupC',len(resC))
Group_all=GroupA+GroupB+GroupC # Concatenate all groups
#Group_all_Unique=np.unique(Group_all)
#print('Length in Group',len(Group_all))
#print('Length unique in GroupC',len(Group_all_Unique))
#################################

#sys.exit()

############################# Read mirBASE IDs
i=0
miRNA1=list()
mirBASE1=list()
with open(sys.argv[5]) as file:
    for line in file:
        temp=line.split('\n')
        temp2=temp[0].split(',')
        miRNA1.append(temp2[0])
        #print(temp2[0])
        mirBASE1.append(temp)
        #print(i,temp)
        i=i+1
miRNA1_len=len(miRNA1)
#############################


############################ Read -ve CC pair file

miRNA_CC=list()
mRNA_CC=list()
CC_value=list()
P_value=list()
with open(sys.argv[7]) as file:
    for line in file:
        temp=line.split('\t')
        #print(temp[0])
        #print(temp[0],temp[1],temp[2],temp[3],temp[4])
        miRNA_CC.append(temp[0])
        mRNA_CC.append(temp[1])
        CC_value.append(temp[3])
        P_value.append(temp[4])
miRNA_CC_nonRed,miRNA_CC_freq = numpy.unique(miRNA_CC,return_counts=True)
miRNA11=list()
mirBASE11=list()
###########################
#sys.exit()

########################### Extracting mRNA of each group
miRNA_CC_GpA=list()#
mRNA_CC_GpA=list()#

temp2=sys.argv[7]+'.GpA.csv'
f = open(temp2, "w")
for i in range(0,len(GroupA)):#Extract GroupA:miRNA,mRNA,GO_ID
    for j in range(0,len(miRNA_CC)):
        if GroupA[i]==miRNA_CC[j]:
           temp=mRNA_CC[j].split('.')
           #print(i,miRNA_CC[j],temp[0],CC_value[j],P_value[j])
           temp3=miRNA_CC[j]+','+temp[0]+','+CC_value[j],','+P_value[j]#+'\n'
           miRNA_CC_GpA.append(miRNA_CC[j])#
           mRNA_CC_GpA.append(temp[0])#
           f.writelines(temp3)
f.close()

temp2=sys.argv[7]+'.GpA.mRNA-miRNA.tsv'
f = open(temp2, "w")
mRNA_GpA_nonRed,mRNA_GpA_freq = numpy.unique(mRNA_CC_GpA,return_counts=True)#
for i in range(0,len(mRNA_GpA_nonRed)):#
    temp=''
    for j in range(0,len(mRNA_CC_GpA)):
        if mRNA_GpA_nonRed[i]==mRNA_CC_GpA[j]:
           temp=temp+miRNA_CC_GpA[j]+';'
    #print(mRNA_GpA_nonRed[i],temp)
    temp3=mRNA_GpA_nonRed[i]+'\t'+temp+'\n'
    f.writelines(temp3)
f.close()
#######################################

####################################### GroupB mRNA
miRNA_CC_GpB=list()#
mRNA_CC_GpB=list()#
temp2=sys.argv[7]+'.GpB.csv'
f = open(temp2, "w")
for i in range(0,len(GroupB)):#Extract GroupA:miRNA,mRNA,GO_ID
    for j in range(0,len(miRNA_CC)):
        if GroupB[i]==miRNA_CC[j]:
           temp=mRNA_CC[j].split('.')
         #  print(miRNA_CC[j],temp[0],CC_value[j],P_value[j])
           temp3=miRNA_CC[j]+','+temp[0]+','+CC_value[j],','+P_value[j]#+'\n'
           miRNA_CC_GpB.append(miRNA_CC[j])#
           mRNA_CC_GpB.append(temp[0])#
           f.writelines(temp3)
f.close()
temp2=sys.argv[7]+'.GpB.mRNA-miRNA.tsv'
f = open(temp2, "w")
mRNA_GpB_nonRed,mRNA_GpB_freq = numpy.unique(mRNA_CC_GpB,return_counts=True)#
for i in range(0,len(mRNA_GpB_nonRed)):#
    temp=''
    for j in range(0,len(mRNA_CC_GpB)):
        if mRNA_GpB_nonRed[i]==mRNA_CC_GpB[j]:
           temp=temp+miRNA_CC_GpB[j]+';'
    #print(mRNA_GpA_nonRed[i],temp)
    temp3=mRNA_GpB_nonRed[i]+'\t'+temp+'\n'
    f.writelines(temp3)
f.close()
##################################





#temp2=sys.argv[7]+'.GpC.csv'
#f = open(temp2, "w")
for i in range(0,len(GroupC)):#Extract GroupA:miRNA,mRNA,GO_ID
    #print(GroupC[i])
    for j in range(0,len(miRNA_CC)):
       # print(miRNA_CC[j])
        if GroupC[i]==miRNA_CC[j]:
           print(GroupC[i])
    #       temp=mRNA_CC[j].split('.')
     #      #print(miRNA_CC[j],temp[0],CC_value[j],P_value[j])
      #     temp3=miRNA_CC[j]+','+temp[0]+','+CC_value[j],','+P_value[j]#+'\n'
       #    f.writelines(temp3)
#f.close()


print('miRNA-mRNA,E_value,NCC,P_value file(Total): ',sys.argv[7])
print('(GroupA) miRNA-mRNA,E_value,NCC,P_value file: ',sys.argv[7]+'.GpA.csv')
print('(GroupA) mRNA-miRNA file: ',sys.argv[7]+'.GpA.mRNA-miRNA.csv')
print('(GroupB) miRNA-mRNA,E_value,NCC,P_value file: ',sys.argv[7]+'.GpB.csv')
print('(GroupA) mRNA-miRNA file: ',sys.argv[7]+'.GpB.mRNA-miRNA.csv')
###########################


sys.exit()




########################### NAT_miRNA
Nat_miRNA=list()
with open(sys.argv[8]) as file:
     for line in file:
         temp=(line.split('\n'))
         #print(temp[0])
         Nat_miRNA.append(temp[0])
###########################
#print('NAT given miRNA of mRNA: ',len(Nat_miRNA))

########################### Detect mRNA target of NAT's list
temp2=sys.argv[8]+'.NegCCPair.csv'
f = open(temp2, "w")
for i in range(0,len(Nat_miRNA)):
    #print(Nat_miRNA[i])
    for j in range(0,len(miRNA_CC)):
        if Nat_miRNA[i]==miRNA_CC[j]:
           temp=mRNA_CC[j].split('.')
           #print(miRNA_CC[j],temp[0])
           temp3=miRNA_CC[j]+','+temp[0]+'\n'
           f.writelines(temp3)
f.close()
#print('NAT miRNA-mRNA -veCC pair generated file is: ',temp2)
###########################


########################### Extract Groups GO
GO_ID=list()
GO_mRNA=list()
with open(sys.argv[9]) as file:
     for line in file:
         temp=(line.split('\n'))
         temp2=temp[0].split(',')
         #print(temp2[0],temp2[1])
         GO_ID.append(temp2[0])
         GO_mRNA.append(temp2[1])
#print('No of mRNA-GO pairs given by Elena: ', len(GO_ID))
#print('No of unique GO Ids in Elena list: ',len(np.unique(GO_ID)))
#print('No of unique mRNAs in Elena list: ',len(np.unique(GO_mRNA)))

#sys.exit()
###########################
print('-----------------------------------------------------------------------------')
frameinfo = getframeinfo(currentframe())
print('(A) script name, line no to open for groupA,B and C (miRNA_ID,mRNA_ID,GO_ID):: ',frameinfo.filename,',', frameinfo.lineno+3)


GP='A'
#GP='B'
#GP='C'
#GP='X'
print('-----------------------------------------------------------------------------')


########################### Extract GroupA mRNA and GO
if GP=='A':
   temp3=sys.argv[9]+'.GO.GroupA.csv'
   f = open(temp3, "w")
   GroupA_mRNA=list()
   GroupA_mRNA_GO=list()
   for i in range(0,len(GroupA)):#Extract GroupA:miRNA,mRNA,GO_ID
       #if i==3:
        #   break
       for j in range(0,len(miRNA_CC)):
           if GroupA[i]==miRNA_CC[j]:
              temp=mRNA_CC[j].split('.')
              GroupA_mRNA.append(temp[0])
              for k in range(0,len(GO_mRNA)):
                  if GO_mRNA[k]==temp[0]:
                     temp2=GroupA[i]+','+GO_mRNA[k]+','+GO_ID[k]+'\n'
                     #print(GO_mRNA[k])
                     GroupA_mRNA_GO.append(GO_ID[k])
                     f.writelines(temp2)
   f.close()
   GroupA_mRNA_uniq=np.unique(GroupA_mRNA)#Unique list
   temp3=sys.argv[9]+'.mRNA.GroupA.csv'
   f = open(temp3, "w")
   for i in range(0,len(GroupA_mRNA_uniq)):
       temp2=GroupA_mRNA_uniq[i]+'\n'
       #print(temp2)
       f.writelines(temp2)
   f.close()
   print('Total number of mRNA id GroupA: ',len(GroupA_mRNA_uniq))
   print('Total number of mRNA-GO pairs in GroupA: ',len(GroupA_mRNA_GO))
   print('Total number of unique GO in pairs of GroupA: ',len(np.unique(GroupA_mRNA_GO)))

###########################


#sys.exit()



########################### GO list making for GroupB 
if GP=='B':
    temp3=sys.argv[9]+'.GO.GroupB.csv'
    f = open(temp3, "w")
    GroupB_mRNA=list()
    GroupB_mRNA_GO=list()
    for i in range(0,len(GroupB)):
        for j in range(0,len(miRNA_CC)):
            if GroupB[i]==miRNA_CC[j]:
               temp=mRNA_CC[j].split('.')
               GroupB_mRNA.append(temp[0])
               #print(miRNA_CC[j],temp[0])
               for k in range(0,len(GO_mRNA)):
                   if GO_mRNA[k]==temp[0]:
                      temp2=GroupB[i]+','+GO_mRNA[k]+','+GO_ID[k]+'\n'
                      #print(temp2)
                      GroupB_mRNA_GO.append(GO_ID[k])
                      f.writelines(temp2)
    f.close()
    GroupB_mRNA_uniq=np.unique(GroupB_mRNA)#Unique list
    temp3=sys.argv[9]+'.mRNA.GroupB.csv'
    f = open(temp3, "w")
    for i in range(0,len(GroupB_mRNA_uniq)):
        temp2=GroupB_mRNA_uniq[i]+'\n'
        #print(temp2)
        f.writelines(temp2)
    f.close()
    print('Total number of mRNA id GroupB: ',len(GroupB_mRNA_uniq))
    print('Total number of mRNA-GO pairs in GroupB: ',len(GroupB_mRNA_GO))
    print('Total number of unique GO in pairs of GroupB: ',len(np.unique(GroupB_mRNA_GO)))

###########################

########################### GO list making for GroupC
if GP=='C':
    miRNA_CC=list()
    mRNA_CC=list()
    with open(sys.argv[10]) as file:
        for line in file:
            temp=line.split('\t')
            miRNA_CC.append(temp[0])
            mRNA_CC.append(temp[1])
    miRNA_CC_nonRed,miRNA_CC_freq = numpy.unique(miRNA_CC,return_counts=True)
    miRNA11=list()
    mirBASE11=list()

    temp3=sys.argv[9]+'.GO.GroupC.csv'
    f = open(temp3, "w")
    GroupC_mRNA=list()
    GroupC_mRNA_GO=list()
    for i in range(0,len(GroupC)):
       # if i==5:
        #   break
        for j in range(0,len(miRNA_CC)):
            if GroupC[i]==miRNA_CC[j]:
               temp=mRNA_CC[j].split('.')
               GroupC_mRNA.append(temp[0])
             #  print(miRNA_CC[j],temp[0])
               for k in range(0,len(GO_mRNA)):
                   if GO_mRNA[k]==temp[0]:
                      temp2=GroupC[i]+','+GO_mRNA[k]+','+GO_ID[k]+'\n'
                      #print(temp2)
                      GroupC_mRNA_GO.append(GO_ID[k])
                      f.writelines(temp2)
    f.close()
    GroupC_mRNA_uniq=np.unique(GroupC_mRNA)#Unique list
    temp3=sys.argv[9]+'.mRNA.GroupC.csv'
    f = open(temp3, "w")
    for i in range(0,len(GroupC_mRNA_uniq)):
        temp2=GroupC_mRNA_uniq[i]+'\n'
        #print(temp2)
        f.writelines(temp2)
    f.close()
    print('Total number of mRNA id GroupC: ',len(GroupC_mRNA_uniq))
    print('Total number of mRNA-GO pairs in GroupC: ',len(GroupC_mRNA_GO))
    print('Total number of unique GO in pairs of GroupC: ',len(np.unique(GroupC_mRNA_GO)))


###########################
#sys.exit()


########################### Replace CC miRNA with group miRNA 
print('Open here heatmap for different groups')
#miRNA_CC_nonRed=GroupC
miRNA_CC_nonRed=GroupA
SaveImage='/cfs/klemming/home/a/ahsan2/script_pdc/SE_GroupA.pdf'

#miRNA_CC_nonRed=Nat_miRNA
#SaveImage='/cfs/klemming/home/a/ahsan2/script_pdc/TakeOut/NAT_SE_Group_CC_miRNA.pdf'
##########################



########################## Assign miBASE id to -ve CC miRNA 
for i in range(0,len(miRNA_CC_nonRed)): #make list on nonredundant miRNA CC
    #print(miRNA_CC_nonRed[i],miRNA_CC_freq[i])
    temp=miRNA_CC_nonRed[i]#+','+str(miRNA_CC_freq[i])
    for j in range(0,len(miRNA1)):
        if miRNA_CC_nonRed[i]==miRNA1[j]:
           temp=str(mirBASE1[j][0])#+','+str(miRNA_CC_freq[i])
    #print(miRNA_CC_nonRed[i],temp)
    #print(temp)
    miRNA11.append(miRNA_CC_nonRed[i])
    mirBASE11.append(temp)
############################







#sys.exit()




############################### Read VST of SE
line1=list()
miRNA_name_in_VST=list()
with open(sys.argv[6]) as file:
    for line in file:
        #print(line)
        line1.append(line)
        temp=line.split('\n')
        #print(temp)
        temp2=temp[0].split('\t')
        #print(temp2[0])
        miRNA_name_in_VST.append(temp2[0])
miRNA_name_in_VST_len=len(miRNA_name_in_VST)
###############################


legend_labels = ['GroupA', 'GroupB', 'GroupC', 'Known','Novel']
legend_colors = ['brown', 'gold', 'green', 'black','purple']

legend_handles = [plt.Rectangle((0, 0), 0, 0, color=color, label=label)
                  for color, label in zip(legend_colors, legend_labels)]

############################### Matrix for plots
matrix3=np.zeros((len(miRNA11),8))# Initiating an matrix for heatmap
y_axis=list()
data=list()
data2=list()
for i in range(0, len(miRNA11)):


    ############################ color setting Known Vs Novel
    temp2='green'
    for m in range(0,len(GroupA)):
        if GroupA[m]==miRNA11[i]:
           temp2='brown'
    for m in range(0,len(GroupB)):
        if GroupB[m]==miRNA11[i]:
           temp2='gold'
    data.append(temp2)
    ############################

    ########################### colour setting group
    temp='purple'
    for k in range(0,len(Known)):
        if Known[k]==miRNA11[i]:       
           temp='black'
    data2.append(temp)
    ############################
    
    print_index=mirBASE11[i]
    for j in range(0, len(miRNA_name_in_VST)):
        if  miRNA_name_in_VST[j]==miRNA11[i]:
            #print(miRNA11[i],print_index,miRNA_name_in_VST[j],line1[j])
            temp1=line1[j].split('\n')
            temp2=temp1[0].split('\t')
            #print(temp2)
            S1=numpy.median( [float(temp2[1]),float(temp2[2]),float(temp2[3]) ] )
            S2=numpy.median( [float(temp2[4]),float(temp2[5]),float(temp2[6]) ] )
            S3=numpy.median( [float(temp2[7]),float(temp2[8]),float(temp2[9]) ] )
            S4=numpy.median([float(temp2[10]),float(temp2[11]),float(temp2[12])])
            S5=numpy.median([float(temp2[13]),float(temp2[14]),float(temp2[15])])
            S6=numpy.median( [float(temp2[16]),float(temp2[17]),float(temp2[18])])
            S7=numpy.median([float(temp2[19]),float(temp2[20]) ] )
            S8=numpy.median([float(temp2[21]),float(temp2[22]),float(temp2[23])])
            Mean=np.mean([S1,S2,S3,S4,S5,S6,S7,S8])
            STD=np.std([S1,S2,S3,S4,S5,S6,S7,S8])
            #print(i,Mean,STD,print_index,temp2[0],S1,S2,S3,S4,S5,S6,S7,S8)
            if STD!=0:
               matrix3[i][0]=(S1-Mean)/STD
               matrix3[i][1]=(S2-Mean)/STD
               matrix3[i][2]=(S3-Mean)/STD
               matrix3[i][3]=(S4-Mean)/STD
               matrix3[i][4]=(S5-Mean)/STD
               matrix3[i][5]=(S6-Mean)/STD
               matrix3[i][6]=(S7-Mean)/STD
               matrix3[i][7]=(S8-Mean)/STD
            ##################################
            #Section to add min and max stage
            a=matrix3[i][:]   
            b=(np.where(a == a.max())   )
            c=(b[0]+1)
            b2=(np.where(a == a.min())   )
            c2=(b2[0]+1)
            #print(a)
            min_max=('S'+str(c)+'/'+'S'+str(c2))
            #print(min_max)
            #print_index=print_index+','+min_max# OFF it if u dont require min max ratio
            #print(print_index)
            ##################################
    y_axis.append(print_index)
####################################
#sys.exit()








####################################
ser = pd.Series(data=data, index=y_axis)
ser2 = pd.Series(data=data2, index=y_axis)
GroupAxis=[ser,ser2]



x_axis=['S1','S2','S3','S4','S5','S6','S7','S8']
plt.figure()
df = pd.DataFrame(data = matrix3,index = y_axis,columns = x_axis)
#CLUSTROGRAM
colors = ['#0000FF','#FFFFFF' ,"#ff0000"]
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)

cm=sns.clustermap(df,col_cluster=False, row_cluster=True,cmap=custom_cmap,method='ward', metric='euclidean', figsize=(20,20),row_colors=GroupAxis,dendrogram_ratio=(0.25, 0.15),cbar_pos=[.4, .95, .5, .03],cbar_kws={'orientation': "horizontal"})

#cm.cax.set_visible(False)

cm.ax_heatmap.set_xticklabels(cm.ax_heatmap.get_xmajorticklabels(), fontsize = 26)

#cm.ax_heatmap.set_yticklabels(cm.ax_heatmap.get_ymajorticklabels(), fontsize = 16,rotation=40)
#cm.fig.suptitle('SE',fontsize = 36) 

cm.ax_heatmap.set_title('SE',fontsize = 36)

cm.ax_row_dendrogram.legend(handles=legend_handles, loc='lower left', bbox_to_anchor=(0, 1.02),fontsize = 26)

#plt.savefig('/cfs/klemming/home/a/ahsan2/script_pdc/TakeOut/SE_Group_CC_miRNA.pdf', dpi=800)

#SaveImage
plt.savefig(SaveImage, format='pdf', bbox_inches='tight')
#plt.savefig(SaveImage, format='png', bbox_inches='tight')
#plt.savefig(SaveImage, dpi=800)
####################################

#print(GroupAxis)

print('WON')
