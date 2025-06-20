#!/usr/bin/python3
#####!/opt/cray/pe/python/3.11.5/bin/python
#- Dr Ahsan Z Rizvi
#-------------------------
import numpy
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import pandas as pd
import scipy.cluster.hierarchy as sch
#-------------------------
with open(sys.argv[1]) as file:#GroupA reading
    for line in file:
        GroupA=(line.split(','))
print('Total member in Group A:',len(GroupA))
with open(sys.argv[2]) as file:#GroupB reading
    for line in file:
        GroupB=(line.split(','))
print('Total member in Group B:',len(GroupB))
with open(sys.argv[3]) as file:#GroupC reading
    for line in file:
        GroupC=(line.split(','))
print('Total member in Group C:',len(GroupC))
Known=list()
with open(sys.argv[4]) as file:#Known reading
    for line in file:
        Temp1=(line.split('\n'))
        Known.append(Temp1[0])
print('Total member in Known:',len(Known))
resA = list(set(GroupA) & set(Known))
print('Known in GroupA',len(resA))
resB = list(set(GroupB) & set(Known))
print('Known in GroupB',len(resB))
resC = list(set(GroupC) & set(Known))
print('Known in GroupC',len(resC))
Group_all=GroupA+GroupB+GroupC # Concatenate all groups
print('Length in Group',len(Group_all))
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


############################ Read CC file
#miRNA_CC=list()
#miRNA_CC_value=list()
#with open(sys.argv[7]) as file:
 #   for line in file:
  #      temp=line.split('\t')
   #     #print(temp[:])
    #    miRNA_CC.append(temp[0])
     #   miRNA_CC_value.append(temp[3])
#     #   print(temp[0],temp[3])
#miRNA_CC_nonRed,miRNA_CC_freq = numpy.unique(miRNA_CC,return_counts=True)
#miRNA11=list()
#mirBASE11=list()
###########################
GroupB=GroupB[:-1]




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
    #print(temp3)
    f.writelines(temp3)
f.close()
print('miRNA-mRNA,E_value,NCC,P_value file(Total): ',sys.argv[7])
print('(GroupB) miRNA-mRNA,E_value,NCC,P_value file: ',sys.argv[7]+'.GpB.csv')
print('(GroupB) mRNA-miRNA file: ',sys.argv[7]+'.GpB.mRNA-miRNA.tsv')
###################################




####################################### GroupC mRNA
miRNA_CC_GpC=list()#
mRNA_CC_GpC=list()#
temp2=sys.argv[7]+'.GpC.csv'
f = open(temp2, "w")
for i in range(0,len(GroupC)):#Extract GroupA:miRNA,mRNA,GO_ID
    for j in range(0,len(miRNA_CC)):
        if GroupC[i]==miRNA_CC[j]:
           temp=mRNA_CC[j].split('.')
         #  print(miRNA_CC[j],temp[0],CC_value[j],P_value[j])
           temp3=miRNA_CC[j]+','+temp[0]+','+CC_value[j],','+P_value[j]#+'\n'
           miRNA_CC_GpC.append(miRNA_CC[j])#
           mRNA_CC_GpC.append(temp[0])#
           f.writelines(temp3)
f.close()
temp2=sys.argv[7]+'.GpC.mRNA-miRNA.tsv'
f = open(temp2, "w")
mRNA_GpC_nonRed,mRNA_GpC_freq = numpy.unique(mRNA_CC_GpC,return_counts=True)#
for i in range(0,len(mRNA_GpC_nonRed)):#
    temp=''
    for j in range(0,len(mRNA_CC_GpC)):
        if mRNA_GpC_nonRed[i]==mRNA_CC_GpC[j]:
           temp=temp+miRNA_CC_GpC[j]+';'
    #print(mRNA_GpA_nonRed[i],temp)
    temp3=mRNA_GpC_nonRed[i]+'\t'+temp+'\n'
    #print(temp3)
    f.writelines(temp3)
f.close()
print('(GroupC) miRNA-mRNA,E_value,NCC,P_value file: ',sys.argv[7]+'.GpC.csv')
print('(GroupC) mRNA-miRNA file: ',sys.argv[7]+'.GpC.mRNA-miRNA.tsv')
##################################






#sys.exit()



mature_name=list()
mature_seq=list()
mature_total=list()
with open(sys.argv[10]) as file:
    for line in file:
        temp=line.split(',')
        #print(temp)
        if temp[0][0]=='>':
          temp2=temp[0].split(' ')
          temp3=(temp2[0].split('>'))
         # print(temp3[1])
          mature_name.append(temp3[1])
        else:
          temp4=(temp[0].split('\n'))
          #print(temp4[0])
          mature_seq.append(temp4[0])
          


#sys.exit()

#print(miRNA_CC_nonRed)
#print(miRNA_CC_freq)

#sys.exit()

########################### NAT_miRNA

Nat_miRNA=list()
with open(sys.argv[11]) as file:
     for line in file:
         temp=(line.split('\n'))
         #print(temp[0])
         Nat_miRNA.append(temp[0])


########################### Replace CC miRNA with group miRNA 
#miRNA_CC_nonRed=GroupC
miRNA_CC_nonRed=GroupB
#miRNA_CC_nonRed=Nat_miRNA
#PDF_name='/cfs/klemming/home/a/ahsan2/script_pdc/ZE_GroupC.pdf'
PDF_name='/cfs/klemming/home/a/ahsan2/script_pdc/ZE_GroupB.pdf'
#miRNA_CC_nonRed=Group_all
#PDF_name='/cfs/klemming/home/a/ahsan2/script_pdc/TakeOut/ZE_Group_CC_miRNA.pdf'
##########################

print('Heatmap file is generated at: ',PDF_name)







#sys.exit()




########################## Assign miBASE id to -ve CC miRNA 
for i in range(0,len(miRNA_CC_nonRed)): #make list on nonredundant miRNA CC
    temp=miRNA_CC_nonRed[i]#+','+str(miRNA_CC_freq[i])
    for j in range(0,len(miRNA1)):
        if miRNA_CC_nonRed[i]==miRNA1[j]:
           temp=str(mirBASE1[j][0])#+','+str(miRNA_CC_freq[i])
#    print(miRNA_CC_nonRed[i],temp)
    #print(temp)
    miRNA11.append(miRNA_CC_nonRed[i])
    mirBASE11.append(temp)
############################





#sys.exit()





############################### Read VST of ZE
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


###############################
legend_labels = ['GroupA', 'GroupB', 'GroupC', 'Known','Novel']
#legend_colors = ['brown', 'gold', 'green', 'black','purple']
legend_colors = ['#f6a500', '#967ea7', '#ff6289', '#8f6a2e', '#c7c7c7']

legend_handles = [plt.Rectangle((0, 0), 0, 0, color=color, label=label)
                  for color, label in zip(legend_colors, legend_labels)]
###############################


############################### Matrix for plots
matrix3=np.zeros((len(miRNA11),17))# Initiating an matrix for heatmap


y_axis=list()
data=list()
data2=list()
for i in range(0, len(miRNA11)):

    ############################ color setting Known Vs Novel
    temp2='white'
    for m in range(0,len(GroupA)):
        if GroupA[m]==miRNA11[i]:
           #temp2='brown'
           temp2='#f6a500'
    for m in range(0,len(GroupB)):
        if GroupB[m]==miRNA11[i]:
           #temp2='gold'
           temp2='#967ea7'
    #data.append(temp2)
    for m in range(0,len(GroupC)):
        if GroupC[m]==miRNA11[i]:
           #temp2='green'
           temp2='#ff6289'
    data.append(temp2)
    ############################

    ########################### colour setting group
    #temp='purple'
    temp='#c7c7c7'
    for k in range(0,len(Known)):
        if Known[k]==miRNA11[i]:       
           #temp='black'
           temp='#8f6a2e'
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
            S3=numpy.median( [float(temp2[7]),float(temp2[8]),float(temp2[-5]) ] )
            S4=numpy.median([float(temp2[9]),float(temp2[11]),float(temp2[13]) ] )
            S5=numpy.median([float(temp2[15]),float(temp2[17]),float(temp2[19])])
            S6=numpy.median( [float(temp2[22]),float(temp2[24]),float(temp2[-1])])
            S7=numpy.median([float(temp2[26]),float(temp2[28]),float(temp2[-4]) ] )
            S8=numpy.median([float(temp2[31]),float(temp2[33]),float(temp2[35])])
            S9=numpy.median([float(temp2[37]),float(temp2[39])])
            S10=numpy.median([float(temp2[-14]),float(temp2[-11]),float(temp2[-8]),float(temp2[-3]),float(temp2[-2])  ])

            S11=numpy.median( [float(temp2[10]),float(temp2[12]),float(temp2[14]) ] )
            S12=numpy.median( [float(temp2[15]),float(temp2[18]),float(temp2[20]) ] )
            S13=numpy.median( [float(temp2[23]),float(temp2[25]) ] )
            S14=numpy.median([float(temp2[27]),float(temp2[29]),float(temp2[30]) ] )
            S15=numpy.median([float(temp2[32]),float(temp2[34]),float(temp2[36])])
            S16=numpy.median( [float(temp2[38]),float(temp2[40])])
            S17=numpy.median([float(temp2[-6]),float(temp2[-9]),float(temp2[-7]),float(temp2[-10]),float(temp2[-12]),float(temp2[-23]) ] )



            Mean=np.mean([S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S16,S17])
            STD=np.std([S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S16,S17])
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
               matrix3[i][8]=(S9-Mean)/STD
               matrix3[i][9]=(S10-Mean)/STD
                
               matrix3[i][10]=(S11-Mean)/STD
               matrix3[i][11]=(S12-Mean)/STD
               matrix3[i][12]=(S13-Mean)/STD
               matrix3[i][13]=(S14-Mean)/STD
               matrix3[i][14]=(S15-Mean)/STD
               matrix3[i][15]=(S16-Mean)/STD
               matrix3[i][16]=(S17-Mean)/STD
            ################################## Ulrika question: High/low stages
            #Section to add min and max stage
            #if sum( matrix3[i][:])==0:
             #  print(print_index)
                                               
            a=matrix3[i][0:3]
            b=(np.where(a == a.max())   )
            c=(b[0]+1)
            b2=(np.where(a == a.min())   )
            c2=(b2[0]+1)
            min_max=('SD'+str(c)+'/'+'SD'+str(c2))
            #print_index=print_index+','+min_max# OFF it if u dont require min ma
            #---------------------------------
            a=(matrix3[i][3:10])
            b=(np.where(a == a.max())   )
            c=(b[0]+1)
            b2=(np.where(a == a.min())   )
            c2=(b2[0]+1)
            min_max=('FMG'+str(c)+'/'+'FMG'+str(c2))
            #print_index=print_index+','+min_max# OFF it
            #---------------------------------
            a=(matrix3[i][10:17])
            b=(np.where(a == a.max())   )
            c=(b[0]+1)
            b2=(np.where(a == a.min())   )
            c2=(b2[0]+1)
            min_max=('ZE'+str(c)+'/'+'ZE'+str(c2))
            #print_index=print_index+','+min_max# OFF it
    y_axis.append(print_index)

    ################################ Check GroupA, related question
#    if sum( matrix3[i][:])==0 and  matrix3[i][5]==0:#cross check zero VST
#       print('VST file, count file check--------------------------------------------')
#       #print(print_index,sum( matrix3[i][:]))
#       print_index2=print_index.split(',')
#       #Checkin ZE
#       command='grep -w '+print_index2[0]+' '+sys.argv[6]
#       print(command)
#       os.system(command)
#       #Checkin SE
#       #command='grep '+print_index2[0]+' '+sys.argv[8]
#       #print(command)
#       #os.system(command)
#       #Check count table
#       command='grep -w '+print_index2[0]+' '+sys.argv[9]
#       print(command)
#       os.system(command)
#       tempS=(print_index.split(','))
#       tempS2=(tempS)
#       #print(tempS2[0])
#       loc='/cfs/klemming/projects/supr/uppstore2017145/V2/users/ahsan/Ulrika_proj_1/sRNA_dataZE/*remake.fa'
  #     for m in range(0,len(mature_name)):
   #        if tempS2[0]==mature_name[m]:
    #          print(mature_name[m],mature_seq[m])
     #         command='grep -w -B 1 '+mature_seq[m]+' '+loc
      #        #command='grep -w -c '+mature_seq[m]+' '+loc
       #       print(command)
        #      os.system(command)
              

####################################



#sys.exit()







####################################
ser = pd.Series(data=data, index=y_axis)
ser2 = pd.Series(data=data2, index=y_axis)
#GroupAxis=[ser,ser2]

GroupAxis=[ser2]


x_axis=['Z1-SD','Z2-SD','Z3-SD','Z4-FMG','Z5-FMG','Z6-FMG','Z7-FMG','Z8-FMG','Z9-FMG','Z10-FMG','Z4-ZE','Z5-ZE','Z6-ZE','Z7-ZE','Z8-ZE','Z9-ZE','Z10-ZE']
plt.figure()
df = pd.DataFrame(data = matrix3,index = y_axis,columns = x_axis)

#CLUSTROGRAM
colors = ['#0000FF','#FFFFFF' ,"#ff0000"]


custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)

#cm=sns.clustermap(df,col_cluster=False,row_cluster=True, cmap=custom_cmap ,method='ward', metric='euclidean', figsize=(20,40)  ,row_colors=GroupAxis ,dendrogram_ratio=(0.25, 0.15),cbar_pos=[.4, .95, .5, .03],cbar_kws={'orientation': "horizontal"})

cm=sns.clustermap(df,col_cluster=False,row_cluster=True, cmap=custom_cmap ,method='ward', metric='euclidean', yticklabels=False ,cbar_pos=[.4, .95, .5, .03],cbar_kws={'orientation': "horizontal"})

#figsize=(20,40)

#cm.cax.set_visible(False)
cm.cax.set_title("Z-score")
cm.cax.set_visible(True)

cm.ax_row_dendrogram.set_visible(False)
cm.ax_heatmap.set_xticklabels(cm.ax_heatmap.get_xmajorticklabels(),fontsize = 18,rotation=90)
#cm.set_ylabel(fontsize=14)

#cm.ax_heatmap.set_xticklabels(cm.ax_heatmap.get_xmajorticklabels(), fontsize = 10)
#cm.fig.suptitle('SE',fontsize = 36) 
cm.ax_heatmap.set_title('ZE',fontsize = 36)

#cm.ax_row_dendrogram.legend(handles=legend_handles, loc='lower left', bbox_to_anchor=(0, 1.02),fontsize = 9)

plt.savefig(PDF_name, format='pdf', bbox_inches='tight')
#plt.savefig(PDF_name, dpi=800)
#plt.savefig('/cfs/klemming/home/a/ahsan2/script_pdc/TakeOut/ZE_GroupB_CC_miRNA.pdf', dpi=800)
####################################


print('WON')
