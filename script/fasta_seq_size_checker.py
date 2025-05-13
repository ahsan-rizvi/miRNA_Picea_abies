#!/opt/cray/pe/python/3.11.5/bin/python
#- Dr Ahsan Z Rizvi
#-------------------------
import numpy
import sys
#from rnanorm import TMM
#-------------------------

seq_len=list()
with open(sys.argv[1]) as file:#Read file
    for line in file:
        temp=line[0]
        if temp != '>':
           seq_len.append(len(line)-1)
           #print(line,len(line)-1)


nonRed,freq = numpy.unique(seq_len,return_counts=True)
temp=sys.argv[1].split('/') 
temp2=temp[-1]
temp3=temp2.split('.')
#temp4=temp3[0]+','+str(freq[0])+','+str(freq[1])+','+str(freq[2])+','+str(freq[3])+','+str(freq[4])+','+str(freq[5])
temp4=temp2+','+str(freq[0])+','+str(freq[1])+','+str(freq[2])+','+str(freq[3])+','+str(freq[4])+','+str(freq[5])
#print(nonRed)
print(temp4)


