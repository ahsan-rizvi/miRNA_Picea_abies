#!/opt/cray/pe/python/3.11.5/bin/python
#- Dr Ahsan Z Rizvi
#-------------------------
import numpy
import sys

#-------------------------

with open(sys.argv[1]) as file:#Read file
    for line in file:
        temp=line.split('\n')
        temp2=temp[0]
        temp3=temp2.split(',')
        #print(temp3)
        #sample_name.append(temp3[0])
        temp4=float(temp3[1])+float(temp3[2])+float(temp3[3])+float(temp3[4])+float(temp3[5])+float(temp3[6])
        temp5=(temp4)/1000000
        #print(temp5)

        a19nt=(float(temp3[1])/temp5)
        a20nt=(float(temp3[2])/temp5)
        a21nt=(float(temp3[3])/temp5)
        a22nt=(float(temp3[4])/temp5)
        a23nt=(float(temp3[5])/temp5)
        a24nt=(float(temp3[6])/temp5)
        #print(int(a19nt), int(a20nt),int(a21nt),int(a22nt), int(a23nt), int(a24nt))
        temp6=temp3[0]+','+str(int(a19nt))+','+str(int(a20nt))+','+str(int(a21nt))+','+str(int(a22nt))+','+str(int(a23nt))+','+str(int(a24nt))
        print(temp6)
