#!/sw/apps/conda/latest/rackham_stage/bin/python3
# Dr Ahsan Z Rizvi
#########################
import numpy
import sys
import os
#import time
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
#########################

list1=list()
seq1=list()
count=0
with open(sys.argv[1]) as file:
    for line in file:
        if line[0] != '>':
           N=len(line)
           #print(N)
           line2=line[0:N-1]
           if (N>17) and ( (N-1)<25 ):
              seq1.append(line2) 
           count=count+1

uniq,counts = numpy.unique(seq1,return_counts=True)

    

#####################################
remake_file_name=sys.argv[1]+'.remake.fa'




i=0
with open(remake_file_name, 'w') as f:
    for x in uniq:
        string='>'+'read'+str("{:08d}".format(i)  )+'_x'+str(counts[i])+'\n'
        f.write(string)
        string2=str(uniq[i])+'\n'
        f.write(string2)
        i=i+1


