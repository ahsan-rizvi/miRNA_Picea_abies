#!/usr/bin/python3.11

from upsetplotly import UpSetPlotly
import matplotlib.pyplot as plt
from upsetplot import generate_counts, plot
from upsetplot import from_contents
from upsetplot import UpSet
import sys
###########################

#with open('./DGE_SE.txt', 'r') as file:
 #   data = file.read().replace('\n', ',')
#print(data)


data1=list()
with open('./Express_SE.txt') as file:#Read file
    for line in file:
        temp1=line.split('\n')
        temp2=temp1[0]
        data1.append(temp2)

data2=list()
with open('./Express_ZE.txt') as file:#Read file
    for line in file:
        temp1=line.split('\n')
        temp2=temp1[0]
        data2.append(temp2)

data3=list()
with open('./Picea_mature_known_id.txt') as file:#Read file
    for line in file:
        temp1=line.split('\n')
        temp2=temp1[0]
        data3.append(temp2)



data4=list()
with open('./Express_SD.txt') as file:#Read file
    for line in file:
        temp1=line.split('\n')
        temp2=temp1[0]
        data4.append(temp2)


data5=list()
with open('./Express_FG.txt') as file:#Read file
    for line in file:
        temp1=line.split('\n')
        temp2=temp1[0]
        data5.append(temp2)


animals = from_contents(
    {"Known": data3, "SE": data1, "ZE": data2,  "SD": data4,  "FG": data5}
)


print(animals)
ax_dict = UpSet(animals, subset_size="count",facecolor="darkblue").plot()
plt.savefig('./xyz.pdf', format='pdf', bbox_inches='tight')


