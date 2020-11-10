#!/usr/bin/python
import os
boys=[('PAK_P11_3a', 'PAK_P01_3a'),('LYB_P50_3a', 'PAK_P11_3a'),('VAO_P54_1b', 'LYB_P59_1b'),('AMC_P51_1b', 'VAO_P54_1b'),('VAO_P29_1a', 'LYB_P06_1a'),('VAO_P29_1a', 'AMC_P14_1a'),('LYB_P50_3a', 'PAK_P01_3a'),('AMC_P14_1a', 'LYB_P06_1a'),('LYB_P59_1b', 'AMC_P51_1b')]
for boy in boys:
    tset=[]
    bset=[]
    file=boy[0]+"_givesto_"+boy[1]
    backwards=boy[1]+"_givesto_"+boy[0]
    f=open(file,"r")
    lines=f.readlines()
    f.close()
    for line in lines:
        splitline=line.split("\t")
        tset.append(splitline[1])
    print(file+"\t"+str(len(set(tset))))
    f=open(backwards,"r")
    lines=f.readlines()
    f.close()
    for line in lines:
        splitline=line.split("\t")
        bset.append(splitline[1])
    print(backwards+"\t"+str(len(set(bset))))
    ps="wc -l "+boy[0]+".fas_walker"
    os.system(ps)
    ps="wc -l "+boy[1]+".fas_walker"
    os.system(ps)
