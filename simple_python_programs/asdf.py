#!/usr/bin/python
import sys
import os
import re
def genotyper(inf):
    path="/scicomp/groups/OID/NCHHSTP/DVH/mnz0/db.msh"
    os.system("mash dist "+path+" "+inf+" > tmp.txt")
    f=open("tmp.txt","r")
    lines=f.readlines()
    f.close()
    vallist={}
    for line in lines:
        splitline=line.split("\t")
        vallist[line]=float(splitline[2])
    if len(set(vallist.values()))==1:
        print("i honestly have no idea what your genotype is")
        return 2
    res=(min(vallist,key=vallist.get)).split("\t")
    gtype=(res[0][-5:-3])
    filetype=res[1]
    guess=re.findall(".*_([0-9][a-z]).fas",filetype)[0]
    g=0
    if guess==gtype:
        print("good")
        g+=1
    else:
        print("whoops, i thought "+filetype+" was genotype "+gtype)
        if filetype[-13]=="2" and gtype=="2X":
            g+=1
        else:
            os.system("cat tmp.txt")
    if g==1:
        return 1
    else:
        return 0
    
gd=0
b=0
u=0
for i in os.listdir(os.getcwd()):
    if i[-11:]==".fas_walker":
        a=genotyper(i)
        if a==1:
            gd+=1
        elif a==0:
            b+=1
        elif a==2:
            u+=1
print(gd,b,u)
