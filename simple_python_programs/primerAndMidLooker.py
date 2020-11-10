#!/usr/bin/python
import sys
import os
import re

fle=sys.argv[1]
if fle[-3:]==".gz":
    asdf="cp "+fle+" asdf.gz"
    os.system(asdf)
    os.system("gunzip asdf.gz")
    fle="asdf"
base=re.findall("(.*)R\d.fastq",fle)
r1=str(base[0])+"R1.fastq"
r2=str(base[0])+"R2.fastq"
l=[r1,r2]
looptrack=0
for fle in l:
    pasdf="Now processing "+fle
    print(pasdf)
    looptrack+=1
    f=open(fle,"r")
    lines=f.readlines()
    f.close()
    counter=0
    types={}
    ptypes={}
    for line in lines:
        counter+=1
        if counter%4==1:
            seqid=line
        elif counter%4==2:
            mid=line[:10]
            if looptrack==1:
                primer=line[10:32]
            elif looptrack==2:
                primer=line[10:30] 
            else:
                print("looptrack broken, fix code")
                break
            if primer not in ptypes:
                ptypes[primer]=1
            else:
                ptypes[primer]+=1
            if mid not in types:
                types[mid]=1
            else:
                types[mid]+=1
    f=open("midList.txt","w")
    for entry in types:
        s=str(types[entry])+"\t"+str(entry)+"\n"
        f.write(s)
    f.close()
    g=open("ListPrimer.txt","w")
    for entry in ptypes:
        astq=str(ptypes[entry])+"\t"+str(entry)+"\n"
        g.write(astq)
    g.close()
    sum=0
    f=open("midList.txt","r")
    lines=f.readlines()
    f.close()
    for line in lines:
        splitline=line.split("\t")
        num=int(splitline[0])
        sum+=num
    os.system("sort -n midList.txt | tail -4")
    os.system("wc -l midList.txt")
    print(sum)
    print("mid info above, primer info below")
    f=open("ListPrimer.txt","r")
    lines=f.readlines()
    sum=0
    for line in lines:
        splitline=line.split("\t")
        num=int(splitline[0])
        sum+=num
    os.system("sort -n ListPrimer.txt | tail -4")
    os.system("wc -l ListPrimer.txt")
    print(sum)
