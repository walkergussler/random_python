#!/usr/bin/python
import sys
import os
import re
import operator

# just return mid/primers properly, give warning if bad separation, also return file names

#startLetter=sys.argv[1]
startLetter="B"
baselist=[]
for fle in os.listdir(os.getcwd()):
#    if fle[-3:]==".gz":
#        asdf="cp "+fle+" asdf.gz"
#        os.system(asdf)
#        os.system("gunzip asdf.gz")
#        fle="asdf"
    if fle[-6:]!=".fastq" or fle[:1]!=startLetter:
        continue
    base=re.findall("(.*)R\d.fastq",fle)
    if str(base[0]) not in baselist:
        baselist.append(str(base[0]))
    
for base in baselist:
    r1=base+"R1.fastq"
    r2=base+"R2.fastq"
    l=[r1,r2]
    looptrack=0
    for fle in l:
        pasdf="echo Now processing "+fle
        os.system(pasdf)
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
                if primer not in ptypes:
                    ptypes[primer]=1
                else:
                    ptypes[primer]+=1
                if mid not in types:
                    types[mid]=1
                else:
                    types[mid]+=1
        temp_mid = sorted(types.items(), key=operator.itemgetter(1))
        temp_pri =sorted(ptypes.items(), key=operator.itemgetter(1))
        MIDS=temp_mid.reverse()
        PRIM=temp_pri.reverse()
        print(MIDS[1],MIDS[2],PRIM[1],PRIM[2])
        
















##        f=open("midList.txt","w")
##        for entry in types:
##            s=str(types[entry])+"\t"+str(entry)+"\n"
##            f.write(s)
##        f.close()
##        g=open("ListPrimer.txt","w")
##        for entry in ptypes:
##            astq=str(ptypes[entry])+"\t"+str(entry)+"\n"
##            g.write(astq)
##        g.close()
##        sum=0
##        f=open("midList.txt","r")
##        lines=f.readlines()
##        f.close()
##        for line in lines:
##            splitline=line.split("\t")
##            num=int(splitline[0])
##            sum+=num
##        os.system("sort -n midList.txt | tail -4")
##        os.system("wc -l midList.txt")
##        safa="echo "+str(sum)
##        os.system(safa)
##        os.system("echo mid info above, primer info below")
##        f=open("ListPrimer.txt","r")
##        lines=f.readlines()
##        f.close()
##        sum=0
##        for line in lines:
##            splitline=line.split("\t")
##            num=int(splitline[0])
##            sum+=num
##        os.system("sort -n ListPrimer.txt | tail -4")
##        os.system("wc -l ListPrimer.txt")
##        safew="echo "+str(sum)
##        os.system(safew)
