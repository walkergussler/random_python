#!/usr/bin/python
import sys
import os
import re

def genotyper(inf):
    path="/scicomp/groups/OID/NCHHSTP/DVH/xzy3/walker/db.msh"
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
        return False
    res=(min(vallist,key=vallist.get)).split("\t")
    gtype=(res[0][-5:-3])
    print("It looks like your file is genotype: "+gtype)
    return True

genotyper(sys.argv[1])
