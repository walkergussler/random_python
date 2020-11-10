#!/usr/bin/python
import sys 
import os 
import re 

file=sys.argv[1]
f=open(file,"r") 
lines=f.readlines() 
f.close()
val=0
for line in lines:
    splitline=line.split("\t")
    val+=1
    if float(splitline[2])>.1242:
        print(line)
        break
print(val-2)
