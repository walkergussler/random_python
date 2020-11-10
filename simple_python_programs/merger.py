#!/usr/bin/python
import os
out=[0]*100
c=out
def prog(fil):
    b=[0]*100
    f=open(fil,"r")
    lines=f.readlines()
    f.close()
    counter=0
    for line in lines:
        splitline=line.split("\t")
        val=float(splitline[1])
        b[counter]+=val
        counter+=1
    return b    

for qwer in os.listdir(os.getcwd()):
    if "givesto" in qwer:
        c=prog(qwer)
        for v in range(100):
            out[v]+=c[v]
for entry in out:
    print(entry)
print(len(out))
