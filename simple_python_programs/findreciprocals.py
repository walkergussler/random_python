#!/usr/bin/python
import os
import re

def findrecips():
    tuplist=[]
    reclist=[]
    for f in os.listdir(os.getcwd()):
        if "_givesto_" in f:
            files=re.findall("(.*)_givesto_(.*)",f)
            propertup=(files[0],files[1])        
            revtup=(files[1],files[0])        
            if revtup in tuplist:
                reclist.append(propertup)
            tuplist.append(propertup)
            
