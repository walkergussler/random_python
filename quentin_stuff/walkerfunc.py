#!/usr/bin/env python
#https://networkx.github.io/documentation/networkx-1.10/reference/index.html

#findTransNetMCMC5
#objTransNetPhyloFitW
#calcLabelsPhyloTree
#simulEvol



def calcMapVal(seqs):#calculate_mapval
    numSeq=len(seqs)
    numPos=len(seqs[0])
    mapval=0;a=0;c=0;g=0;t=0;blank=0
    for pos in range(numPos,0,-1):
        for seq in seqs:
            if seq[pos].lower()=="a":
                a+=1
            elif seq[pos].lower()=="c":
                c+=1
            elif seq[pos].lower()=="t":
                t+=1
            elif seq[pos].lower()=="g":
                g+=1
            else:
                blank+=1
                print("blank")
        consensus=max(a,c,t,g,blank)
        if a+c+t+g+blank<numSeq or consensus>=numSeq or blank>0:
            a=0;c=0;g=0;t=0;blank=0;
        else:
            mapval+=1;a=0;c=0;g=0;t=0;blank=0;
    return mapval
        
def objTNQD2(DM,nhubs):#objTransNetQuadDeg2
    nvert=5
	
def sCore(n,k):#sCore
    d1=(n-k)/k+k-1
    d2=(n-k)/k+1
    return (n-k)/k*d1+(k-1)*(n-k)/k*d2+(k-1)*d1*d2

#treeToTransNet
def testprog():
    for a in range(5):
        print("a",a)
        for b in range(5):
            if b==2:
                break
            else:
                print(b)

def main(): #wrapper
    pass
    
if __name__=="__main__":
    testprog()
    #main()
