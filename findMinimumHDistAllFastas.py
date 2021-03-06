import numpy as np
import sys, os
from Bio import SeqIO
from pyseqdist import hamming
from itertools import combinations

def getseqs(input): #get sequences from 2 files
    seqs=[]
    with open(input) as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
            seqs.append(record.seq)
    return seqs

def calcDistanceMatrix(seqs1,seqs2): #calculate distance matrix from the 1-step list
    hdist=hamming(seqs1,seqs2,ignore_gaps=False)
    l=len(seqs1)
    w=len(seqs2)
    arr=np.zeros([l,w])
    for id in range(len(hdist)):
        item=hdist[id]
        arr[:,id]=item[:,0]
    return arr
    
if __name__=="__main__":
    # files=[]
    # for file in os.listdir(os.getcwd()):
        # if file.endswith(".fas"):
            # files.append(file)
    # for f1,f2 in combinations(files,2):
        # short1=os.path.splitext(os.path.basename(f1))[0]
        # short2=os.path.splitext(os.path.basename(f2))[0]
        # seqs1=getseqs(f1)
        # seqs2=getseqs(f2)
        # array=calcDistanceMatrix(seqs1,seqs2)
        # val=np.amin(array)
        # print(short1+","+short2+","+str(val))
    files={}
    for file in os.listdir(os.getcwd()):
        if file.endswith(".fas"):
            files[file]=getseqs(file)
    for f1,f2 in combinations(files,2):
        seqs1=files[f1]        
        seqs2=files[f2]        
        array=calcDistanceMatrix(seqs1,seqs2)
        val=int(np.amin(array))
        print(f1+","+f2+","+str(val))