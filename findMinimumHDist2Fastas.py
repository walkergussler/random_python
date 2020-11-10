print('importing..')
import numpy as np
import sys
from Bio import SeqIO
from pyseqdist import hamming


def getseqs(input): #get sequences from 2 files
    seqs=[]
    input_handle = open(input) 
    for record in SeqIO.parse(input_handle, "fasta"): # for FASTQ use "fastq", for fasta "fasta"
        if len(record.seq) > 0 and len(record.seq) < 50000:
            seqs.append(record.seq)
    input_handle.close()
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
    seqs1=getseqs(sys.argv[1])
    seqs2=getseqs(sys.argv[2])
    print('calculating...')
    array=calcDistanceMatrix(seqs1,seqs2)
    val=np.amin(array)
    print(array)
    print(val)
