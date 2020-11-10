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

def calcDistanceMatrix(finalSeqs): #calculate distance matrix from the 1-step list
    l=len(finalSeqs)
    arr=np.zeros([l,l])
    hdist=hamming(finalSeqs,finalSeqs,ignore_gaps=False)
    for id in range(len(hdist)):
        item=hdist[id]
        arr[:,id]=item[:,0]
    return arr
    
if __name__=="__main__":
    seqs=getseqs(sys.argv[1])
    seqnum=len(seqs)
    array=calcDistanceMatrix(seqs)
    fixedArray=np.eye(len(array))*100+array #add 100 along the diagonal so that a self distance wont return 0 as mindist
    val=int(np.amin(fixedArray))
    av=np.mean(array)
    print(sys.argv[1]+','+str(seqnum)+','+str(val)+','+str(av)[:5])
